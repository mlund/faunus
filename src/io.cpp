#include "io.h"
#include "units.h"
#include "random.h"
#include "group.h"
#include "space.h"
#include <spdlog/spdlog.h>
#include <zstr.hpp>
#include <cereal/archives/binary.hpp>
#include <fstream>
#include <iostream>
#include <range/v3/view.hpp>

namespace Faunus {

std::vector<std::string> IO::loadLinesFromFile(const std::string &filename) {
    if (std::ifstream stream(filename); !stream) {
        throw std::runtime_error("load error: "s + filename);
    } else {
        std::string single_line;
        std::vector<std::string> destination;
        while (std::getline(stream, single_line)) {
            destination.push_back(single_line);
        }
        return destination;
    }
}

/**
 * @param filename Filename
 * @param contents String to write
 * @param mode `std::ios_base::out` (new, default) or `std::ios_base::app` (append)
 */
void IO::writeFile(const std::string &filename, const std::string &contents, std::ios_base::openmode mode) {
    if (std::ofstream stream(filename, mode); stream) {
        stream << contents;
    } else {
        throw std::runtime_error("cannot open file for writing: "s + filename);
    }
}

/**
 * @param strings vector of std::string
 * @param pattern Pattern to search for
 *
 * All lines with matching `pattern` will be remove regardless
 * of where in the string it it found
 */
void IO::strip(std::vector<std::string> &strings, const std::string &pattern) {
    auto unary_predicate = [&](auto &s) { return s.find(pattern) != std::string::npos; };
    strings.erase(std::remove_if(strings.begin(), strings.end(), unary_predicate), strings.end());
}

/**
 * Create an output file stream. If the filename extension is `.gz`, gzip
 * compression is enabled, otherwise a standard `std::ostream` is created.
 *
 * @param filename Name of output file
 * @return pointer to stream; nullptr if it could not be created
 */
std::unique_ptr<std::ostream> IO::openCompressedOutputStream(const std::string &filename) {
    if (filename.substr(filename.find_last_of(".") + 1) == "gz") {
        faunus_logger->trace("enabling gzip compression for {}", filename);
        return std::make_unique<zstr::ofstream>(filename);
    }
    return std::make_unique<std::ofstream>(filename);
}

bool FormatAAM::prefer_charges_from_file = true;

Particle FormatAAM::recordToParticle(const std::string &record) {
    std::stringstream o;
    std::string name;
    o << record;
    o >> name;
    if (auto it = findName(atoms, name); it == atoms.end()) {
        throw std::runtime_error("AAM load error: unknown atom name '" + name + "'.");
    } else {
        Particle particle = *it;
        double radius;
        double molecular_weight;
        int atom_number;
        o >> atom_number >> particle.pos.x() >> particle.pos.y() >> particle.pos.z() >> particle.charge >>
            molecular_weight >> radius;

        // does charge match AtomData?
        if (std::fabs(it->charge - particle.charge) > pc::epsilon_dbl) {
            faunus_logger->warn("charge mismatch on atom {0} {1}: {2} atomlist value", atom_number, name,
                                (prefer_charges_from_file) ? "ignoring" : "using");
            if (not prefer_charges_from_file) {
                particle.charge = it->charge; // let's use atomdata charge
            }
        }
        // does radius match AtomData?
        if (std::fabs(it->sigma - 2 * radius) > pc::epsilon_dbl) {
            faunus_logger->warn("radius mismatch on atom {0} {1}: using atomlist value", atom_number, name);
        }
        return particle;
    }
}
ParticleVector FormatAAM::load(const std::string &filename, bool _keepcharges) {
    prefer_charges_from_file = _keepcharges;
    if (auto lines = IO::loadLinesFromFile(filename); lines.empty()) {
        throw std::runtime_error("empty aam file: "s + filename);
    } else {
        IO::strip(lines, "#"); // all lines beginning w. # are comments
        size_t number_of_atoms = std::stoul(lines[0]);
        if (lines.size() < number_of_atoms + 1) {
            throw std::runtime_error("mismatch in number of atoms in "s + filename);
        }
        ParticleVector positions; // positions are loaded here
        positions.reserve(number_of_atoms);
        auto lines_atoms = lines | ranges::views::slice(1uL, 1uL + number_of_atoms);
        std::transform(lines_atoms.begin(), lines_atoms.end(), std::back_inserter(positions),
                       [](auto &i) { return recordToParticle(i); });
        assert(positions.size() == number_of_atoms);
        return positions;
    }
}

void FormatAAM::save(const std::string &file, const ParticleVector &particles) {
    std::ostringstream o;
    o << particles.size() << "\n";
    for (size_t i = 0; i < particles.size(); i++) {
        o << p2s(particles[i], i);
    }
    IO::writeFile(file, o.str());
}
std::string FormatAAM::p2s(const Particle &particle, int zero_based_index) {
    std::ostringstream o;
    o.precision(5);
    auto &atom = Faunus::atoms.at(particle.id);
    o << atom.name << " " << zero_based_index + 1 << " " << particle.pos.transpose() << " " << particle.charge << " "
      << atom.mw << " " << atom.sigma / 2 << "\n";
    return o.str();
}

std::string FormatPQR::writeCryst1(const Point &box_length, const Point &angle) {
    return fmt::format("CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n", box_length.x(),
                       box_length.y(), box_length.z(), angle.x(), angle.y(), angle.z());
}

bool FormatPQR::readCrystalRecord(const std::string &record, Point &box_length) {
    std::stringstream o(record);
    std::string key;
    o >> key;
    if (key == "CRYST1") {
        o >> box_length.x() >> box_length.y() >> box_length.z();
        return true;
    } else {
        return false;
    }
}

bool FormatPQR::readAtomRecord(const std::string &record, Particle &particle, double &radius) {
    std::stringstream o(record);
    std::string key;
    o >> key;
    if (key == "ATOM" or key == "HETATM") {
        int atom_index;
        int res_index;
        std::string atom_name;
        std::string res_name;
        o >> atom_index >> atom_name;
        if (auto it = findName(Faunus::atoms, atom_name); it != Faunus::atoms.end()) {
            particle = *it;
            o >> res_name >> res_index >> particle.pos.x() >> particle.pos.y() >> particle.pos.z() >> particle.charge >>
                radius;
            return true;
        } else {
            throw std::runtime_error("PQR load error: unknown atom name '" + atom_name + "'.");
        }
    }
    return false;
}

/**
 * This will load coordinates from a PQR stream into
 * a particle vector. Atoms are identified by the
 * `ATOM` and `HETATM` keywords and the PQR atom
 * name is used to identify the particle from `AtomMap`.
 * If the unit cell dimension, given by `CRYST1`, is
 * found a non-zero vector is returned.
 *
 * @param stream PQR (file)stream to load from
 * @param destination Target particle vector to *append* to
 * @return Vector with unit cell dimensions. Zero length if not found.
 * @note This is a lazy and unsafe loader that doesn't properly
 *       respect the PQR fixed column format. Has been tested to
 *       work with files from VMD, pdb2pqr, and Faunus.
 */
Point FormatPQR::load(std::istream &stream, ParticleVector &destination, bool keep_charges) {
    Point box_length(0, 0, 0);
    if (stream) {
        bool log_charge = true;
        bool log_radius = true; // emit warnings only once
        std::string record;
        int atom_index = 0;
        while (std::getline(stream, record)) {
            double radius;
            Particle particle;
            if (readAtomRecord(record, particle, radius)) { // read ATOM or HETATOM record
                atom_index++;
                AtomData &atom_data = Faunus::atoms.at(particle.id);
                // does charge match AtomData?
                if (std::fabs(atom_data.charge - particle.charge) > pc::epsilon_dbl) {
                    if (log_charge) {
                        faunus_logger->warn("charge mismatches in PQR file: using {0} values",
                                            keep_charges ? "PQR" : "atomlist");
                        log_charge = false;
                    }
                    if (keep_charges) {
                        faunus_logger->debug("charge mismatch on atom {0} {1}: using {2} (PQR) over {3} (atomlist)",
                                             atom_data.name, atom_index, particle.charge, atom_data.charge);
                    } else {
                        faunus_logger->debug("charge mismatch on atom {0} {1}: using {2} (atomlist) over {3} (PQR)",
                                             atom_data.name, atom_index, atom_data.charge, particle.charge);
                        particle.charge = atom_data.charge;
                    }
                }

                // does radius match AtomData?
                if (std::fabs(atom_data.sigma - 2 * radius) > pc::epsilon_dbl) {
                    if (log_radius) {
                        faunus_logger->warn("radius mismatches in PQR file: using atomlist values");
                        log_radius = false;
                    }
                    faunus_logger->info("radius mismatch on atom {0} {1}: using {2} (atomlist) over {3} (PQR)",
                                        atom_data.name, atom_index, atom_data.sigma / 2, radius);
                }
                destination.push_back(particle);
            } else { // is it a CRYST1 record? (unit cell information)
                readCrystalRecord(record, box_length);
            }
        }
    }
    return box_length;
}

/**
 * @param filename PQR filename to load from
 * @param particles Target particle vector to *append* to
 * @return Vector with unit cell dimensions. Zero length if not found.
 */
Point FormatPQR::load(const std::string &filename, ParticleVector &particles, bool keep_charges) {
    if (std::ifstream stream(filename); !stream) {
        throw std::runtime_error("load error: "s + filename);
    } else {
        return load(stream, particles, keep_charges);
    }
}

/**
 * @param filename PQR trajectory filename
 * @param traj Destination vector of particle vector pointers
 *
 * Each trajectory must be separated by an "END" record. Blank lines are
 * ignored.
 */
void FormatPQR::loadTrajectory(const std::string &filename, std::vector<ParticleVector> &traj) {
    if (std::ifstream stream(filename); !stream) {
        throw std::runtime_error("load error: "s + filename);
    } else {
        traj.clear();
        traj.resize(1); // prepare first frame
        std::string record;
        while (std::getline(stream, record)) {
            Particle particle;
            double radius = 0.0;
            if (readAtomRecord(record, particle, radius)) {
                traj.back().push_back(particle);
            } else if (record.find("END") == 0) {         // if END record, advance to next frame
                traj.back().shrink_to_fit();              // attempt to clean up
                traj.push_back(ParticleVector());         // prepare next frame
                traj.back().reserve(traj.front().size()); // reserve memory
            }
        }
        if (traj.back().empty()) { // delete empty frame after last "END" record
            traj.pop_back();
        }
        if (traj.empty()) {
            faunus_logger->warn("pqr trajectory {} is empty", filename);
        }
    }
}

/**
 * @param stream Output stream
 * @param particles Particle vector
 * @param box_length Unit cell dimensions (default: [0,0,0], not printed)
 * @param atoms_per_residue Number of atoms in each residue (default: automatic)
 */
void FormatPQR::save(std::ostream &stream, const ParticleVector &particles, const Point &box_length,
                     int atoms_per_residue) {
    if (stream and !particles.empty()) {
        if (box_length.norm() > pc::epsilon_dbl) {
            stream << writeCryst1(box_length);
        }
        int residue_cnt = 1;
        int atom_cnt = 1;
        for (const auto &particle : particles) {
            AtomData &d = atoms.at(particle.id);        // atom properties as defined in topology
            auto pos = particle.pos + 0.5 * box_length; // move origin to corner of box (faunus used the middle)
            stream << fmt::format("ATOM  {:5d} {:4} {:4}{:5d}    {:8.3f} {:8.3f} {:8.3f} {:.3f} {:.3f}\n", atom_cnt++,
                                  d.name, d.name, residue_cnt, pos.x(), pos.y(), pos.z(), particle.charge,
                                  0.5 * d.sigma);
            if (d.name == "CTR" or d.name == "HCTR") {
                residue_cnt++;
            } else if (atom_cnt % atoms_per_residue == 0) {
                residue_cnt++;
            }
        }
        stream << "END\n";
    }
}

/**
 * @brief Write vector of groups to output stream
 * @param stream Output stream
 * @param groups Vector of groups
 * @param box_dimensions Box dimensions
 *
 * The residue number follows the group index and inactive particles will
 * have zero charge; zero radius; and positioned in the corner of the box.
 */
void FormatPQR::save(std::ostream &stream, const Tgroup_vector &groups, const Point &box_dimensions) {
    if (stream) {
        if (box_dimensions.norm() > pc::epsilon_dbl) {
            stream << writeCryst1(box_dimensions);
        }
        int residue_cnt = 1;
        int atom_cnt = 1;
        for (const auto &group : groups) {                                                 // loop over all molecules
            for (auto particle = group.begin(); particle != group.trueend(); particle++) { // loop over particles
                double scale = (particle < group.end()) ? 1.0 : 0.0;                       // zero if inactive particle
                auto pos = scale * (particle->pos + 0.5 * box_dimensions);                 // origin to corner of box
                stream << fmt::format("ATOM  {:5d} {:4} {:4}{:5d}    {:8.3f} {:8.3f} {:8.3f} {:.3f} {:.3f}\n", atom_cnt,
                                      particle->traits().name, particle->traits().name, residue_cnt, pos.x(), pos.y(),
                                      pos.z(), scale * particle->charge, scale * particle->traits().sigma / 2.0);
                atom_cnt++;
            }
            residue_cnt++;
        }
        stream << "END\n";
    }
}

/**
 * @param filename Output PQR filename
 * @param particles Particle vector
 * @param box_length Unit cell dimensions (default: [0,0,0], not printed)
 * @param atoms_per_residue Number of atoms in each residue (default: automatic)
 * @exception std::runtime_error if file cannot be opened.
 */
void FormatPQR::save(const std::string &filename, const ParticleVector &particles, const Point &box_length,
                     int atoms_per_residue) {
    if (!particles.empty()) {
        if (std::ofstream stream(filename); stream) {
            save(stream, particles, box_length, atoms_per_residue);
        } else {
            throw std::runtime_error("output error: "s + filename);
        }
    }
}

/**
 * @param filename Output PQR filename
 * @param particles Particle vector
 * @param box_length Unit cell dimensions (default: [0,0,0], not printed)
 * @param n Number of atoms in each residue (default: automatic)
 * @exception std::runtime_error if file cannot be opened.
 */
void FormatPQR::save(const std::string &filename, const Tgroup_vector &groups, const Point &box_length) {
    if (not groups.empty()) {
        if (std::ofstream stream(filename); stream) {
            save(stream, groups, box_length);
        } else {
            throw std::runtime_error("write error: "s + filename);
        }
    }
}

void FormatXYZ::save(const std::string &file, const ParticleVector &particles, const Point &box) {
    std::ostringstream o;
    o << particles.size() << "\n" << box.transpose() << "\n";
    for (const auto &particle : particles) {
        o << atoms.at(particle.id).name << " " << particle.pos.transpose() << "\n";
    }
    IO::writeFile(file, o.str());
}

TEST_CASE("[Faunus] FormatPQR") {
    using doctest::Approx;
    // Space object with two salt pairs, i.e. four particles
    Space spc;
    SpaceFactory::makeNaCl(spc, 2, R"( {"type": "cuboid", "length": [20,30,40]} )"_json);
    double d = 0;

    // fix positions
    for (int i = 0; i < 4; i++) {
        spc.p[i].pos = {d, d + 0.1, d + 0.2};
        spc.p[i].charge = double(i);
        d += 0.5;
    }

    // write PQR stream
    std::ostringstream out;
    FormatPQR::save(out, spc.p, spc.geo.getLength());

    // read from stream
    std::stringstream in(out.str());
    spc.p.clear();
    Point box_length = FormatPQR::load(in, spc.p, true);
    CHECK(box_length.x() == 20);
    CHECK(box_length.y() == 30);
    CHECK(box_length.z() == 40);
    CHECK(spc.p.size() == 4);

    // check if positions are restored
    d = 0;
    for (int i = 0; i < 4; i++) {
        spc.p[i].pos = spc.p[i].pos - box_length / 2;
        CHECK(spc.p[i].pos.x() == Approx(d));
        CHECK(spc.p[i].pos.y() == Approx(d + 0.1));
        CHECK(spc.p[i].pos.z() == Approx(d + 0.2));
        CHECK(spc.p[i].charge == double(i));
        d += 0.5;
    }
}

/**
 * Loads coordinates from XYZ file into particle vector. Remaining atom properties
 * are taken from the defined atom list; error if the atom name is unknown.
 *
 * @param filename Filename
 * @param particles Destination particle vector
 * @param append True means append to `particles` (default). If `false`,`particles` is first cleared.
 */

void FormatXYZ::load(const std::string &filename, ParticleVector &particles, bool append) {
    if (std::ifstream stream(filename); !stream) {
        throw std::runtime_error("cannot open file: "s + filename);
    } else {
        std::string comment;
        std::string atom_name;
        size_t number_of_atoms;
        stream >> number_of_atoms;     // first line is number of atoms
        std::getline(stream, comment); // ">>" token doesn't gobble new line
        std::getline(stream, comment); // read comment line
        if (append == false) {
            particles.clear();
        }
        particles.reserve(particles.size() + number_of_atoms);
        for (size_t i = 0; i < number_of_atoms; i++) {
            stream >> atom_name;
            if (auto it = findName(Faunus::atoms, atom_name); it != Faunus::atoms.end()) {
                Particle particle = *it;
                stream >> particle.pos.x() >> particle.pos.y() >> particle.pos.z();
                particles.push_back(particle);
            } else {
                throw std::runtime_error("XYZ load error: unknown atom name '" + atom_name + "'.");
            }
        }
        assert(particles.size() == number_of_atoms);
    }
}

Particle FormatGRO::recordToParticle(const std::string &record) {
    std::stringstream o;
    std::string atom_name;
    o << record.substr(10, 5) << record.substr(20, 8) << record.substr(28, 8) << record.substr(36, 8);
    o >> atom_name;
    if (auto it = findName(Faunus::atoms, atom_name); it != Faunus::atoms.end()) {
        Particle particle = *it; // copy all atom properties, except positions
        o >> particle.pos.x() >> particle.pos.y() >> particle.pos.z();
        particle.pos *= 1.0_nm; // GRO files use nanometers
        return particle;
    } else {
        throw std::runtime_error("GRO file: unknown atom name: " + atom_name);
    }
}

/**
 * @param filename Filename
 * @returns Particle vector
 *
 * The first line in a GRO file is a comment; the
 * second the number of atoms; and the last the
 * box dimensions. All lengths in nanometers.
 */
ParticleVector FormatGRO::load(const std::string &filename) {
    constexpr size_t header_size = 2; // two lines before positions
    if (auto lines = IO::loadLinesFromFile(filename); lines.size() <= header_size) {
        throw std::runtime_error("GRO file seems empty: " + filename);
    } else {
        size_t number_of_atoms = std::stoul(lines[1]); // second line = number of atoms
        if (lines.size() < header_size + number_of_atoms) {
            throw std::runtime_error("Mismatch in number of atoms in GRO file");
        } else {
            ParticleVector positions; // positions are loaded here
            positions.reserve(number_of_atoms);
            std::transform(lines.begin(), lines.begin() + number_of_atoms, std::back_inserter(positions),
                           [](auto &i) { return recordToParticle(i); });
            return positions;
        }
    }
}

void FormatGRO::save(const std::string &filename, const Space &spc) {
    if (std::ofstream stream(filename); !stream) {
        std::runtime_error("write error: "s + filename);
    } else {
        int number_of_residues = 1;
        int number_of_atoms = 1;
        Point boxlength = spc.geo.getLength();
        stream << "Generated by Faunus -- https://github.com/mlund/faunus\n" << spc.numParticles() << "\n";
        for (const auto &group : spc.groups) { // loop over groups
            for (auto &particle : group) {     // loop over active particles
                auto &atom_name = Faunus::atoms.at(particle.id).name;
                Point pos = (particle.pos + 0.5 * boxlength) / 1.0_nm; // shift origin and convert to nm
                stream << fmt::format("{:5d}{:5}{:5}{:5d}{:8.3f}{:8.3f}{:8.3f}\n", number_of_residues, atom_name,
                                      atom_name, number_of_atoms++, pos.x(), pos.y(), pos.z());
            }
            number_of_residues++;
        }
        stream << boxlength.transpose() / 1.0_nm << "\n";
    }
}

// ========== XTCTrajectoryFrame ==========

XTCTrajectoryFrame::XTCTrajectoryFrame(int number_of_atoms) { initNumberOfAtoms(number_of_atoms); }

XTCTrajectoryFrame::XTCTrajectoryFrame(const TrajectoryFrame &frame) {
    initNumberOfAtoms(frame.coordinates.size());
    importFrame(frame);
}

void XTCTrajectoryFrame::operator=(const TrajectoryFrame &frame) {
    if (frame.coordinates.size() != number_of_atoms) {
        throw std::runtime_error("wrong number of particles to be assign into the XTC frame");
    }
    importFrame(frame);
}

void XTCTrajectoryFrame::importFrame(const TrajectoryFrame &frame) {
    importTimestamp(frame.step, frame.timestamp);
    importBox(frame.box);
    importCoordinates(frame.coordinates, 0.5 * frame.box);
}

void XTCTrajectoryFrame::importTimestamp(const int step, const float time) {
    xtc_step = step;
    xtc_time = time / 1.0_ps;
}

void XTCTrajectoryFrame::importBox(const Point &box) {
    // empty box tensor
    XTCMatrix xtc_box_matrix = XTCMatrix::Zero();
    // only XYZ dimensions in nanometers on diagonal, as floats
    xtc_box_matrix.diagonal() = (box / 1.0_nm).cast<XTCFloat>();
    // copy underlaying eigen structure (1D array, row-major) to the C-style 2D array
    std::copy(xtc_box_matrix.data(), xtc_box_matrix.data() + DIM * DIM, &(xtc_box[0][0]));
}

void XTCTrajectoryFrame::importCoordinates(const PointVector &coordinates, const Point &offset) {
    // setNumberOfAtoms(coordinates.size());
    if (coordinates.size() != number_of_atoms) {
        // to avoid mistakes, the number_of_atoms is immutable
        throw std::runtime_error("wrong number of particles to be saved in the XTC frame");
    }
    for (int i = 0; i < number_of_atoms; ++i) {
        // coordinates shifted by an offset (i.e., box / 2), in nanometers, as floats
        const XTCVector xtc_pos = ((coordinates[i] + offset) / 1.0_nm).cast<XTCFloat>();
        // copy underlaying eigen structure (1D array) to the correct place in C-style 2D array
        std::copy(xtc_pos.data(), xtc_pos.data() + DIM, xtc_coordinates.get()[i]);
    }
}

void XTCTrajectoryFrame::exportFrame(TrajectoryFrame &frame) const {
    exportTimestamp(frame.step, frame.timestamp);
    exportBox(frame.box);
    exportCoordinates(frame.coordinates, 0.5 * frame.box);
}

void XTCTrajectoryFrame::exportTimestamp(int &step, float &time) const {
    step = xtc_step;
    time = xtc_time * 1.0_ps;
}

void XTCTrajectoryFrame::exportBox(Point &box) const {
    XTCMatrix xtc_box_matrix = Eigen::Map<const XTCTrajectoryFrame::XTCMatrix>(&(xtc_box[0][0]));
    if (xtc_box_matrix.diagonal().asDiagonal().toDenseMatrix() != xtc_box_matrix) {
        throw std::runtime_error("cannot load non-orthogonal box");
    }
    box = Point(xtc_box_matrix.diagonal().cast<double>() * 1.0_nm);
}

void XTCTrajectoryFrame::exportCoordinates(PointVector &coordinates, const Point &offset) const {
    if (coordinates.size() != number_of_atoms) {
        throw std::runtime_error("wrong number of particles in the loaded XTC frame");
    }
    for (size_t i = 0; i < number_of_atoms; ++i) {
        XTCVector xtc_atom_coordinates(xtc_coordinates.get()[i]);
        coordinates[i] = Point(xtc_atom_coordinates.cast<double>() * 1.0_nm) - offset;
    }
}

void XTCTrajectoryFrame::initNumberOfAtoms(int new_number_of_atoms) {
    assert(new_number_of_atoms >= 0);
    if (number_of_atoms != new_number_of_atoms) {
        number_of_atoms = new_number_of_atoms;
        xtc_coordinates = std::unique_ptr<XDRfile::rvec[]>(new XDRfile::rvec[number_of_atoms]);
    }
}

// ========== TrajectoryFrame ==========

TrajectoryFrame::TrajectoryFrame(const Point &box, const PointVector &coordinates, int step, float timestamp)
    : box(box), coordinates(coordinates), step(step), timestamp(timestamp) {}

TrajectoryFrame::TrajectoryFrame(const XTCTrajectoryFrame &xtc_frame) {
    coordinates.resize(xtc_frame.number_of_atoms);
    xtc_frame.exportFrame(*this);
}

void TrajectoryFrame::operator=(const XTCTrajectoryFrame &xtc_frame) { xtc_frame.exportFrame(*this); }

TEST_CASE("XTCFrame") {
    using doctest::Approx;
    TrajectoryFrame frame({10.0, 12.0, 8.0}, {{1.0, 2.0, -1.0}, {-4.0, 4.0, 2.0}}, 10, 0.2);

    SUBCASE("To XTCFrame") {
        XTCTrajectoryFrame xtc_frame(frame);
        CHECK_EQ(xtc_frame.xtc_step, 10);
        CHECK_EQ(xtc_frame.xtc_time, Approx(0.2 / 1._ps));
        CHECK_EQ(xtc_frame.xtc_box[1][1], Approx(12. / 1._nm));
        CHECK_EQ(xtc_frame.xtc_box[0][1], 0);
        REQUIRE_EQ(xtc_frame.number_of_atoms, 2);
        CHECK_EQ(xtc_frame.xtc_coordinates.get()[1][2], Approx((2.0 + 4.0) / 1._nm));

        SUBCASE("From XTCFrame") {
            TrajectoryFrame new_frame(xtc_frame);
            CHECK_EQ(new_frame.step, frame.step);
            CHECK_EQ(new_frame.timestamp, Approx(frame.timestamp));
            CHECK(new_frame.box.isApprox(frame.box, 1e-6));
            CHECK(new_frame.coordinates[0].isApprox(frame.coordinates[0], 1e-6));
            CHECK(new_frame.coordinates[1].isApprox(frame.coordinates[1], 1e-6));
        }

        SUBCASE("From XTCFrame Iterator") {
            TrajectoryFrame new_frame;
            PointVector coordinates{2};
            REQUIRE_EQ(coordinates.size(), 2);
            xtc_frame.exportFrame(new_frame.step, new_frame.timestamp, new_frame.box, coordinates.begin(),
                                  coordinates.end());
            CHECK(coordinates[0].isApprox(frame.coordinates[0], 1e-6));
            CHECK(coordinates[1].isApprox(frame.coordinates[1], 1e-6));
            PointVector coordinates_too_small{1}, coordinates_too_big{3};
            REQUIRE_LT(coordinates_too_small.size(), 2);
            REQUIRE_GT(coordinates_too_big.size(), 2);
            CHECK_THROWS_AS(xtc_frame.exportFrame(new_frame.step, new_frame.timestamp, new_frame.box,
                                                  coordinates_too_small.begin(), coordinates_too_small.end()),
                            std::runtime_error);
            CHECK_THROWS_AS(xtc_frame.exportFrame(new_frame.step, new_frame.timestamp, new_frame.box,
                                                  coordinates_too_big.begin(), coordinates_too_big.end()),
                            std::runtime_error);
        }
    }
}

// ========== XTCReader ==========

XTCReader::XTCReader(const std::string &filename) : filename(filename) {
    int number_of_atoms;
    if (XDRfile::read_xtc_natoms(filename.c_str(), &number_of_atoms) == XDRfile::exdrOK) {
        xtc_frame = std::make_shared<XTCTrajectoryFrame>(number_of_atoms);
        xdrfile = XDRfile::xdrfile_open(filename.c_str(), "r");
    }
    if (!xtc_frame || !xdrfile) {
        throw std::runtime_error(fmt::format("xtc file {} could not be opened", filename));
    }
}

XTCReader::~XTCReader() { XDRfile::xdrfile_close(xdrfile); }

int XTCReader::getNumberOfCoordinates() { return xtc_frame->number_of_atoms; }

bool XTCReader::readFrame() {
    return_code = XDRfile::read_xtc(xdrfile, xtc_frame->number_of_atoms, &xtc_frame->xtc_step, &xtc_frame->xtc_time,
                                    xtc_frame->xtc_box, xtc_frame->xtc_coordinates.get(), &xtc_frame->precision);
    if (return_code != XDRfile::exdrENDOFFILE && return_code != XDRfile::exdrOK) {
        throw std::runtime_error(fmt::format("xtc file {} could not be read (error code {})", filename, return_code));
    }
    return return_code == XDRfile::exdrOK;
}

bool XTCReader::read(TrajectoryFrame &frame) {
    bool is_ok = readFrame();
    if (is_ok) {
        frame = *xtc_frame;
    }
    return is_ok;
}

// ========== XTCWriter ==========

XTCWriter::XTCWriter(const std::string &filename) : filename(filename) {
    xdrfile = XDRfile::xdrfile_open(filename.c_str(), "w");
    if (!xdrfile) {
        throw std::runtime_error(fmt::format("xtc file {} could not be opened", filename));
    }
}

XTCWriter::~XTCWriter() { XDRfile::xdrfile_close(xdrfile); }

void XTCWriter::writeFrameAt(int step, float time) {
    return_code = XDRfile::write_xtc(xdrfile, xtc_frame->number_of_atoms, step, time, xtc_frame->xtc_box,
                                     xtc_frame->xtc_coordinates.get(), xtc_frame->precision);
    if (return_code != XDRfile::exdrOK) {
        throw std::runtime_error(
            fmt::format("xtc file {} could not be written (error code {})", filename, return_code));
    }
}

void XTCWriter::writeFrame() {
    return_code = XDRfile::write_xtc(xdrfile, xtc_frame->number_of_atoms, xtc_frame->xtc_step, xtc_frame->xtc_time,
                                     xtc_frame->xtc_box, xtc_frame->xtc_coordinates.get(), xtc_frame->precision);
    if (return_code != XDRfile::exdrOK) {
        throw std::runtime_error(
            fmt::format("xtc file {} could not be written (error code {})", filename, return_code));
    }
}

void XTCWriter::write(const TrajectoryFrame &frame) {
    if (!xtc_frame) {
        xtc_frame = std::make_shared<XTCTrajectoryFrame>(frame.coordinates.size());
    }
    *xtc_frame = frame;
    writeFrame();
    step_counter = frame.step + 1;
}

void XTCWriter::writeNext(const TrajectoryFrame &frame) {
    if (!xtc_frame) {
        xtc_frame = std::make_shared<XTCTrajectoryFrame>(frame.coordinates.size());
    }
    *xtc_frame = frame;
    writeFrameAt(step_counter, step_counter * time_delta);
    ++step_counter;
}

ParticleVector fastaToParticles(const std::string &fasta_sequence, double bond_length, const Point &origin) {
    ParticleVector particles;                  // particle vector
    auto ids = fastaToAtomIds(fasta_sequence); // convert letters to atom ids
    std::transform(ids.begin(), ids.end(), std::back_inserter(particles), [&](auto &id) {
        Particle particle = Faunus::atoms.at(id);
        particle.pos = particles.empty() ? origin : particles.back().pos + ranunit(random) * bond_length;
        return particle;
    });
    return particles;
}

ParticleVector loadStructure(const std::string &file, bool keep_charges) {
    try {
        ParticleVector particles;
        if (std::string suffix = file.substr(file.find_last_of(".") + 1); suffix == "aam") {
            particles = FormatAAM::load(file, keep_charges);
        } else if (suffix == "pqr") {
            FormatPQR::load(file, particles, keep_charges);
        } else if (suffix == "xyz") {
            FormatXYZ::load(file, particles);
        } else if (suffix == "gro") {
            particles = FormatGRO::load(file);
        } else {
            throw std::runtime_error("unknown format");
        }
        if (particles.empty()) {
            throw std::runtime_error("empty structure");
        }
        return particles;
    } catch (std::exception &e) {
        throw std::runtime_error("structure load error: "s + e.what());
    }
}

/**
 * @param fasta_sequence FASTA sequence, capital letters.
 * @return vector of verified and existing atom id's
 */
std::vector<int> fastaToAtomIds(const std::string &fasta_sequence) {
    const std::map<char, std::string> map = {{'A', "ALA"},
                                             {'R', "ARG"},
                                             {'N', "ASN"},
                                             {'D', "ASP"},
                                             {'C', "CYS"},
                                             {'E', "GLU"},
                                             {'Q', "GLN"},
                                             {'G', "GLY"},
                                             {'H', "HIS"},
                                             {'I', "ILE"},
                                             {'L', "LEU"},
                                             {'K', "LYS"},
                                             {'M', "MET"},
                                             {'F', "PHE"},
                                             {'P', "PRO"},
                                             {'S', "SER"},
                                             {'T', "THR"},
                                             {'W', "TRP"},
                                             {'Y', "TYR"},
                                             {'V', "VAL"},
                                             // faunus specific codes
                                             {'n', "NTR"},
                                             {'c', "CTR"},
                                             {'a', "ANK"}};

    std::vector<std::string> names;
    names.reserve(fasta_sequence.size());

    for (auto letter : fasta_sequence) {                   // loop over letters
        if (auto it = map.find(letter); it != map.end()) { // is it in map?
            names.push_back(it->second);
        } else {
            throw std::runtime_error("invalid FASTA letter '" + std::string(1, letter) + "'");
        }
    }
    return Faunus::names2ids(atoms, names);
}

FormatSpaceTrajectory::FormatSpaceTrajectory(std::ostream &ostream) {
    if (ostream) {
        output_archive = std::make_unique<cereal::BinaryOutputArchive>(ostream);
    }
}
FormatSpaceTrajectory::FormatSpaceTrajectory(std::istream &istream) {
    if (istream) {
        input_archive = std::make_unique<cereal::BinaryInputArchive>(istream);
    }
}
void FormatSpaceTrajectory::load(Space &) { assert(input_archive != nullptr); }
void FormatSpaceTrajectory::save(const Space &) { assert(output_archive != nullptr); }
} // namespace Faunus
