#include "io.h"
#include "units.h"
#include "random.h"
#include "spdlog/spdlog.h"
#include <fstream>
#include <iostream>

namespace Faunus {

bool IO::readFile(const std::string &file, std::vector<std::string> &v) {
    std::ifstream f(file);
    if (f) {
        std::string s;
        while (getline(f, s))
            v.push_back(s);
        f.close();
        return true;
    }
    faunus_logger->warn("cannot read file: {0}", file);
    return false;
}

/**
 * @param file Filename
 * @param s String to write
 * @param mode `std::ios_base::out` (new, default) or `std::ios_base::app` (append)
 */
bool IO::writeFile(const std::string &file, const std::string &s, std::ios_base::openmode mode) {
    std::ofstream f(file, mode);
    if (f) {
        f << s;
        return true;
    }
    return false;
}

/**
 * @param string_vector vector of std::string
 * @param pattern Pattern to search for
 */
void IO::strip(std::vector<std::string> &string_vector, const std::string &pattern) {
    auto iter = string_vector.begin();
    while (iter < string_vector.end())
        if ((*iter).find(pattern) != std::string::npos)
            string_vector.erase(iter);
        else
            ++iter;
}

int FormatXTC::getNumAtoms() { return natoms_xtc; }

bool FormatAAM::keepcharges = true;

Particle &FormatAAM::s2p(const std::string &record, Particle &particle) {
    std::stringstream o;
    std::string name;
    double radius, mw;
    int num;
    o << record;
    o >> name;
    auto it = findName(atoms, name);
    if (it == atoms.end())
        throw std::runtime_error("AAM load error: unknown atom name '" + name + "'.");

    particle = *it;
    o >> num >> particle.pos.x() >> particle.pos.y() >> particle.pos.z() >> particle.charge >> mw >> radius;

    // does charge match AtomData?
    if (std::fabs(it->charge - particle.charge) > pc::epsilon_dbl) {
        faunus_logger->warn("charge mismatch on atom {0} {1}: {2} atomlist value", num, name,
                            (keepcharges) ? "ignoring" : "using");
        if (not keepcharges)
            particle.charge = it->charge; // let's use atomdata charge
    }

    // does radius match AtomData?
    if (std::fabs(it->sigma - 2 * radius) > pc::epsilon_dbl)
        faunus_logger->warn("radius mismatch on atom {0} {1}: using atomlist value", num, name);
    return particle;
}
bool FormatAAM::load(const std::string &filename, ParticleVector &target, bool _keepcharges) {
    keepcharges = _keepcharges;
    std::vector<std::string> v;
    target.clear();
    if (IO::readFile(filename, v)) {
        IO::strip(v, "#");
        unsigned int n = atoi(v[0].c_str());
        target.resize(n);
        for (unsigned int i = 1; i <= n; i++) {
            assert(i < v.size());
            s2p(v.at(i), target.at(i - 1));
        }
        return true;
    }
    return false;
}
bool FormatAAM::save(const std::string &file, const ParticleVector &particle_vector) {
    std::ostringstream o;
    o << particle_vector.size() << std::endl;
    for (size_t i = 0; i < particle_vector.size(); i++)
        o << p2s(particle_vector[i], i);
    return IO::writeFile(file, o.str());
}
std::string FormatAAM::p2s(const Particle &particle, int zero_based_index) {
    std::ostringstream o;
    o.precision(5);
    auto &atom = Faunus::atoms.at(particle.id);
    o << atom.name << " " << zero_based_index + 1 << " " << particle.pos.transpose() << " " << particle.charge << " "
      << atom.mw << " " << atom.sigma / 2 << std::endl;
    return o.str();
}

bool FormatXTC::open(std::string filename) {
    if (xd != NULL)
        close();
    xd = xdrfile_open(&filename[0], "r");
    if (xd != NULL) {
        int rc = read_xtc_natoms(&filename[0], &natoms_xtc); // get number of atoms
        if (rc == exdrOK) {
            x_xtc = new rvec[natoms_xtc]; // resize coordinate array
            return true;
        }
    } else
        faunus_logger->warn("xtc file could not be opened");
    return false;
}

void FormatXTC::close() {
    xdrfile_close(xd);
    xd = NULL;
    delete[] x_xtc;
}

FormatXTC::FormatXTC(double len) {
    prec_xtc = 1000.;
    time_xtc = step_xtc = 0;
    setbox(len);
    xd = NULL;
    x_xtc = NULL;
}

FormatXTC::~FormatXTC() { close(); }

void FormatXTC::setbox(double x, double y, double z) {
    assert(x > 0 && y > 0 && z > 0);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            xdbox[i][j] = 0;
    xdbox[0][0] = 0.1 * x; // corners of the
    xdbox[1][1] = 0.1 * y; // rectangular box
    xdbox[2][2] = 0.1 * z; // in nanometers!
}

void FormatXTC::setbox(double box_length) { setbox(box_length, box_length, box_length); }

void FormatXTC::setbox(const Point &box_length) { setbox(box_length.x(), box_length.y(), box_length.z()); }

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
    } else
        return false;
}

bool FormatPQR::readAtomRecord(const std::string &record, Particle &particle, double &radius) {
    std::stringstream o(record);
    std::string key;
    o >> key;
    if (key == "ATOM" or key == "HETATM") {
        int atom_index, res_index;
        std::string atom_name, res_name;
        o >> atom_index >> atom_name;
        auto it = findName(Faunus::atoms, atom_name);
        if (it == atoms.end())
            throw std::runtime_error("PQR load error: unknown atom name '" + atom_name + "'.");
        particle = *it;
        o >> res_name >> res_index >> particle.pos.x() >> particle.pos.y() >> particle.pos.z() >> particle.charge >>
            radius;
        return true;
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
        bool log_charge = true, log_radius = true; // emit warnings only once
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
                    if (keep_charges)
                        faunus_logger->debug("charge mismatch on atom {0} {1}: using {2} (PQR) over {3} (atomlist)",
                                             atom_data.name, atom_index, particle.charge, atom_data.charge);
                    else {
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
            } else // is it a CRYST1 record? (unit cell information)
                readCrystalRecord(record, box_length);
        }
    }
    return box_length;
}

/**
 * @param filename PQR filename to load from
 * @param destination Target particle vector to *append* to
 * @return Vector with unit cell dimensions. Zero length if not found.
 */
Point FormatPQR::load(const std::string &filename, ParticleVector &destination, bool keep_charges) {
    std::ifstream f(filename);
    return load(f, destination, keep_charges);
}

/**
 * @param filename PQR trajectory filename
 * @param traj Destination vector of particle vector pointers
 *
 * Each trajectory must be separated by an "END" record. Blank lines are
 * ignored.
 */
void FormatPQR::loadTrajectory(const std::string &filename, std::vector<ParticleVector> &traj) {
    std::ifstream file(filename);
    if (file) {
        std::string record;
        traj.clear();
        traj.resize(1); // prepare first frame
        while (std::getline(file, record)) {
            Particle particle;
            double radius = 0;
            if (readAtomRecord(record, particle, radius)) {
                traj.back().push_back(particle);
            } else if (record.find("END") == 0) {         // if END record, advance to next frame
                traj.back().shrink_to_fit();              // attempt to clean up
                traj.push_back(ParticleVector());         // prepare next frame
                traj.back().reserve(traj.front().size()); // reserve memory
            }
        }
        if (traj.empty())
            faunus_logger->warn("pqr trajectory {} is empty", filename);

    } else
        throw std::runtime_error("pqr file not found: " + filename);
}

/**
 * @param stream Output stream
 * @param particles Particle vector
 * @param box_length Unit cell dimensions (default: [0,0,0], not printed)
 * @param n Number of atoms in each residue (default: automatic)
 */
bool FormatPQR::save(std::ostream &stream, const ParticleVector &particles, Point box_length, int n) {
    if (stream and !particles.empty()) {
        int residue_cnt = 1, atom_cnt = 1;
        if (box_length.norm() > pc::epsilon_dbl)
            stream << writeCryst1(box_length);
        for (auto &i : particles) {
            AtomData &d = atoms.at(i.id);         // atom properties as defined in topology
            Point pos = i.pos + 0.5 * box_length; // move origin to corner of box (faunus used the middle)
            stream << fmt::format("ATOM  {:5d} {:4} {:4}{:5d}    {:8.3f} {:8.3f} {:8.3f} {:.3f} {:.3f}\n", atom_cnt++,
                                  d.name, d.name, residue_cnt, pos.x(), pos.y(), pos.z(), i.charge, 0.5 * d.sigma);
            if (d.name == "CTR")
                residue_cnt++;
            else if (atom_cnt % n == 0)
                residue_cnt++;
        }
        stream << "END\n";
        return true;
    }
    return false;
}

/**
 * @param filename Output PQR filename
 * @param particles Particle vector
 * @param box_length Unit cell dimensions (default: [0,0,0], not printed)
 * @param n Number of atoms in each residue (default: automatic)
 */
bool FormatPQR::save(const std::string &filename, const ParticleVector &particles, Point box_length, int n) {
    if (not particles.empty()) {
        std::ofstream file(filename);
        return save(file, particles, box_length, n);
    }
    return false;
}

bool FormatXYZ::save(const std::string &file, const ParticleVector &particles, const Point &box_length) {
    std::ostringstream o;
    o << particles.size() << "\n" << box_length.transpose() << "\n";
    for (auto &i : particles)
        o << atoms.at(i.id).name << " " << i.pos.transpose() << "\n";
    return IO::writeFile(file, o.str());
}

/**
 * Loads coordinates from XYZ file into particle vector. Remaining atom properties
 * are taken from the defined atom list; a warning is issued of the atom name
 * is unknown.
 *
 * @param filename Filename
 * @param particle_vector Destination particle vector
 * @param append True means append to `p` (default). If `false`,`p` is first cleared.
 */

bool FormatXYZ::load(const std::string &filename, ParticleVector &particle_vector, bool append) {
    std::ifstream f(filename);
    if (f) {
        if (append == false)
            particle_vector.clear();
        Particle particle;
        std::string comment, name;
        size_t n;
        f >> n;
        particle_vector.reserve(particle_vector.size() + n);
        std::getline(f, comment); // ">>" token doesn't gobble new line
        std::getline(f, comment); // read comment line
        for (size_t i = 0; i < n; i++) {
            f >> name;
            auto it = findName(atoms, name);
            if (it == atoms.end())
                throw std::runtime_error("XYZ load error: unknown atom name '" + name + "'.");
            particle = *it;
            f >> particle.pos.x() >> particle.pos.y() >> particle.pos.z();
            particle_vector.push_back(particle);
        }
        if (!particle_vector.empty())
            return true;
    }
    return false;
}

std::string FormatMXYZ::p2s(const Particle &, int) {
    assert(false && "under construction");
    std::ostringstream o;
    o.precision(5);
    // o << a.transpose() << " " << a.dir.transpose() << " " << a.patchdir.transpose() << "\n";
    return o.str();
}

Particle &FormatMXYZ::s2p(const std::string &s, Particle &a) {
    assert(false && "under construction");
    std::stringstream o;
    o << s;
    // o >> a.x() >> a.y() >> a.z() >> a.dir.x() >> a.dir.y() >> a.dir.z() >> a.patchdir.x() >> a.patchdir.y() >>
    //    a.patchdir.z();
    // a.init();
    return a;
}

bool FormatMXYZ::save(const std::string &file, const ParticleVector &p, const Point &len, int time) {
    std::ostringstream o;
    o << p.size() << "\n"
      << "sweep " << time << "; box " << len.transpose() << "\n";
    for (size_t i = 0; i < p.size(); i++)
        o << p2s(p[i], i);
    return IO::writeFile(file, o.str(), std::ios_base::app);
}

bool FormatMXYZ::load(const std::string &file, ParticleVector &p, Point &len) {
    std::stringstream o;
    std::vector<std::string> v;
    if (IO::readFile(file, v)) {
        IO::strip(v, "#");
        size_t n = atoi(v[0].c_str());
        if (p.size() != n)
            faunus_logger->error("mxyz load error: number of particles in xyz file {0} does not match input file ({1})",
                                 n, p.size());
        o << v[1].erase(0, v[1].find_last_of("x") + 1);
        o >> len.x() >> len.y() >> len.z();
        for (size_t i = 2; i < n + 2; i++)
            s2p(v.at(i), p.at(i - 2));
        return true;
    }
    return false;
}

void FormatGRO::s2p(const std::string &s, Particle &dst) {
    std::stringstream o;
    std::string name;
    o << s.substr(10, 5) << s.substr(20, 8) << s.substr(28, 8) << s.substr(36, 8);
    o >> name;
    auto it = findName(atoms, name);
    if (it != atoms.end()) {
        o >> dst.pos.x() >> dst.pos.y() >> dst.pos.z();
        dst = *it;
    } else
        throw std::runtime_error("gro: unknown atom name");
    dst.pos = 10 * dst.pos; // nm->angstrom
}

/**
 * @param file Filename
 * @param particle_vector Destination particle vector
 */
bool FormatGRO::load(const std::string &file, ParticleVector &particle_vector) {
    particle_vector.clear();
    v.resize(0);
    if (IO::readFile(file, v)) {
        int last = atoi(v[1].c_str()) + 1;
        for (int i = 2; i <= last; i++) {
            Particle a;
            s2p(v[i], a);
            particle_vector.push_back(a);
        }
        return true;
    }
    return false;
}

std::vector<Particle> fastaToParticles(const std::string &fasta_sequence, double bond_length, const Point &origin) {
    ParticleVector particles; // particle vector
    particles.reserve(fasta_sequence.size());
    Particle p;                                      // single particle
    p.pos = origin;                                  // first atom placed at `origin`
    for (auto id : fastaToAtomIds(fasta_sequence)) { // loop over atom ids
        p.id = id;
        p.charge = Faunus::atoms[id].charge;
        if (not particles.empty())
            p.pos = particles.back().pos + ranunit(random) * bond_length;
        particles.push_back(p);
    }
    return particles;
}

/**
 * @param file filename to load (aam, pqr, xyz, ...)
 * @param particle_vector destination particle vector
 * @param append if true, expand dst vector
 * @param keep_charges if true, ignore AtomData charges
 * @return true if successfully loaded
 */
bool loadStructure(const std::string &file, ParticleVector &particle_vector, bool append, bool keep_charges) {
    if (append == false)
        particle_vector.clear();
    std::string suffix = file.substr(file.find_last_of(".") + 1);
    if (suffix == "aam")
        FormatAAM::load(file, particle_vector, keep_charges);
    if (suffix == "pqr")
        FormatPQR::load(file, particle_vector, keep_charges);
    if (suffix == "xyz")
        FormatXYZ::load(file, particle_vector);
    return !particle_vector.empty();
}

/**
 * @param fasta_sequence FASTA sequence, capital letters.
 * @return vector of verified and existing atom id's
 */
std::vector<int> fastaToAtomIds(const std::string &fasta_sequence) {
    std::map<char, std::string> map = {{'A', "ALA"},
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

    for (auto c : fasta_sequence) { // loop over letters
        auto it = map.find(c);      // is it in map?
        if (it == map.end())
            throw std::runtime_error("Invalid FASTA letter '" + std::string(1, c) + "'");
        else
            names.push_back(it->second);
    }
    return Faunus::names2ids(atoms, names);
}

} // namespace Faunus
