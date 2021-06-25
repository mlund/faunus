#include "io.h"
#include "units.h"
#include "random.h"
#include "group.h"
#include "space.h"
#include <spdlog/spdlog.h>
#include <zstr.hpp>
#include <cereal/archives/binary.hpp>

void Faunus::StructureFileReader::handleChargeMismatch(Particle& particle, const int atom_index) {
    if (fabs(particle.traits().charge - particle.charge) > pc::epsilon_dbl) {
        faunus_logger->warn("charge mismatch on atom {0} {1}: {2} atomlist value", atom_index, particle.traits().name,
                            (prefer_charges_from_file) ? "ignoring" : "using");
        if (not prefer_charges_from_file) {
            particle.charge = particle.traits().charge; // let's use atomdata charge
        }
    }
}
void Faunus::StructureFileReader::handleRadiusMismatch(const Particle& particle, const double radius,
                                                       const int atom_index) {
    if (std::fabs(particle.traits().sigma - 2.0 * radius) > pc::epsilon_dbl) {
        faunus_logger->warn("radius mismatch on atom {0} {1}: using atomlist value", atom_index,
                            particle.traits().name);
    }
}
namespace Faunus {

/**
 * Create an output file stream. If the filename extension is `.gz`, gzip
 * compression is enabled, otherwise a standard `std::ostream` is created.
 *
 * @param filename Name of output file
 * @param throw_on_error Throw exception if file cannot be opened (default: false)
 * @return pointer to stream; nullptr if it could not be created
 */
std::unique_ptr<std::ostream> IO::openCompressedOutputStream(const std::string& filename, bool throw_on_error) {
    try { // any neater way to check if path is writable?
        std::ofstream f;
        f.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        f.open(filename);
    } catch (std::exception& e) { // reaching here, file cannot be created
        if (throw_on_error) {
            throw std::runtime_error("could not writeKeyValuePairs to file "s + filename);
        }
        return std::make_unique<std::ofstream>(); // empty object
    }
    if (filename.substr(filename.find_last_of('.') + 1) == "gz") {
        faunus_logger->trace("enabling gzip compression for {}", filename);
        return std::make_unique<zstr::ofstream>(filename); // compressed
    } else {
        return std::make_unique<std::ofstream>(filename); // uncompressed
    }
}

TEST_CASE("[Faunus] openCompressedOutputStream") {
    using doctest::Approx;
    CHECK_THROWS(IO::openCompressedOutputStream("/../file", true));
    CHECK_NOTHROW(IO::openCompressedOutputStream("/../file"));
}

bool PQRTrajectoryReader::readAtomRecord(const std::string& record, Particle& particle, double& radius) {
    std::stringstream o(record);
    std::string key;
    o >> key;
    if (key == "ATOM" or key == "HETATM") {
        int atom_index;
        int res_index;
        std::string atom_name;
        std::string res_name;
        o >> atom_index >> atom_name;
        const auto atom = findAtomByName(atom_name);
        particle = atom;
        o >> res_name >> res_index >> particle.pos.x() >> particle.pos.y() >> particle.pos.z() >> particle.charge >>
            radius;
        return true;
    }
    return false;
}

/**
 * @param filename PQR trajectory filename
 * @param destination Destination vector of particle vector pointers
 *
 * Each trajectory must be separated by an "END" record. Blank lines are
 * ignored.
 */
void PQRTrajectoryReader::loadTrajectory(const std::string& filename, std::vector<ParticleVector>& destination) {
    if (std::ifstream stream(filename); stream) {
        destination.clear();
        destination.resize(1); // prepare first frame
        std::string record;
        while (std::getline(stream, record)) {
            Particle particle;
            double radius = 0.0;
            if (readAtomRecord(record, particle, radius)) {
                destination.back().push_back(particle);
            } else if (record.find("END") == 0) {                       // if END record, advance to next frame
                destination.back().shrink_to_fit();                     // attempt to clean up
                destination.push_back(ParticleVector());                // prepare next frame
                destination.back().reserve(destination.front().size()); // reserve memory
            }
        }
        if (destination.back().empty()) { // delete empty frame after last "END" record
            destination.pop_back();
        }
        if (destination.empty()) {
            faunus_logger->warn("pqr trajectory {} is empty", filename);
        }
    } else {
        throw std::runtime_error("cannot open file");
    }
}

TEST_CASE("[Faunus] StructureFileReader and Writer") {
    using doctest::Approx;
    // Space object with two salt pairs, i.e. four particles
    Space spc;
    SpaceFactory::makeNaCl(spc, 2, R"( {"type": "cuboid", "length": [20,30,40]} )"_json);

    // set positions
    double displacement = 0.0;
    for (int i = 0; i < 4; i++) {
        spc.p[i].pos = {displacement, displacement + 0.1, displacement + 0.2};
        spc.p[i].charge = double(i);
        displacement += 0.5;
    }

    auto io_roundtrip = [&](StructureFileReader& reader, StructureFileWriter& writer) {
        std::stringstream stream;
        writer.save(stream, spc.p.begin(), spc.p.end(), spc.geo.getLength());
        stream.seekg(0); // rewind
        reader.load(stream);
        if (reader.box_dimension_support) {
            CHECK(reader.box_length.x() == Approx(spc.geo.getLength().x()));
            CHECK(reader.box_length.y() == Approx(spc.geo.getLength().y()));
            CHECK(reader.box_length.z() == Approx(spc.geo.getLength().z()));
        }
        CHECK(reader.particles.size() == 4);
        for (int i = 0; i < 4; i++) {
            CHECK(reader.particles[i].pos.x() == Approx(spc.p[i].pos.x()));
            CHECK(reader.particles[i].pos.y() == Approx(spc.p[i].pos.y()));
            CHECK(reader.particles[i].pos.z() == Approx(spc.p[i].pos.z()));
            if (reader.particle_charge_support) {
                CHECK(reader.particles[i].charge == Approx(spc.p[i].charge));
            }
        }
    };

    SUBCASE("PQR roundtrip") {
        PQRReader reader;
        PQRWriter writer;
        io_roundtrip(reader, writer);
    }

    SUBCASE("XYZ roundtrip") {
        XYZReader reader;
        XYZWriter writer;
        io_roundtrip(reader, writer);
    }

    SUBCASE("AAM roundtrip") {
        AminoAcidModelReader reader;
        AminoAcidModelWriter writer;
        io_roundtrip(reader, writer);
    }

    SUBCASE("Gromacs roundtrip") {
        GromacsReader reader;
        GromacsWriter writer;
        io_roundtrip(reader, writer);
        writer.save("/Users/mikael/test.gro", spc.p.begin(), spc.p.end(), spc.geo.getLength());
    }
}

// ========== XTCTrajectoryFrame ==========

XTCTrajectoryFrame::XTCTrajectoryFrame(int number_of_atoms) { initNumberOfAtoms(number_of_atoms); }

XTCTrajectoryFrame::XTCTrajectoryFrame(const TrajectoryFrame& frame) {
    initNumberOfAtoms(frame.coordinates.size());
    importFrame(frame);
}

void XTCTrajectoryFrame::operator=(const TrajectoryFrame& frame) {
    if (frame.coordinates.size() != number_of_atoms) {
        throw std::runtime_error("wrong number of particles to be assign into the XTC frame");
    }
    importFrame(frame);
}

void XTCTrajectoryFrame::importFrame(const TrajectoryFrame& frame) {
    importTimestamp(frame.step, frame.timestamp);
    importBox(frame.box);
    importCoordinates(frame.coordinates, 0.5 * frame.box);
}

void XTCTrajectoryFrame::importTimestamp(const int step, const float time) {
    xtc_step = step;
    xtc_time = time / 1.0_ps;
}

void XTCTrajectoryFrame::importBox(const Point& box) {
    // empty box tensor
    XTCMatrix xtc_box_matrix = XTCMatrix::Zero();
    // only XYZ dimensions in nanometers on diagonal, as floats
    xtc_box_matrix.diagonal() = (box / 1.0_nm).cast<XTCFloat>();
    // copy underlaying eigen structure (1D array, row-major) to the C-style 2D array
    std::copy(xtc_box_matrix.data(), xtc_box_matrix.data() + DIM * DIM, &(xtc_box[0][0]));
}

void XTCTrajectoryFrame::importCoordinates(const PointVector& coordinates, const Point& offset) {
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

void XTCTrajectoryFrame::exportFrame(TrajectoryFrame& frame) const {
    exportTimestamp(frame.step, frame.timestamp);
    exportBox(frame.box);
    exportCoordinates(frame.coordinates, 0.5 * frame.box);
}

void XTCTrajectoryFrame::exportTimestamp(int& step, float& time) const {
    step = xtc_step;
    time = xtc_time * 1.0_ps;
}

void XTCTrajectoryFrame::exportBox(Point& box) const {
    XTCMatrix xtc_box_matrix = Eigen::Map<const XTCTrajectoryFrame::XTCMatrix>(&(xtc_box[0][0]));
    if (xtc_box_matrix.diagonal().asDiagonal().toDenseMatrix() != xtc_box_matrix) {
        throw std::runtime_error("cannot load non-orthogonal box");
    }
    box = Point(xtc_box_matrix.diagonal().cast<double>() * 1.0_nm);
}

void XTCTrajectoryFrame::exportCoordinates(PointVector& coordinates, const Point& offset) const {
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

TrajectoryFrame::TrajectoryFrame(const Point& box, const PointVector& coordinates, int step, float timestamp)
    : box(box), coordinates(coordinates), step(step), timestamp(timestamp) {}

TrajectoryFrame::TrajectoryFrame(const XTCTrajectoryFrame& xtc_frame) {
    coordinates.resize(xtc_frame.number_of_atoms);
    xtc_frame.exportFrame(*this);
}

void TrajectoryFrame::operator=(const XTCTrajectoryFrame& xtc_frame) { xtc_frame.exportFrame(*this); }

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

XTCReader::XTCReader(const std::string& filename) : filename(filename) {
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

bool XTCReader::read(TrajectoryFrame& frame) {
    bool is_ok = readFrame();
    if (is_ok) {
        frame = *xtc_frame;
    }
    return is_ok;
}

// ========== XTCWriter ==========

XTCWriter::XTCWriter(const std::string& filename) : filename(filename) {
    xdrfile = XDRfile::xdrfile_open(filename.c_str(), "w");
    if (xdrfile == nullptr) {
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

void XTCWriter::write(const TrajectoryFrame& frame) {
    if (!xtc_frame) {
        xtc_frame = std::make_shared<XTCTrajectoryFrame>(frame.coordinates.size());
    }
    *xtc_frame = frame;
    writeFrame();
    step_counter = frame.step + 1;
}

void XTCWriter::writeNext(const TrajectoryFrame& frame) {
    if (!xtc_frame) {
        xtc_frame = std::make_shared<XTCTrajectoryFrame>(frame.coordinates.size());
    }
    *xtc_frame = frame;
    writeFrameAt(step_counter, step_counter * time_delta);
    ++step_counter;
}

ParticleVector fastaToParticles(const std::string& fasta_sequence, double bond_length, const Point& origin) {
    ParticleVector particles;                  // particle vector
    auto ids = fastaToAtomIds(fasta_sequence); // convert letters to atom ids
    std::transform(ids.begin(), ids.end(), std::back_inserter(particles), [&](auto& id) {
        Particle particle = Faunus::atoms.at(id);
        particle.pos = particles.empty() ? origin : particles.back().pos + ranunit(random) * bond_length;
        return particle;
    });
    return particles;
}

ParticleVector loadStructure(const std::string& filename, bool prefer_charges_from_file) {
    ParticleVector particles;
    std::unique_ptr<StructureFileReader> reader;
    const auto suffix = filename.substr(filename.find_last_of('.') + 1);
    if (suffix == "aam") {
        reader = std::make_unique<AminoAcidModelReader>();
    } else if (suffix == "pqr") {
        reader = std::make_unique<PQRReader>();
    } else if (suffix == "xyz") {
        reader = std::make_unique<XYZReader>();
    } else if (suffix == "gro") {
        reader = std::make_unique<GromacsReader>();
    }
    try {
        if (reader) {
            reader->prefer_charges_from_file = prefer_charges_from_file;
            particles = reader->load(filename);
        } else {
            throw std::runtime_error("unknown format");
        }
        if (particles.empty()) {
            throw std::runtime_error("empty structure");
        }
        return particles;
    } catch (std::exception& e) { throw std::runtime_error(fmt::format("{} load error -> {}", filename, e.what())); }
}

/**
 * @param fasta_sequence FASTA sequence, capital letters.
 * @return vector of verified and existing atom id's
 */
std::vector<int> fastaToAtomIds(const std::string& fasta_sequence) {
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

    for (const auto letter : fasta_sequence) {             // loop over letters
        if (auto it = map.find(letter); it != map.end()) { // is it in map?
            names.push_back(it->second);
        } else {
            throw std::runtime_error("invalid FASTA letter '" + std::string(1, letter) + "'");
        }
    }
    return Faunus::names2ids(atoms, names);
}

FormatSpaceTrajectory::FormatSpaceTrajectory(std::ostream& ostream) {
    if (ostream) {
        output_archive = std::make_unique<cereal::BinaryOutputArchive>(ostream);
    }
}
FormatSpaceTrajectory::FormatSpaceTrajectory(std::istream& istream) {
    if (istream) {
        input_archive = std::make_unique<cereal::BinaryInputArchive>(istream);
    }
}
void FormatSpaceTrajectory::load([[maybe_unused]] Space& spc) { assert(input_archive != nullptr); }
void FormatSpaceTrajectory::save([[maybe_unused]] const Space& spc) { assert(output_archive != nullptr); }

// ------------------------

size_t StructureFileReader::getNumberOfAtoms(const std::string& line) {
    try {
        return std::stoul(line);
    } catch (std::exception& e) { throw std::invalid_argument("invalid number of particles"); }
}

void StructureFileReader::getNextLine(std::istream& stream, std::string& line) {
    std::getline(stream, line);
    if (!line.empty()) {
        if (line.substr(0, 1) == "#") {
            comments.push_back(line);
            getNextLine(stream, line);
        }
    }
}

ParticleVector& StructureFileReader::load(std::istream& stream) {
    stream.exceptions(std::istream::failbit | std::istream::badbit);
    loadHeader(stream);
    particles.clear();
    if (expected_number_of_particles > 0) {
        particles.reserve(expected_number_of_particles);
    }
    while (!stream.eof()) {
        try {
            particles.emplace_back(loadParticle(stream));
        } catch (std::istream::failure& e) {
            break; // end of stream reached
        } catch (std::exception& e) { throw std::runtime_error(e.what()); }
    }
    particles.shrink_to_fit();
    loadFooter(stream);
    if (expected_number_of_particles > 0 && expected_number_of_particles != particles.size()) {
        throw std::out_of_range(
            fmt::format("expected {} particle records but found {}", expected_number_of_particles, particles.size()));
    }
    if (particles.empty()) {
        std::cout << "no particles loaded" << std::endl;
        faunus_logger->warn("no particles loaded");
    }
    return particles;
}

ParticleVector& StructureFileReader::load(const std::string& filename) {
    try {
        if (std::ifstream stream(filename); stream) {
            return load(stream);
        }
        throw std::runtime_error("cannot open file");
    } catch (std::exception& e) { throw std::runtime_error(e.what()); }
}

void StructureFileReader::loadFooter([[maybe_unused]] std::istream& stream) {}

// -----------------------------

void AminoAcidModelWriter::saveHeader(std::ostream& stream, int number_of_particles) const {
    stream << number_of_particles << "\n";
}
void AminoAcidModelWriter::saveParticle(std::ostream& stream, const Particle& particle) {
    const auto& traits = particle.traits();
    auto scale = static_cast<double>(particle_is_active) / 1.0_angstrom;
    stream << fmt::format("{:7} {:7d} {:>13.6E} {:>13.6E} {:>13.6E} {:>13.6E} {:9.3f} {:9.3f}\n", traits.name,
                          particle_index + 1, scale * particle.pos[0], scale * particle.pos[1], scale * particle.pos[2],
                          scale * particle.charge, scale * traits.mw, 0.5 * traits.sigma);
}

// -----------------------------

void AminoAcidModelReader::loadHeader(std::istream& stream) {
    std::string line;
    getNextLine(stream, line); // all lines starting w. "#" are skipped
    expected_number_of_particles = getNumberOfAtoms(line);
}

Particle AminoAcidModelReader::loadParticle(std::istream& stream) {
    std::string line;
    getNextLine(stream, line);
    std::stringstream record(line);
    record.exceptions(std::ios::failbit);
    std::string particle_name;
    record >> particle_name;
    int particle_index;
    double radius;
    double molecular_weight;
    const auto& atomdata = Faunus::findAtomByName(particle_name);
    auto particle = Particle(atomdata);
    try {
        record >> particle_index >> particle.pos.x() >> particle.pos.y() >> particle.pos.z() >> particle.charge >>
            molecular_weight >> radius;
    } catch (std::exception& e) { throw std::runtime_error("syntax error in: " + line); }
    particle.pos *= 1.0_angstrom;
    handleChargeMismatch(particle, particle_index);
    handleRadiusMismatch(particle, radius, particle_index);
    return particle;
}

AminoAcidModelReader::AminoAcidModelReader() {
    particle_charge_support = true;
    particle_radius_support = true;
}

// -----------------------------

void XYZWriter::saveHeader(std::ostream& stream, int number_of_particles) const {
    stream << number_of_particles << "\n"
           << "Generated by Faunus - https://github.com/mlund/faunus\n";
}
void XYZWriter::saveParticle(std::ostream& stream, const Particle& particle) {
    stream << particle.traits().name << " " << particle.pos.transpose() << "\n";
}

// -----------------------------

void XYZReader::loadHeader(std::istream& stream) {
    std::string line;
    try {
        std::getline(stream, line);
        expected_number_of_particles = std::stoi(line);
    } catch (std::exception& e) { throw std::invalid_argument("atom count expected on first line"); }
    try {
        auto flags = stream.exceptions();     // backup state flags
        stream.exceptions(std::ios::goodbit); // do not throw if empty comment line
        std::getline(stream, line);
        stream.exceptions(flags); // restore original state flags
        if (!line.empty()) {
            comments.push_back(line);
            faunus_logger->debug("xyz-file comment: {}", line);
        }
    } catch (std::exception& e) { throw std::invalid_argument("cannot load comment"); }
}

Particle XYZReader::loadParticle(std::istream& stream) {
    std::string line;
    getNextLine(stream, line);
    std::stringstream record(line);
    record.exceptions(std::ios::failbit);
    std::string atom_name;
    record >> atom_name;
    const auto& atom = Faunus::findAtomByName(atom_name);
    auto particle = Particle(atom);
    record >> particle.pos.x() >> particle.pos.y() >> particle.pos.z();
    particle.pos *= 1.0_angstrom; // xyz files are commonly in Ångstroms
    return particle;
}
// -----------------------------

/**
 * Read line by line until ATOM or HETATM entry is found or EOF is reached
 * @throw if eof is reached or if syntax error in ATOM or HETATM entry
 */
Particle PQRReader::loadParticle(std::istream& stream) {
    std::string line;
    std::string type;
    std::getline(stream, line); // throws if eof or empty
    std::stringstream record(line);
    record.exceptions(std::ios::failbit); // throws if syntax error
    record >> type;
    if (type == "ATOM" || type == "HETATM") {
        double radius;
        std::size_t atom_index;
        std::size_t residue_index;
        std::string atom_name;
        std::string residue_name;
        record >> atom_index >> atom_name;
        Particle particle = Faunus::findAtomByName(atom_name);
        record >> residue_name >> residue_index >> particle.pos.x() >> particle.pos.y() >> particle.pos.z() >>
            particle.charge >> radius;
        particle.pos -= 0.5 * box_length; //
        handleChargeMismatch(particle, atom_index);
        handleRadiusMismatch(particle, radius, atom_index);
        return particle;
    }
    return loadParticle(stream); // keep trying until eof...
}

/**
 * Load line by line until eof or first ATOM/HETATM is reached
 * @throw if eof
 */
void PQRReader::loadHeader(std::istream& stream) {
    std::string line;
    std::string type;
    while (!stream.eof()) {
        auto original_position = stream.tellg();
        try {
            std::getline(stream, line);
        } catch (std::istream::failure& e) { return; }
        if (!line.empty()) {
            std::stringstream record(line);
            record.exceptions(std::ios::failbit); // throws if error
            try {
                record >> type;
                if (type == "REMARK") {
                    comments.push_back(line);
                } else if (type == "CRYST1") {
                    record >> box_length.x() >> box_length.y() >> box_length.z();
                    box_length *= 1.0_angstrom;
                    faunus_logger->debug("read box dimensions: {} {} {}", box_length.x(), box_length.y(),
                                         box_length.z());
                } else if (type == "ATOM" || type == "HETATM") {
                    stream.seekg(original_position); // rewind and stop
                    break;
                }
            } catch (std::exception& e) { throw std::runtime_error("header -> " + line + " ("s + e.what() + ")"s); }
        }
    }
}
PQRReader::PQRReader() {
    particle_charge_support = true;
    particle_radius_support = true;
    box_dimension_support = true;
}

// -----------------------------

void PQRWriter::saveParticle(std::ostream& stream, const Particle& particle) {
    const auto scale = static_cast<double>(particle_is_active);
    const auto position = scale * (particle.pos + 0.5 * box_dimensions); // origin to corner of box
    auto atom_name = particle.traits().name;
    if (atom_name.size() > 4) {
        atom_name.erase(4);
    }
    if (group_name.size() > 3) {
        group_name.erase(3);
    } else if (group_name.empty()) {
        group_name = "n/a";
    }
    const std::string chain = "A";
    if (style == PQR) {
        stream << fmt::format("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   "
                              "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n",
                              "ATOM", particle_index + 1, atom_name, "A", group_name, chain, group_index + 1, "0",
                              position.x(), position.y(), position.z(), scale * particle.charge,
                              scale * particle.traits().sigma * 0.5);

    } else if (style == PDB) { // see https://cupnet.net/pdb-format
        const double occupancy = 0.0;
        const double temperature_factor = 1.0;
        const std::string element_symbol = atom_name;
        const std::string charge = "0";
        stream << fmt::format("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   "
                              "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n",
                              "ATOM", particle_index + 1, atom_name, "A", group_name, chain, group_index + 1, "0",
                              position.x(), position.y(), position.z(), occupancy, temperature_factor, element_symbol,
                              charge);
    } else { // legacy PQR
        stream << fmt::format("ATOM  {:5d} {:4} {:4}{:5d}    {:8.3f} {:8.3f} {:8.3f} {:.3f} {:.3f}\n",
                              particle_index + 1, atom_name, group_name, group_index + 1, position.x(), position.y(),
                              position.z(), scale * particle.charge, scale * particle.traits().sigma * 0.5);
    }
}

void PQRWriter::saveHeader(std::ostream& stream, [[maybe_unused]] int number_of_particles) const {
    if (box_dimensions.squaredNorm() > pc::epsilon_dbl) {
        const Point angle = {90.0, 90.0, 90.0};
        const Point box = box_dimensions / 1.0_angstrom;
        stream << fmt::format("CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n", box.x(), box.y(),
                              box.z(), angle.x(), angle.y(), angle.z());
    }
}

void PQRWriter::saveFooter(std::ostream& stream) const { stream << "END\n"; }
PQRWriter::PQRWriter(PQRWriter::Style style) : style(style) {}

// -----------------------------

void StructureFileWriter::saveFooter([[maybe_unused]] std::ostream& stream) const {}

void StructureFileWriter::saveGroup(std::ostream& stream, const Group<Particle>& group) {
    group_name = group.traits().name;
    for (auto particle = group.begin(); particle != group.trueend(); particle++) { // loop over particles
        particle_is_active = particle < group.end();
        saveParticle(stream, *particle);
        particle_index++;
    }
    group_index++;
}

std::shared_ptr<StructureFileWriter> createStructureFileWriter(const std::string& suffix) {
    std::shared_ptr<StructureFileWriter> writer;
    if (suffix == "aam") {
        writer = std::make_shared<AminoAcidModelWriter>();
    } else if (suffix == "pqr") {
        writer = std::make_shared<PQRWriter>();
    } else if (suffix == "xyz") {
        writer = std::make_shared<XYZWriter>();
    } else if (suffix == "gro") {
        writer = std::make_shared<GromacsWriter>();
    }
    return writer;
}

// -----------------------

void GromacsWriter::saveHeader(std::ostream& stream, int number_of_particles) const {
    stream << "# Generated by Faunus - https://github.com/mlund/faunus\n" << number_of_particles << "\n";
}

void GromacsWriter::saveFooter(std::ostream& stream) const {
    if (box_dimensions.squaredNorm() > pc::epsilon_dbl) {
        stream << box_dimensions.transpose() / 1.0_nm << "\n";
    }
}

void GromacsWriter::saveParticle(std::ostream& stream, const Particle& particle) {
    const auto& atom_name = particle.traits().name;
    auto scale = static_cast<double>(particle_is_active);
    Point position = scale * (particle.pos + 0.5 * box_dimensions) / 1.0_nm; // shift origin and convert to nm
    stream << fmt::format("{:5d}{:5}{:5}{:5d}{:8.3f}{:8.3f}{:8.3f}\n", group_index, atom_name, atom_name,
                          particle_index, position.x(), position.y(), position.z());
}

// -----------------------

GromacsReader::GromacsReader() { box_dimension_support = true; }

/**
 * This also loads the *footer* as box information is placed after the coordinates,
 * but we need it to offset while loading particles.
 */
void GromacsReader::loadHeader(std::istream& stream) {
    std::string line;
    std::getline(stream, line);
    comments.push_back(line);
    std::getline(stream, line);
    try {
        expected_number_of_particles = std::stoi(line);
    } catch (std::exception& e) { throw std::invalid_argument("atom count expected on second line"); }
    loadBoxInformation(stream);
}

void GromacsReader::loadBoxInformation(std::istream& stream) {
    auto initial_position = stream.tellg();
    stream.seekg(-2, std::istream::end); // jump to end and skip possible newline
    for (int i = stream.tellg(); i > 0; i--) {
        if (stream.peek() != '\n') {
            stream.seekg(-1, std::istream::cur); // one step backwards
        } else {
            stream.get(); // gobble newline
            std::string line;
            std::getline(stream, line);
            std::stringstream o(line);
            o.exceptions(std::stringstream::failbit);
            try {
                o >> box_length.x() >> box_length.y() >> box_length.z();
                box_length *= 1.0_nm;
            } catch (...) { faunus_logger->debug("no box information found in gro file"); }
            break;
        }
    }
    stream.seekg(initial_position);
}

Particle GromacsReader::loadParticle(std::istream& stream) {
    if (particles.size() >= expected_number_of_particles) {
        throw std::istream::failure("cannot read beyond expected number of particles");
    }
    std::string record;
    std::getline(stream, record);
    std::stringstream o;
    o << record.substr(10, 5) << record.substr(20, 8) << record.substr(28, 8) << record.substr(36, 8);
    std::string atom_name;
    o >> atom_name;
    Particle particle = Particle(findAtomByName(atom_name));
    o >> particle.pos.x() >> particle.pos.y() >> particle.pos.z();
    particle.pos = particle.pos * 1.0_nm - 0.5 * box_length;
    return particle;
}
} // namespace Faunus
