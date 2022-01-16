#include "io.h"
#include "units.h"
#include "random.h"
#include "group.h"
#include "space.h"
#include "multipole.h"
#include <spdlog/spdlog.h>
#include <zstr.hpp>
#include <cereal/archives/binary.hpp>

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
            throw std::runtime_error("could not open file "s + filename);
        }
        return nullptr;
    }
    if (filename.substr(filename.find_last_of('.') + 1) == "gz") {
        faunus_logger->trace("enabling gzip compression for {}", filename);
        return std::make_unique<zstr::ofstream>(filename); // compressed
    } else {
        return std::make_unique<std::ofstream>(filename); // uncompressed
    }
}

TEST_CASE("[Faunus] openCompressedOutputStream") {
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
                destination.emplace_back();                             // prepare next frame
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

// ========== XTCTrajectoryFrame ==========

XTCTrajectoryFrame::XTCTrajectoryFrame(int number_of_atoms) { initNumberOfAtoms(number_of_atoms); }

XTCTrajectoryFrame::XTCTrajectoryFrame(const TrajectoryFrame& frame) {
    initNumberOfAtoms(static_cast<int>(frame.coordinates.size()));
    importFrame(frame);
}

XTCTrajectoryFrame& XTCTrajectoryFrame::operator=(const TrajectoryFrame& frame) {
    if (frame.coordinates.size() != number_of_atoms) {
        throw std::runtime_error("wrong number of particles to be assign into the XTC frame");
    }
    importFrame(frame);
    return *this;
}

void XTCTrajectoryFrame::importFrame(const TrajectoryFrame& frame) {
    importTimestamp(frame.step, frame.timestamp);
    importBox(frame.box);
    importCoordinates(frame.coordinates, 0.5 * frame.box);
}

void XTCTrajectoryFrame::importTimestamp(const int step, const float time) {
    xtc_step = step;
    xtc_time = time / static_cast<float>(1.0_ps);
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
    time = xtc_time * static_cast<float>(1.0_ps);
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
    for (int i = 0; i < number_of_atoms; ++i) {
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
        xtc_frame = std::make_unique<XTCTrajectoryFrame>(number_of_atoms);
        xdrfile = XDRfile::xdrfile_open(filename.c_str(), "r");
    }
    if (!xtc_frame || (xdrfile == nullptr)) {
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
    if (readFrame()) {
        frame = *xtc_frame;
        return true;
    }
    return false;
}

// ========== XTCWriter ==========

XTCWriter::XTCWriter(const std::string& filename) : xdrfile(XDRfile::xdrfile_open(filename.c_str(), "w")), filename(filename) {
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

void XTCWriter::write(const TrajectoryFrame& frame) {
    if (!xtc_frame) {
        xtc_frame = std::make_unique<XTCTrajectoryFrame>(frame.coordinates.size());
    }
    *xtc_frame = frame;
    writeFrame();
    step_counter = frame.step + 1;
}

void XTCWriter::writeNext(const TrajectoryFrame& frame) {
    if (!xtc_frame) {
        xtc_frame = std::make_unique<XTCTrajectoryFrame>(frame.coordinates.size());
    }
    *xtc_frame = frame;
    writeFrameAt(step_counter, step_counter * time_delta);
    ++step_counter;
}

ParticleVector fastaToParticles(std::string_view fasta_sequence, double bond_length, const Point& origin) {
    ParticleVector particles;                  // particle vector
    auto ids = fastaToAtomIds(fasta_sequence); // convert letters to atom ids
    std::transform(ids.begin(), ids.end(), std::back_inserter(particles), [&](auto& id) {
        Particle particle = Faunus::atoms.at(id);
        particle.pos = particles.empty() ? origin : particles.back().pos + randomUnitVector(random) * bond_length;
        return particle;
    });
    return particles;
}

std::shared_ptr<StructureFileWriter> createStructureFileWriter(const std::string& suffix) {
    std::shared_ptr<StructureFileWriter> writer;
    if (suffix == "pqr") {
        writer = std::make_shared<PQRWriter>();
    } else if (suffix == "aam") {
        writer = std::make_shared<AminoAcidModelWriter>();
    } else if (suffix == "xyz") {
        writer = std::make_shared<XYZWriter>();
    } else if (suffix == "gro") {
        writer = std::make_shared<GromacsWriter>();
    } else if (suffix == "pdb") {
        writer = std::make_shared<PQRWriter>(PQRWriter::Style::PDB);
    }
    return writer;
}

// -----------------------------

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
std::vector<AtomData::index_type> fastaToAtomIds(std::string_view fasta_sequence) {
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

    for (const auto letter : fasta_sequence) { // loop over letters
        if (letter == '*') {
            break;
        }
        if (letter == ' ' || letter == '\n') {
            continue;
        }
        if (auto it = map.find(letter); it != map.end()) { // is it in map?
            names.push_back(it->second);
        } else {
            throw std::runtime_error("invalid FASTA letter '" + std::string(1, letter) + "'");
        }
    }
    return Faunus::names2ids(Faunus::atoms, names);
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

void StructureFileReader::handleChargeMismatch(Particle& particle, const int atom_index) const {
    if (std::fabs(particle.traits().charge - particle.charge) > pc::epsilon_dbl) {
        faunus_logger->warn("charge mismatch on atom {0} {1}: {2} atomlist value", atom_index, particle.traits().name,
                            (prefer_charges_from_file) ? "ignoring" : "using");
        if (not prefer_charges_from_file) {
            particle.charge = particle.traits().charge; // let's use atomdata charge
        }
    }
}
void StructureFileReader::handleRadiusMismatch(const Particle& particle, const double radius, const int atom_index) {
    if (std::fabs(particle.traits().sigma - 2.0 * radius) > pc::epsilon_dbl) {
        faunus_logger->warn("radius mismatch on atom {0} {1}: using atomlist value", atom_index,
                            particle.traits().name);
    }
}

size_t StructureFileReader::getNumberOfAtoms(const std::string& line) {
    try {
        return std::stoul(line);
    } catch (std::exception& e) { throw std::invalid_argument("invalid number of particles"); }
}

/**
 * @param stream Stream to operate on
 * @param destination Destination string for loaded line
 * @param comment_identifiers Array of characters used to identify comment lines, e.g. "#", ";" etc.
 * @throw If line cannot be loaded
 *
 * If the *first* character in the loaded line matches any character in `comment_identifiers`, the line
 * will be ignored and the next will be loaded.
 */
void StructureFileReader::getNextLine(std::istream& stream, std::string& destination,
                                      const std::string& comment_identifiers) {
    assert(stream.exceptions() & std::ios::failbit); // check that we throw upon failure
    while (true) {
        std::getline(stream, destination);
        if (!destination.empty()) {
            const bool is_comment = comment_identifiers.find_first_of(destination.front()) != std::string::npos;
            if (is_comment) {
                const auto first_non_whitespace = destination.substr(1).find_first_not_of(' ') + 1;
                comments.emplace_back(destination.begin() + first_non_whitespace, destination.end());
                continue;
            }
            break;
        }
    }
}

ParticleVector& StructureFileReader::load(std::istream& stream) {
    comments.clear();
    particles.clear();
    stream.exceptions(std::istream::failbit | std::istream::badbit);
    loadHeader(stream);
    std::for_each(comments.begin(), comments.end(),
                  [](auto& comment) { faunus_logger->debug("comment line: {}", comment); });
    if (expected_number_of_particles > 0) {
        particles.reserve(expected_number_of_particles);
    }
    while (true) {
        try {
            particles.emplace_back(loadParticle(stream));
        } catch (std::istream::failure& e) {
            break; // end of stream reached
        } catch (std::exception& e) { throw std::runtime_error(e.what()); }
    }
    particles.shrink_to_fit();
    loadFooter(stream);
    checkLoadedParticles();
    if (particle_charge_support) {
        faunus_logger->debug("net charge of loaded structure: {:.3}e",
                             Faunus::monopoleMoment(particles.begin(), particles.end()));
    }
    return particles;
}

void StructureFileReader::checkLoadedParticles() const {
    if (expected_number_of_particles > 0 && expected_number_of_particles != particles.size()) {
        throw std::out_of_range(
            fmt::format("expected {} particle records but found {}", expected_number_of_particles, particles.size()));
    }
    if (particles.empty()) {
        std::cout << "no particles loaded" << std::endl;
        faunus_logger->warn("no particles loaded");
    }
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

const std::string StructureFileWriter::generated_by_faunus_comment =
    "Generated by Faunus - https://github.com/mlund/faunus";

void StructureFileWriter::saveFooter([[maybe_unused]] std::ostream& stream) const {}

void StructureFileWriter::saveGroup(std::ostream& stream, const Group& group) {
    group_name = group.traits().name;
    for (auto particle = group.begin(); particle != group.trueend(); particle++) {
        particle_is_active = particle < group.end();
        saveParticle(stream, *particle);
        particle_index++;
    }
    group_index++;
}

TEST_CASE("[Faunus] StructureFileReader and StructureFileWriter") {
    using doctest::Approx;
    // Space object with two salt pairs, i.e. four particles
    Space spc;
    SpaceFactory::makeNaCl(spc, 2, R"( {"type": "cuboid", "length": [20,30,40]} )"_json);

    CHECK(spc.particles.size() == 4);

    // set positions
    double displacement = 0.0;
    for (int i = 0; i < 4; i++) {
        spc.particles[i].pos = {displacement, displacement + 0.1, displacement + 0.2};
        spc.particles[i].charge = double(i);
        displacement += 0.5;
    }

    auto io_roundtrip = [&](StructureFileReader& reader, StructureFileWriter& writer) {
        std::stringstream stream;
        writer.save(stream, spc.particles.begin(), spc.particles.end(), spc.geometry.getLength());
        stream.seekg(0); // rewind
        reader.load(stream);
        if (reader.box_dimension_support) {
            CHECK(reader.box_length.x() == Approx(spc.geometry.getLength().x()));
            CHECK(reader.box_length.y() == Approx(spc.geometry.getLength().y()));
            CHECK(reader.box_length.z() == Approx(spc.geometry.getLength().z()));
        }
        CHECK(reader.particles.size() == 4);
        for (int i = 0; i < 4; i++) {
            CHECK(reader.particles[i].pos.x() == Approx(spc.particles[i].pos.x()));
            CHECK(reader.particles[i].pos.y() == Approx(spc.particles[i].pos.y()));
            CHECK(reader.particles[i].pos.z() == Approx(spc.particles[i].pos.z()));
            if (reader.particle_charge_support) {
                CHECK(reader.particles[i].charge == Approx(spc.particles[i].charge));
            }
        }
        CHECK(!reader.comments.empty());
        CHECK(reader.comments.front() == "Generated by Faunus - https://github.com/mlund/faunus");
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
    }
}

// -----------------------------

void AminoAcidModelWriter::saveHeader(std::ostream& stream, int number_of_particles) const {
    stream << fmt::format("# {}\n{}\n", generated_by_faunus_comment, number_of_particles);
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
    getNextLine(stream, line, "#"); // all lines starting w. "#" are skipped
    expected_number_of_particles = getNumberOfAtoms(line);
}

Particle AminoAcidModelReader::loadParticle(std::istream& stream) {
    std::string line;
    getNextLine(stream, line, "#");
    std::stringstream record(line);
    record.exceptions(std::ios::failbit);
    std::string particle_name;
    record >> particle_name;
    int particle_index = 0;
    double radius = 0.0;
    double molecular_weight = 0.0;
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
    stream << fmt::format("{}\n{}\n", number_of_particles, generated_by_faunus_comment);
}
void XYZWriter::saveParticle(std::ostream& stream, const Particle& particle) {
    const auto scale = static_cast<double>(particle_is_active);
    stream << particle.traits().name << " " << scale * particle.pos.transpose() << "\n";
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
        }
    } catch (std::exception& e) { throw std::invalid_argument("cannot load comment"); }
}

Particle XYZReader::loadParticle(std::istream& stream) {
    std::string line;
    getNextLine(stream, line, ";#");
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
    while (true) {
        std::getline(stream, line); // throws if eof or empty
        std::stringstream record(line);
        record.exceptions(std::ios::failbit); // throws if syntax error
        record >> type;
        if (type == "ATOM" || type == "HETATM") {
            double radius;     // NOLINT(cppcoreguidelines-init-variables)
            int atom_index;    // NOLINT(cppcoreguidelines-init-variables)
            int residue_index; // NOLINT(cppcoreguidelines-init-variables)
            std::string atom_name;
            std::string residue_name;
            record >> atom_index >> atom_name;
            auto particle = Particle(Faunus::findAtomByName(atom_name));
            record >> residue_name >> residue_index >> particle.pos.x() >> particle.pos.y() >> particle.pos.z() >>
                particle.charge >> radius;
            particle.pos = particle.pos * 1.0_angstrom - 0.5 * box_length;
            handleChargeMismatch(particle, atom_index);
            handleRadiusMismatch(particle, radius, atom_index);
            return particle;
        }
    }
}

/**
 * Load line by line until eof or first ATOM/HETATM is reached
 * @throw if eof
 */
void PQRReader::loadHeader(std::istream& stream) {
    std::string line;
    std::string type;
    while (true) {
        const auto original_position = stream.tellg();
        try {
            std::getline(stream, line);
        } catch (std::istream::failure& e) { return; }
        if (!line.empty()) {
            std::stringstream record(line);
            record.exceptions(std::ios::failbit); // throws if error
            try {
                record >> type;
                if (type == "REMARK") {
                    const auto first_non_white_space = line.substr(6).find_first_not_of(' ') + 6;
                    comments.emplace_back(line.begin() + first_non_white_space, line.end());
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
    const auto& atom_name = particle.traits().name;
    const auto scale = static_cast<double>(particle_is_active) / 1.0_angstrom;
    const auto position = scale * (particle.pos + 0.5 * box_dimensions); // origin to corner of box
    const auto chain = "A"s;
    const auto charge = "0"s;
    const auto occupancy = 0.0;
    const auto temperature_factor = 1.0;
    const std::string element_symbol = atom_name;

    if (group_name.empty()) {
        group_name = "n/a";
    }

    switch (style) {
    case Style::PQR:
        stream << fmt::format("{:6s}{:5d} {:^4.4s}{:1s}{:3.3s} {:1s}{:4d}{:1s}   "
                              "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n",
                              "ATOM", particle_index + 1, atom_name, "A", group_name, chain, group_index + 1, "0",
                              position.x(), position.y(), position.z(), scale * particle.charge,
                              scale * particle.traits().sigma * 0.5);
        break;
    case Style::PDB: // see https://cupnet.net/pdb-format
        stream << fmt::format("{:6s}{:5d} {:^4.4s}{:1s}{:3.3s} {:1s}{:4d}{:1s}   "
                              "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2.2s}{:2s}\n",
                              "ATOM", particle_index + 1, atom_name, "A", group_name, chain, group_index + 1, "0",
                              position.x(), position.y(), position.z(), occupancy, temperature_factor, element_symbol,
                              charge);
        break;
    case Style::PQR_LEGACY:
        stream << fmt::format("ATOM  {:5d} {:4.4} {:4.3}{:5d}    {:8.3f} {:8.3f} {:8.3f} {:.3f} {:.3f}\n",
                              particle_index + 1, atom_name, group_name, group_index + 1, position.x(), position.y(),
                              position.z(), scale * particle.charge, scale * particle.traits().sigma * 0.5);
        break;
    default:
        throw std::runtime_error("unimplemented style");
    }
}

void PQRWriter::saveHeader(std::ostream& stream, [[maybe_unused]] int number_of_particles) const {
    stream << "REMARK " << generated_by_faunus_comment << "\n";
    if (box_dimensions.squaredNorm() > pc::epsilon_dbl) {
        const Point angle = {90.0, 90.0, 90.0};
        const Point box = box_dimensions / 1.0_angstrom;
        stream << fmt::format("CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n", box.x(), box.y(),
                              box.z(), angle.x(), angle.y(), angle.z());
    }
}

void PQRWriter::saveFooter(std::ostream& stream) const { stream << "END\n"; }

PQRWriter::PQRWriter(PQRWriter::Style style) : style(style) {}

// -----------------------

void GromacsWriter::saveHeader(std::ostream& stream, const int number_of_particles) const {
    stream << fmt::format("{}\n{:5d}\n", generated_by_faunus_comment, number_of_particles);
}

void GromacsWriter::saveFooter(std::ostream& stream) const {
    if (box_dimensions.squaredNorm() > pc::epsilon_dbl) {
        stream << fmt::format("{:10.5f}{:10.5f}{:10.5f}\n", box_dimensions.x() / 1.0_nm, box_dimensions.y() / 1.0_nm,
                              box_dimensions.z() / 1.0_nm);
    }
}

void GromacsWriter::saveParticle(std::ostream& stream, const Particle& particle) {
    const auto& atom_name = particle.traits().name;
    const auto scale = static_cast<double>(particle_is_active) / 1.0_nm; // zero if inactive
    Point position = scale * (particle.pos + 0.5 * box_dimensions);      // shift origin
    stream << fmt::format("{:5d}{:5.5}{:>5.5}{:5d}{:8.3f}{:8.3f}{:8.3f}\n", group_index + 1, group_name, atom_name,
                          particle_index + 1, position.x(), position.y(), position.z());
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
    expected_number_of_particles = getNumberOfAtoms(line);
    loadBoxInformation(stream);
}

void GromacsReader::loadBoxInformation(std::istream& stream) {
    auto initial_position = stream.tellg();
    stream.seekg(-2, std::istream::end); // jump to end and skip possible newline and the very end
    while (stream.tellg() > 0) {         // read file backwards, letter by letter
        if (stream.peek() != '\n') {
            stream.seekg(-1, std::ios::cur); // one step backwards
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
    std::string atom_name;
    std::string record;
    std::getline(stream, record);
    std::stringstream o;
    o.exceptions(std::stringstream::failbit);
    o << record.substr(10, 5) << " " << record.substr(20, 8) << " " << record.substr(28, 8) << " "
      << record.substr(36, 8);
    o >> atom_name;
    auto particle = Particle(findAtomByName(atom_name));
    o >> particle.pos.x() >> particle.pos.y() >> particle.pos.z();
    particle.pos = particle.pos * 1.0_nm - 0.5 * box_length; // origin --> middle of the box
    return particle;
}

// -----------------------

/**
 * @brief Advance stream position past all initial lines starting with ";" or ">"
 * @todo See future literal `z`: https://en.cppreference.com/w/cpp/language/integer_literal
 */
void CoarseGrainedFastaFileReader::loadHeader(std::istream& stream) {
    std::string non_comment_line;
    getNextLine(stream, non_comment_line, ">;");
    const auto position_of_non_comment_line = -static_cast<std::ios::pos_type>(non_comment_line.size() + 1U);
    stream.seekg(position_of_non_comment_line, std::ios::cur); // rewind
}

char CoarseGrainedFastaFileReader::getFastaLetter(std::istream& stream) {
    while (true) {
        const auto letter = static_cast<char>(stream.get());
        switch (letter) {
        case ' ': // ignore white space
            continue;
        case '\n': // ignore newlines
            continue;
        case '*': // signals end of sequence -> stop
            throw std::istream::failure("end of sequence (*) reached");
        case '>': // stop and warn if multiple sequences
            faunus_logger->warn("multiple FASTA sequences detected; using only the first");
            throw std::istream::failure("stopping after first sequence");
        default:
            return letter;
        }
    }
}

Particle CoarseGrainedFastaFileReader::loadParticle(std::istream& stream) {
    assert(stream.exceptions() & std::ios::failbit); // check that we throw upon failure
    const auto letter = getFastaLetter(stream);
    const auto atomid = fastaToAtomIds(std::string(1, letter)).at(0);
    auto particle = Particle(Faunus::atoms.at(atomid));
    particle.pos = new_particle_position;
    new_particle_position += Faunus::randomUnitVector(random) * bond_length;
    return particle;
}

CoarseGrainedFastaFileReader::CoarseGrainedFastaFileReader(const double bond_length,
                                                           const Point& initial_particle_position)
    : bond_length(bond_length), new_particle_position(initial_particle_position) {}

void CoarseGrainedFastaFileReader::setBondLength(const double bond_length) {
    CoarseGrainedFastaFileReader::bond_length = bond_length;
}

std::string CoarseGrainedFastaFileReader::loadSequence(std::istream& stream) {
    stream.exceptions(std::istream::failbit | std::istream::badbit);
    loadHeader(stream);
    std::string sequence;
    while (true) {
        try {
            sequence.push_back(getFastaLetter(stream));
        } catch (std::istream::failure& e) {
            break; // end of stream reached
        } catch (std::exception& e) { throw std::runtime_error(e.what()); }
    }
    return sequence;
}

TEST_CASE("[Faunus] CoarseGrainedFastaFileReader") {
    Faunus::atoms = R"([{ "ARG": { } }, { "LYS": { } } ])"_json.get<decltype(atoms)>();
    const auto bond_length = 7.0;

    SUBCASE("single sequence") {
        CoarseGrainedFastaFileReader reader(bond_length, {10, 20, 30});
        std::stringstream stream(">MCHU - Calmodulin\n; another comment\nKR\n");
        reader.load(stream);
        CHECK(reader.particles.size() == 2);
        CHECK(reader.particles[0].id == 1);
        CHECK(reader.particles[1].id == 0);
        CHECK(reader.particles[0].pos.x() == doctest::Approx(10));
        CHECK(reader.particles[0].pos.y() == doctest::Approx(20));
        CHECK(reader.particles[0].pos.z() == doctest::Approx(30));
        Point distance = reader.particles[0].pos - reader.particles[1].pos;
        CHECK(distance.norm() == doctest::Approx(bond_length));
        CHECK(reader.comments.size() == 2);
        CHECK(reader.comments.at(0) == "MCHU - Calmodulin");
        CHECK(reader.comments.at(1) == "another comment");
    }
    SUBCASE("single sequence - end by star") {
        CoarseGrainedFastaFileReader reader(bond_length);
        std::stringstream stream(">MCHU - Calmodulin\nR*K\n");
        reader.load(stream);
        CHECK(reader.particles.size() == 1);
        CHECK(reader.comments.size() == 1);
        CHECK(reader.comments.at(0) == "MCHU - Calmodulin");
    }
    SUBCASE("multisequence") {
        CoarseGrainedFastaFileReader reader(bond_length);
        std::stringstream stream(">MCHU - Calmodulin\nKKK\n>RK*\n");
        reader.load(stream);
        CHECK(reader.particles.size() == 3);
        CHECK(reader.comments.size() == 1);
        CHECK(reader.comments.at(0) == "MCHU - Calmodulin");
    }
    SUBCASE("unknown particle") {
        CoarseGrainedFastaFileReader reader(bond_length);
        std::stringstream stream(">MCHU - Calmodulin\nRKE*\n");
        CHECK_THROWS(reader.load(stream));
    }
    SUBCASE("unknown fasta letter") {
        CoarseGrainedFastaFileReader reader(bond_length);
        std::stringstream stream(">MCHU - Calmodulin\nRKx*\n");
        CHECK_THROWS(reader.load(stream));
    }
    SUBCASE("no comment") {
        CoarseGrainedFastaFileReader reader(bond_length);
        std::stringstream stream("RKK*\n");
        reader.load(stream);
        CHECK(reader.particles.size() == 3);
        CHECK(reader.comments.empty());
    }
}

// -----------------------

} // namespace Faunus
