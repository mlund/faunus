#pragma once

#include "particle.h"
#include "spdlog/spdlog.h"
#include "units.h"
#include <fstream>

namespace cereal {
class BinaryOutputArchive;
class BinaryInputArchive;
} // namespace cereal

namespace Faunus {

class Space;

template <typename T = Particle> class Group;

#ifndef __cplusplus
#define __cplusplus
#endif
#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif

namespace XDRfile {
#include <xdrfile_trr.h>
#include <xdrfile_xtc.h>
} // namespace XDRfile

namespace IO {

/**
 * @brief Open (gzip compressed) output stream
 */
std::unique_ptr<std::ostream> openCompressedOutputStream(const std::string&, bool throw_on_error = false);

/**
 * Write a map to an output stream as key-value pairs
 * @tparam TKey
 * @tparam TValue
 * @param stream Output stream
 * @param data
 */
template <typename TKey, typename TValue>
void writeKeyValuePairs(std::ostream& stream, const std::map<TKey, TValue>& data, const std::string& sep = " ",
                        const std::string& end = "\n") {
    if (stream) {
        for (const auto& [key, value] : data) {
            stream << key << sep << value << end;
        }
    }
}

/**
 * Write a map to a file as key-value pairs
 * @tparam TKey
 * @tparam TValue
 * @param filename
 * @param data
 */
template <typename TKey, typename TValue>
void writeKeyValuePairs(const std::string& filename, const std::map<TKey, TValue>& data) {
    if (!data.empty()) {
        std::ofstream file(filename);
        writeKeyValuePairs(file, data);
    }
}

} // namespace IO

/**
 * Base class to load simple structure files such as XYZ, AAM, PQR etc.
 */
class StructureFileReader {
  private:
    virtual void loadHeader(std::istream& stream) = 0;       //!< Gobble entire header
    virtual void loadFooter(std::istream& stream);           //!< Gobble entire header
    virtual Particle loadParticle(std::istream& stream) = 0; //!< Load single particle
    void checkLoadedParticles() const;                       //!< Checks if expected number of particles were loaded

  protected:
    static size_t getNumberOfAtoms(const std::string& line); //!< Helper function to extract N
    void getNextLine(std::istream& stream, std::string& destination,
                     const std::string& comment_identifiers); //!< Helper function to forward stream
    size_t expected_number_of_particles = 0;

    void handleChargeMismatch(Particle& particle, int atom_index) const; //!< Policy if charge mismatch
    static void handleRadiusMismatch(const Particle& particle, double radius,
                                     int atom_index); //!< Policy if radius mismatch

  public:
    ParticleVector particles;
    Point box_length = Point::Zero();
    std::vector<std::string> comments;
    bool prefer_charges_from_file = true; //!< If applicable, prefer charges from AAM file over `AtomData`

    ParticleVector& load(std::istream& stream); //!< Load entire stream and populate data
    ParticleVector& load(const std::string& filename);
    virtual ~StructureFileReader() = default;

    bool box_dimension_support = false;
    bool particle_charge_support = false;
    bool particle_radius_support = false;
};

/**
 * Reads FASTA file and generate one particle per amino acid
 *
 * Positions are generated as a random walk beginning from the
 * given `initial_particle_position`, advancing in steps of `bond_length`
 */
class CoarseGrainedFastaFileReader : public StructureFileReader {
  private:
    void loadHeader(std::istream& stream) override;
    Particle loadParticle(std::istream& stream) override;
    double bond_length = 0.0;
    Point new_particle_position = Point::Zero();

  public:
    explicit CoarseGrainedFastaFileReader(double bond_length, const Point& initial_particle_position = Point(0, 0, 0));
    void setBondLength(double bond_length);
    std::string loadSequence(std::istream& stream);
    static char getFastaLetter(std::istream& stream);
};

class AminoAcidModelReader : public StructureFileReader {
  private:
    void loadHeader(std::istream& stream) override;
    Particle loadParticle(std::istream& stream) override;

  public:
    AminoAcidModelReader();
};

class PQRReader : public StructureFileReader {
  private:
    void loadHeader(std::istream& stream) override;
    Particle loadParticle(std::istream& stream) override;

  public:
    PQRReader();
};

/**
 * @brief XYZ file loader
 *
 * This format has number of particles at the first line;
 * comment on the second line; and positions of named particles on the following lines.
 *
 * Example:
 *
 *     3
 *     Three point water model (comment)
 *     OW  2.30  6.28  1.13
 *     HW  1.37  6.26  1.50
 *     HW  2.31  5.89  0.21
 */
class XYZReader : public StructureFileReader {
  private:
    void loadHeader(std::istream& stream) override;
    Particle loadParticle(std::istream& stream) override;
};

class GromacsReader : public StructureFileReader {
  private:
    void loadBoxInformation(std::istream& stream); //!< Load box dimensions (stream position is preserved)
    void loadHeader(std::istream& stream) override;
    Particle loadParticle(std::istream& stream) override;

  public:
    GromacsReader();
};

/**
 * Base class to writeKeyValuePairs simple structure files such as XYZ, AAM, PQR etc.
 */
class StructureFileWriter {
  private:
    virtual void saveHeader(std::ostream& stream, int number_of_particles) const = 0; //!< Write header
    virtual void saveFooter(std::ostream& stream) const; //!< Called when all particles have been written
    virtual void saveParticle(std::ostream& stream, const Particle& particle) = 0; //!< Write single particle
    void saveGroup(std::ostream& stream, const Group<Particle>& group);            //!< Write entire group

    template <class ParticleIter> void saveParticles(std::ostream& stream, ParticleIter begin, ParticleIter end) {
        group_index = 0;
        particle_index = 0;
        std::for_each(begin, end, [&](const auto& particle) {
            saveParticle(stream, particle);
            particle_index++;
        });
    }

  protected:
    static const std::string generated_by_faunus_comment;
    bool particle_is_active = true;
    std::string group_name;
    std::size_t particle_index = 0;
    std::size_t group_index = 0;
    Point box_dimensions = Point::Zero();

  public:
    template <class ParticleIter>
    void save(std::ostream& stream, ParticleIter begin, ParticleIter end, const Point& box_length) {
        if (auto number_of_particles = std::distance(begin, end); number_of_particles > 0) {
            box_dimensions = box_length;
            saveHeader(stream, number_of_particles);
            saveParticles(stream, begin, end);
            saveFooter(stream);
        }
    }

    template <typename Range> void save(std::ostream& stream, const Range& groups, const Point& box_length) {
        group_index = 0;
        particle_index = 0;
        box_dimensions = box_length;

        int number_of_particles = 0;
        std::for_each(groups.begin(), groups.end(),
                      [&](const auto& group) { number_of_particles += group.capacity(); });
        saveHeader(stream, number_of_particles);

        std::for_each(groups.begin(), groups.end(), [&](const auto& group) { saveGroup(stream, group); });
        saveFooter(stream);
    }

    template <class... Args> void save(const std::string& filename, const Args&... args) {
        if (std::ofstream stream(filename); stream) {
            faunus_logger->debug("writing to {}", filename);
            save(stream, args...);
        } else {
            throw std::runtime_error("writeKeyValuePairs error: "s + filename);
        }
    }

    virtual ~StructureFileWriter() = default;
};

class AminoAcidModelWriter : public StructureFileWriter {
  private:
    void saveHeader(std::ostream& stream, int number_of_particles) const override;
    void saveParticle(std::ostream& stream, const Particle& particle) override;
};

class XYZWriter : public StructureFileWriter {
  private:
    void saveHeader(std::ostream& stream, int number_of_particles) const override;
    void saveParticle(std::ostream& stream, const Particle& particle) override;
};

class PQRWriter : public StructureFileWriter {
  private:
    void saveHeader(std::ostream& stream, int number_of_particles) const override;
    void saveFooter(std::ostream& stream) const override;
    void saveParticle(std::ostream& stream, const Particle& particle) override;

  public:
    enum Style { PQR_LEGACY, PDB, PQR }; //!< PQR style (for ATOM records)
    Style style = PQR_LEGACY;
    explicit PQRWriter(Style style = PQR_LEGACY);
};

class GromacsWriter : public StructureFileWriter {
  private:
    void saveHeader(std::ostream& stream, int number_of_particles) const override;
    void saveFooter(std::ostream& stream) const override;
    void saveParticle(std::ostream& stream, const Particle& particle) override;
};

namespace PQRTrajectoryReader {
bool readAtomRecord(const std::string& record, Particle& particle, double& radius); //!< Read ATOM or HETATOM record
void loadTrajectory(const std::string& filename, std::vector<ParticleVector>& destination); //!< Load trajectory
} // namespace PQRTrajectoryReader

struct TrajectoryFrame;

/**
 * @brief Base data structure for native XTC format as used in Gromacs and supported by the C library. Import methods
 * do the data conversion from the Faunus native format to the XTC format, and export methods do the oposite.
 *
 * By convention, XTC has coordinates' origin in a corner (main box's coordinates are always positive), while Faunus
 * has coordinates' origin in the center of the simulation box. During conversion the corresponding offset is
 * subtracted, or added, respectively.
 *
 * XTC format uses floats (Faunus doubles) and dimensions are in nanometers (Faunus ångströms). XTC library requires
 * raw 2D C-style arrays for coordinates in row-major format. XTC tensor of the simulation box is converted
 * to an XYZ point pressuming orthogonal geometry bacause of current limitation of Faunus. XTC format does not support
 * variable number of coordinates (atoms) between frames.
 */
struct XTCTrajectoryFrame {
    XDRfile::matrix xtc_box;                          //!< box tensor; only diagonal elements are used
    std::unique_ptr<XDRfile::rvec[]> xtc_coordinates; //!< C-style array of particle coordinates
    int xtc_step = 0;                                 //!< current frame number
    float xtc_time = 0.0;                             //!< current time (unit?)
    int number_of_atoms = 0;                          //!< number of coordinates (atoms) in each frame
    float precision = 1000.0;                         //!< output precision
    /**
     * @brief Creates an empty XTC trajectory frame for given number of coordinates (atoms).
     * @param number_of_atoms  number of coordinates (atoms)
     */
    XTCTrajectoryFrame(int number_of_atoms);
    /**
     * @brief Creates an XTC trajectory frame from TrajectoryFrame and converts data accordingly.
     * @param frame  source trajectory frame
     */
    XTCTrajectoryFrame(const TrajectoryFrame& frame);
    /**
     * @brief Creates an XTC trajectory frame from scalar parameters and input iterator and converts data accordingly.
     * @tparam begin_iterator
     * @tparam end_iterator
     * @param[in] step  frame step
     * @param[in] time  timestamp in picoseconds
     * @param[in] box  box dimensions (xyz) in nanometers
     * @param[in] coordinates_begin  input iterator with coordinates in nanometers
     * @param[in] coordinates_end  input iterator's end
     */
    template <class begin_iterator, class end_iterator>
    XTCTrajectoryFrame(const int step, const float time, const Point box, begin_iterator coordinates_begin,
                       end_iterator coordinates_end) {
        initNumberOfAtoms(std::distance(coordinates_begin, coordinates_end));
        importFrame(step, time, box, coordinates_begin, coordinates_end);
    }
    /**
     * @brief Copies TrajectoryFrame and converts data. Calls importFrame.
     *
     * The number of coordinates is immutable to prevent mistakes.
     * @param frame  source frame
     * @throw std::runtime_error  when the number of coordinates does not match
     */
    void operator=(const TrajectoryFrame& frame);
    /**
     * @brief Imports data from a TrajectoryFrame.
     * @param frame  source frame
     * @throw std::runtime_error  when the number of coordinates does not match
     */
    void importFrame(const TrajectoryFrame& frame);
    /**
     * @brief Imports data from scalar parameters and an input iterator over coordinates.
     * @tparam begin_iterator
     * @tparam end_iterator
     * @param[in] step  frame step
     * @param[in] time  frame timestamp
     * @param[in] box  box dimensions (xyz) in nanometers
     * @param[in] coordinates_begin  input iterator with coordinates in nanometers
     * @param[in] coordinates_end  input iterator's end
     * @throw std::runtime_error  when the number of coordinates does not match
     */
    template <class begin_iterator, class end_iterator>
    void importFrame(const int step, const float time, const Point box, begin_iterator coordinates_begin,
                     end_iterator coordinates_end) {
        importTimestamp(step, time);
        importBox(box);
        importCoordinates(coordinates_begin, coordinates_end, 0.5 * box);
    }
    /**
     * @brief Exports data from a TrajectoryFrame.
     * @param frame  target frame
     * @throw std::runtime_error  when the number of coordinates does not match
     */
    void exportFrame(TrajectoryFrame& frame) const;
    /**
     * @brief Exports data to scalar paramers and output iterator over atomic coordinates.
     * @tparam begin_iterator
     * @tparam end_iterator
     * @param[out] step  frame step
     * @param[out] time  frame timestamp
     * @param[out] box  box dimensions (xyz) in nanometers
     * @param[out] coordinates_begin  output iterator with coordinates in nanometers
     * @param[out] coordinates_end  output iterator's end
     * @throw std::runtime_error  when the number of coordinates does not match
     */
    template <class begin_iterator, class end_iterator>
    void exportFrame(int& step, float& time, Point& box, begin_iterator coordinates_begin,
                     end_iterator coordinates_end) const {
        exportTimestamp(step, time);
        exportBox(box);
        exportCoordinates(coordinates_begin, coordinates_end, 0.5 * box);
    }

  protected:
    typedef float XTCFloat;
    typedef Eigen::Matrix<XTCFloat, DIM, DIM, Eigen::RowMajor> XTCMatrix; //<! eigen equivalent of the box tensor
    typedef Eigen::Matrix<XTCFloat, DIM, 1> XTCVector; //<! eigen equivalent of the single set of coordinates
    /**
     * @brief Imports and converts step and timestamp.
     * @param[in] step  frame step
     * @param[in] time  frame timestamp
     */
    void importTimestamp(const int step, const float time);
    /**
     * @brief Imports and converts simulation box dimensions.
     * @param[in] box  simulation box dimensions in nanometers (xyz)
     */
    void importBox(const Point& box);
    /**
     * @brief Imports and converts atomic coordinates. Offset is added to all coordinates to account different
     * coordinates' origin.
     * @param[in] coordinates  atomic coordinates in nanometers
     * @param[in] offset  offset in nanometers to add to all coordinates upon conversion
     */
    void importCoordinates(const PointVector& coordinates, const Point& offset);
    /**
     * @brief Imports and converts atomic coordinates. Offset is added to all coordinates to account different
     * coordinates' origin.
     * @tparam begin_iterator
     * @tparam end_iterator
     * @param[in] begin  input iterator with coordinates in nanometers
     * @param[in] end  input iterator's end
     * @param[in] offset  offset in nanometers to add to all coordinates upon conversion
     */
    template <class begin_iterator, class end_iterator>
    void importCoordinates(begin_iterator begin, end_iterator end, const Point& offset) const {
        // comparison of i is probably faster than prior call of std::distance
        int i = 0;
        for (auto coordinates_it = begin; coordinates_it != end; ++coordinates_it) {
            if (i >= number_of_atoms) {
                throw std::runtime_error("number of particles too high to be saved in the XTC frame");
            }
            const Point pos = (*coordinates_it + offset) / 1.0_nm;
            const XTCVector xtc_pos = pos.cast<XTCFloat>();
            std::copy(xtc_pos.data(), xtc_pos.data() + DIM, xtc_coordinates.get()[i++]);
        }
        if (i < number_of_atoms) {
            throw std::runtime_error("number of particles too low to be saved in the XTC frame");
        }
    }
    /**
     * @brief Exports and converts step and timestamp.
     * @param[out] step  frame step
     * @param[out] time  frame timestamp
     */
    void exportTimestamp(int& step, float& time) const;
    /**
     * @brief Exports and converts simulation box dimensions.
     * @param[out] box  simulation box dimensions in nanometers (xyz)
     */
    void exportBox(Point& box) const;
    /**
     * @brief Exports and converts atomic coordinates. Offset is subtracted from all coordinates to account different
     * coordinates' origin.
     * @param[out] coordinates  atomic coordinates in nanometers
     * @param[in] offset  offset in nanometers to subtract from all coordinates upon conversion
     * @throw std::runtime_error  when the source box is not orthogonal
     */
    void exportCoordinates(PointVector& coordinates, const Point& offset) const;
    /**
     * @brief Exports and converts atomic coordinates. Offset is subtracted to all coordinates to account different
     * coordinates' origin.
     * @tparam begin_iterator
     * @tparam end_iterator
     * @param[out] begin  output iterator with coordinates in nanometers
     * @param[out] end  output iterator's end
     * @param[in] offset  offset in nanometers to subtract from all coordinates upon conversion
     */
    template <class begin_iterator, class end_iterator>
    void exportCoordinates(begin_iterator begin, end_iterator end, const Point& offset) const {
        // comparison of i is probably faster than prior call of std::distance
        int i = 0;
        for (auto coordinates_it = begin; coordinates_it != end; ++coordinates_it) {
            if (i >= number_of_atoms) {
                throw std::runtime_error("number of particles in the loaded XTC frame too low");
            }
            XTCVector xtc_atom_coordinates(xtc_coordinates.get()[i++]);
            *coordinates_it = Point(xtc_atom_coordinates.cast<double>() * 1.0_nm) - offset;
        }
        if (i < number_of_atoms) {
            throw std::runtime_error("number of particles in the loaded XTC frame too high");
        }
    }

  private:
    /**
     * @brief Set number of coordinates (atoms) and allocate memory for them.
     */
    void initNumberOfAtoms(int);
};

/**
 * @brief Simple data structure to store a trajectory frame in native Faunus format.
 *
 * All coordinates are in ånström, origin is placed into the geometric center of the simulation box and timestamps are
 * in picoseconds.
 */
struct TrajectoryFrame {
    Point box;               //!< simulation box
    PointVector coordinates; //!< coordinates of particles
    int step = 0;            //!< frame number
    float timestamp = 0.0;   //!< frame timestamp

    TrajectoryFrame() = default;
    TrajectoryFrame(const Point& box, const PointVector& coordinates, int step, float timestamp);
    /**
     * @brief Creates a new trajectory frame based on an XTC frame. Data are converted as necessary. Convenient wrapper
     * around XTCTrajectoryFrame::exportFrame.
     * @param xtc_frame  source XTC trajectory frame
     */
    TrajectoryFrame(const XTCTrajectoryFrame& xtc_frame);
    /**
     * @brief Assignes an XTC frame. Data are converted as necessary. However, the number of coordinates has to be
     * the same in both (source and target) frames. Convinient wrapper around XTCTrajectoryFrame::exportFrame.
     * @param xtc_frame  source XTC trajectory frame
     */
    void operator=(const XTCTrajectoryFrame& xtc_frame);
};

/**
 * @brief Reads frames from an XTC file (GROMACS compressed trajectory file format). It is a wrapper around
 * C function calls.
 *
 * Frames are stored into a TrajectoryFrame structure or as a list of positions in an output iterator. The class
 * is responsible for I/O operations, not data conversion. For details about data conversion XTCTrajectoryFrame.
 */
class XTCReader {
    int return_code = XDRfile::exdrOK;   //!< last return code of a C function
    XDRfile::XDRFILE* xdrfile = nullptr; //!< file handle
    //! data structure for C functions; the number of coordinates is immutable
    std::unique_ptr<XTCTrajectoryFrame> xtc_frame;

  public:
    std::string filename; //!< name of the trajectory file, mainly for error reporting
    /**
     * @param filename  a name of the XTC file to open
     */
    XTCReader(const std::string& filename);
    ~XTCReader();
    /**
     * @brief Returns number of coordinates (atoms) in each frame. Immutable during object lifetime.
     * @return number of coordinates (atoms)
     */
    int getNumberOfCoordinates();
    /**
     * @brief Reads the next frame in the trajectory from a file
     * @param[out] frame   target frame
     * @return true on success, false at the end of file
     * @throw std::runtime_error  when other I/O error occures
     */
    bool read(TrajectoryFrame& frame);
    /**
     * @brief Reads the next frame in the trajectory.
     * @tparam begin_iterator
     * @tparam end_iterator
     * @param[out] step  frame step
     * @param[out] time  frame timestamp in picoseconds
     * @param[out] box  box dimensions (xyz) in nanometers
     * @param[out] coordinates_begin  output iterator to store coordinates
     * @param[out] coordinates_end  output iterator's end
     * @return true on success, false at the end of file
     * @throw std::runtime_error  when other I/O error occures
     */
    template <class begin_iterator, class end_iterator>
    bool read(int& step, float& time, Point& box, begin_iterator coordinates_begin, end_iterator coordinates_end) {
        bool is_ok = readFrame();
        if (is_ok) {
            xtc_frame->exportFrame(step, time, box, coordinates_begin, coordinates_end);
        }
        return is_ok;
    }

  protected:
    /**
     * @brief Actual wrapper around C function that reads the next frame into xtc_frame.
     * @return true on success, false at the end of file
     * @throw std::runtime_error  when other I/O error occures
     */
    bool readFrame();
};

/**
 * @brief Writes frames into an XTC file (GROMACS compressed trajectory file format). It is a wrapper around
 * C function calls.
 *
 * The frames can be provided as a TrajectoryFrame structure or as a list of positions in an input iterator. The class
 * is responsible for I/O operations, not data conversion.
 */
class XTCWriter {
    int return_code = XDRfile::exdrOK;             //!< last return code of a C function
    XDRfile::XDRFILE* xdrfile = nullptr;           //!< file handle
    std::unique_ptr<XTCTrajectoryFrame> xtc_frame; //!< data structure for C functions;
                                                   //!< the number of coordinates is immutable
    int step_counter = 0;                          //!< frame counter for automatic increments
    float time_delta = 1.0_ps;                     //!< timestamp of a frame is computed as step * time_delta

  public:
    std::string filename; //!< name of the trajectory file, mainly for error reporting
    /**
     * @param filename  a name of the XTC file to open
     */
    XTCWriter(const std::string& filename);
    ~XTCWriter();
    /**
     * @brief Writes a frame into the file.
     * @param[in] frame  frame to be written
     * @throw std::runtime_error  when other I/O error occures
     */
    void write(const TrajectoryFrame& frame);
    /**
     * @brief Writes a next frame into the file using own automatic counter for step and timestamp.
     * The corresponding values in the frame are ignored.
     * @param[in] frame  frame to be written
     * @throw std::runtime_error  when other I/O error occures
     */
    void writeNext(const TrajectoryFrame& frame);
    /**
     * @brief Writes a next frame into the file using own automatic counter for step and timestamp.
     * @tparam begin_iterator
     * @tparam end_iterator
     * @param[in] box  dimensions of the cubic box (xyz)
     * @param[in] coordinates_begin  input iterator with coordinates (not particles)
     * @param[in] coordinates_end  input iterator's end
     */
    template <class begin_iterator, class end_iterator>
    void writeNext(const Point& box, begin_iterator coordinates_begin, end_iterator coordinates_end) {
        if (!xtc_frame) {
            auto number_of_atoms = std::distance(coordinates_begin, coordinates_end);
            xtc_frame = std::make_unique<XTCTrajectoryFrame>(number_of_atoms);
        }
        xtc_frame->importFrame(step_counter, step_counter * time_delta, box, coordinates_begin, coordinates_end);
        writeFrame();
        ++step_counter;
    }

  protected:
    /**
     * @brief Actual wrapper around C function that writes the current frame into xtc_frame.
     * @throw std::runtime_error  when other I/O error occures
     */
    void writeFrame();
    /**
     * @brief Actual wrapper around C function that writes the current frame into xtc_frame overriding step
     * and timestamp.
     * @throw std::runtime_error  when other I/O error occures
     */
    void writeFrameAt(int step, float time);
};

std::vector<int> fastaToAtomIds(const std::string&); //!< Convert FASTA sequence to atom id sequence

/**
 * @brief Create particle vector from FASTA sequence with equally spaced atoms
 *
 * Particle positions are generated as a random walk, starting at `origin`,
 * propagating in `bond_length` steps.
 */
ParticleVector fastaToParticles(const std::string& fasta_sequence, double bond_length = 7.0,
                                const Point& origin = {0, 0, 0});

/**
 * @brief Load structure file into particle vector
 * @param filename filename to load (aam, pqr, xyz, gro)
 * @param prefer_charges_from_file if true, ignore AtomData charges
 * @throws Throws exception if nothing was loaded or if unknown suffix
 * @returns particles destination particle vector (will be overwritten)
 */
ParticleVector loadStructure(const std::string& filename, bool prefer_charges_from_file = true);

/**
 * @brief Create structure writer
 * @param suffix Filename suffix (pqr, pdb, aam, xyz, gro)
 * @return Shared pointer to write instance; empty if unknown suffix
 */
std::shared_ptr<StructureFileWriter> createStructureFileWriter(const std::string& suffix);

/**
 * @brief Placeholder for Space Trajectory
 *
 * The idea is that the format handles both input and
 * output streams that may of may not be compressed.
 */
class FormatSpaceTrajectory {
  private:
    std::unique_ptr<cereal::BinaryOutputArchive> output_archive;
    std::unique_ptr<cereal::BinaryInputArchive> input_archive;

  public:
    FormatSpaceTrajectory(std::ostream& ostream);
    FormatSpaceTrajectory(std::istream& istream);
    void load(Space&);       //!< Load single frame from stream
    void save(const Space&); //!< Save single frame from stream
};

} // namespace Faunus
