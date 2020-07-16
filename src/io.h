#pragma once

#include "units.h"
#include "core.h"
#include "particle.h"
#include "spdlog/spdlog.h"
#include <cereal/archives/binary.hpp>
#include <fstream>
#include <range/v3/distance.hpp>

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
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"
} // namespace XDRfile

namespace IO {

/**
 * @brief Read lines from file into vector
 * @param filename Filename to read
 * @param destination Reference to vector to load into
 * @return True if file was opened
 */
bool readFile(const std::string &filename, std::vector<std::string> &destination); //!< Read lines from file into vector

bool writeFile(const std::string &filename, const std::string &s,
               std::ios_base::openmode mode = std::ios_base::out); //!< Write string to file

void strip(std::vector<std::string> &string_vector, const std::string &pattern); //!< Strip lines matching a pattern

/**
 * @brief Open (gzip compressed) output stream
 */
std::unique_ptr<std::ostream> openCompressedOutputStream(const std::string &);

/**
 * Write a map to an output stream as key-value pairs
 * @tparam TKey
 * @tparam TValue
 * @param stream Output stream
 * @param data
 */
template <typename TKey, typename TValue>
void write(std::ostream &stream, const std::map<TKey, TValue> &data, const std::string &sep = " ",
           const std::string &end = "\n") {
    if (stream) {
        for (auto [key, value] : data) {
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
template <typename TKey, typename TValue> void write(const std::string &filename, const std::map<TKey, TValue> &data) {
    if (!data.empty()) {
        std::ofstream file(filename);
        write(file, data);
    }
}

} // namespace IO

/**
 * @brief Read/write AAM file format
 *
 * The AAM format is a simple format for loading particle positions
 * charges, radii and molecular weights. The structure is as follows:
 *
 * - Lines beginning with # are ignored and can be placed anywhere
 * - The first non-# line gives the number of particles
 * - Every subsequent line gives atom information in the format:
 *
 *   `name number x y z charge  weight radius`
 *
 * - Positions and radii should be in angstroms.
 * - Currently, data in the number field is ignored.
 * - No particular spacing is required.
 *
 * Example:
 *
 *     2
 *     Na    1     10.234 5.4454 -2.345  +1    22.0   1.7
 *     Cl    2    5.011     1.054  20.02   -1   35.0   2.0
 */
class FormatAAM {
  private:
    static bool prefer_charges_from_file; // true of we prefer charges from AAM file over AtomData
    static std::string p2s(const Particle &, int);
    static Particle &s2p(const std::string &, Particle &); // convert string line to particle

  public:
    static bool load(const std::string &, ParticleVector &, bool = true);
    static bool save(const std::string &, const ParticleVector &);
};

/**
 * @brief Create and read PQR files
 * @date December 2007
 */
class FormatPQR {
  private:
    static std::string writeCryst1(const Point &, const Point & = {90, 90, 90}); //!< Write CRYST1 record
    static bool readCrystalRecord(const std::string &, Point &);                 //!< Read CRYST1 record
    static bool readAtomRecord(const std::string &, Particle &, double &);       //!< Read ATOM or HETATOM record

  public:
    typedef std::vector<Group<Particle>> Tgroup_vector;
    static Point load(std::istream &, ParticleVector &, bool);                      //!< Load PQR from stream
    static Point load(const std::string &, ParticleVector &, bool);                 //!< Load PQR from file
    static void loadTrajectory(const std::string &, std::vector<ParticleVector> &); //!< Load trajectory
    static bool save(std::ostream &, const ParticleVector &, const Point & = Point(0, 0, 0),
                     int = 1e9); //!< Save PQR file
    static bool save(const std::string &, const ParticleVector &, const Point & = Point(0, 0, 0),
                     int = 1e9); //!< Save PQR file

    static bool save(std::ostream &, const Tgroup_vector &, const Point & = Point(0, 0, 0));      //!< Save PQR file
    static bool save(const std::string &, const Tgroup_vector &, const Point & = Point(0, 0, 0)); //!< Save PQR file
};

/**
 * @brief XYZ format
 * @date June 2013
 *
 * Saves particles as a XYZ file. This format has number of particles at
 * the first line; comment on second line (box dimensions) followed by positions of all
 * particles xyz position on each line
 */
struct FormatXYZ {
    static bool save(const std::string &, const ParticleVector &, const Point & = {0, 0, 0});           //!< Save XYZ
    static bool load(const std::string &filename, ParticleVector &particle_vector, bool append = true); //!< Load XYZ
};

/**
 * @brief Gromacs GRO format
 * @date December 2007
 */
class FormatGRO {
  private:
    static Particle recordToParticle(const std::string &); //!< Parse record to particle

  public:
    /**
     * @brief Load GRO file into particle vector
     * @param file Filename
     * @returns Destination particle vector
     */
    static ParticleVector load(const std::string &);
    static bool save(const std::string &filename, const Space &spc);
};

/**
 * @brief GROMACS xtc compressed trajectory file format
 *
 * Saves simulation frames to a Gromacs xtc trajectory file including
 * box information if applicable.
 *
 * @date June 2007-2011, Prague / Malmo
 * @note Alternative pure C++ version(?):
 * http://loos.sourceforge.net/xtc_8hpp_source.html
 */
class FormatXTC {
  private:
    XDRfile::matrix box;                        //!< box dimensions
    XDRfile::XDRFILE *xdrfile = nullptr;        //!< file handle
    std::shared_ptr<XDRfile::rvec> coordinates; //!< vector of particle coordinates
    int step_counter = 0;                       //!< Current number of time steps
    float timestamp = 0.0;                      //!< Current time (unit?)
    int number_of_atoms = 0;                    //!< Atoms in trajectory
    float precision = 1000.0;                   //!< Output precision

  public:
    int getNumAtoms();

    /**
     * @brief Load a single frame into cuboid
     *
     * This will read a single frame from the xtc file (must be open) into
     * a Cuboid container. The box dimensions are retrieved for the frame and transfered
     * to the container. Coordinates are copied into both the particle vector "p" and the
     * "trial" vector. In doing so, positions are converted from nm to angstroms and the
     * coordinate system is shifted so that origin is on the middle of the box. As a safefy
     * measure we do a container collision check to see if all particles are within
     * the Cuboid boundaries.
     *
     * @note The container particle vector *must* match the number of particles
     *       in the xtc file. If not
     *       an error message will be issued and the function will abort.
     */
    bool loadNextFrame(Space &spc, bool setbox = true, bool apply_periodic_boundaries = false);

    /**
     * This will take an arbitrary particle vector and add it
     * to an xtc file. If the file is already open, coordinates will
     * be added, while a new file is created if not.
     * Coordinates are shifted and converted to nanometers.
     * Box dimensions for the frame must be manually
     * set by the `setLength()` function before calling this.
     */
    template <class begin_iterator, class end_iterator /** particle vector iterator */>
    bool save(const std::string &file, begin_iterator begin, end_iterator end) {
        if (begin != end) {
            if (!xdrfile) {
                xdrfile = XDRfile::xdrfile_open(&file[0], "w");
            }
            if (xdrfile) {
                std::shared_ptr<XDRfile::rvec> coords(new XDRfile::rvec[ranges::distance(begin, end)],
                                                      std::default_delete<XDRfile::rvec[]>());
                size_t N = 0;
                for (auto particle = begin; particle != end; ++particle) {
                    coords.get()[N][0] = particle->pos.x() / 1.0_nm + box[0][0] * 0.5;
                    coords.get()[N][1] = particle->pos.y() / 1.0_nm + box[1][1] * 0.5;
                    coords.get()[N][2] = particle->pos.z() / 1.0_nm + box[2][2] * 0.5;
                    N++;
                }
                XDRfile::write_xtc(xdrfile, N, step_counter++, timestamp++, box, coords.get(), precision);
                return true;
            }
        }
        return false;
    }

    /**
     * This will open an xtc file for reading. The number of atoms in each frame
     * is saved and memory for the coordinate array is allocated.
     */

    FormatXTC(double);
    ~FormatXTC();
    bool open(std::string);
    void close();
    void setLength(const Point &);
};

std::vector<int> fastaToAtomIds(const std::string &); //!< Convert FASTA sequence to atom id sequence

/**
 * @brief Create particle vector from FASTA sequence with equally spaced atoms
 *
 * Particle positions are generated as a random walk, starting at `origin`,
 * propagating in `bond_length` steps.
 */
ParticleVector fastaToParticles(const std::string &fasta_sequence, double bond_length = 7.0,
                                const Point &origin = {0, 0, 0});

/**
 * @brief Load structure file into particle vector
 */
bool loadStructure(const std::string &, std::vector<Particle> &, bool, bool = true);

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
    FormatSpaceTrajectory(std::ostream &ostream);
    FormatSpaceTrajectory(std::istream &istream);
    void load(Space &);       //!< Load single frame from stream
    void save(const Space &); //!< Save single frame from stream
};

} // namespace Faunus
