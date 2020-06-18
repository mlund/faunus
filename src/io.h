#pragma once

#include "core.h"
#include "particle.h"
#include "spdlog/spdlog.h"
#include <cereal/archives/binary.hpp>
#include <fstream>
#include <range/v3/distance.hpp>

namespace Faunus {

template <typename T = Particle> class Group;

#ifndef __cplusplus
#define __cplusplus
#endif
#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"

namespace IO {

bool readFile(const std::string &, std::vector<std::string> &); //!< Read lines from file into vector

bool writeFile(const std::string &file, const std::string &s,
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
void write(std::ostream &stream, const std::map<TKey, TValue> &data, const std::string &sep = " ", const std::string &end = "\n") {
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
    static bool keepcharges; // true of we prefer charges from AAM file over AtomData
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
    static bool save(std::ostream &, const ParticleVector &, Point = Point(0, 0, 0), int = 1e9);      //!< Save PQR file
    static bool save(const std::string &, const ParticleVector &, Point = Point(0, 0, 0), int = 1e9); //!< Save PQR file

    static bool save(std::ostream &, const Tgroup_vector &, Point = Point(0, 0, 0));      //!< Save PQR file
    static bool save(const std::string &, const Tgroup_vector &, Point = Point(0, 0, 0)); //!< Save PQR file
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
 * @todo Non cubic dimensions
 */
class FormatGRO {
  private:
    std::vector<std::string> v;
    void s2p(const std::string &, Particle &);

  public:
    double boxlength; //!< Box side length (cubic so far)

    /**
     * @brief Load GRO file into particle vector
     * @param file Filename
     * @param p Destination particle vector
     */
    bool load(const std::string &, ParticleVector &);

    template <class Tspace> bool static save(const std::string &filename, const Tspace &spc) {
        if (std::ofstream file(filename); bool(file)) {
            int nres = 1, natom = 1;
            Point boxlength = spc.geo.getLength();
            file << "Generated by Faunus -- https://github.com/mlund/faunus\n" << spc.numParticles() << "\n";
            for (auto &group : spc.groups) { // loop over groups
                for (auto &i : group) {      // loop over active particles
                    auto &name = Faunus::atoms.at(i.id).name;
                    Point pos = 0.1 * (i.pos + 0.5 * boxlength); // shift origin and convert to nm
                    file << fmt::format("{:5d}{:5}{:5}{:5d}{:8.3f}{:8.3f}{:8.3f}\n", nres, name, name, natom++, pos.x(),
                                        pos.y(), pos.z());
                }
                nres++;
            }
            file << 0.1 * boxlength.transpose() << "\n";
            return true;
        }
        return false;
    }
};

/**
 * @brief GROMACS xtc compressed trajectory file format
 *
 * Saves simulation frames to a Gromacs xtc trajectory file including
 * box information if applicable. Molecules with periodic boundaries
 * can be saved as "whole" by adding their groups to the public g-vector.
 *
 * @date June 2007-2011, Prague / Malmo
 * @note Alternative pure C++ version(?):
 * http://loos.sourceforge.net/xtc_8hpp_source.html
 */
class FormatXTC {
  private:
    XDRFILE *xd = nullptr; //!< file handle
    matrix xdbox;          //!< box dimensions
    rvec *x_xtc;           //!< vector of particle coordinates
    float time_xtc, prec_xtc = 1000;
    int natoms_xtc, step_xtc;

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
     *       You may want to transfer the new box size to the pair potential if
     *       periodic boundaries are used.
     */
    template <class Tspace> bool loadNextFrame(Tspace &c, bool setbox = true, bool applypbc = false) {
        if (xd != nullptr) {
            if (natoms_xtc == (int)c.p.size()) {
                int rc = read_xtc(xd, natoms_xtc, &step_xtc, &time_xtc, xdbox, x_xtc, &prec_xtc);
                if (rc == 0) {
                    // Geometry::Chameleon *geo = dynamic_cast<Geometry::Chameleon *>(&c.geo);
                    // if (geo == nullptr or geo->type not_eq Geometry::CUBOID)
                    //    throw std::runtime_error("Cuboid-like geometry required");
                    Point len_half = 0.5 * c.geo.getLength();
                    if (setbox)
                        c.geo.setLength(Point(10.0 * xdbox[0][0], 10.0 * xdbox[1][1], 10.0 * xdbox[2][2]));
                    for (size_t i = 0; i < c.p.size(); i++) {
                        c.p[i].pos.x() = 10.0 * x_xtc[i][0];
                        c.p[i].pos.y() = 10.0 * x_xtc[i][1];
                        c.p[i].pos.z() = 10.0 * x_xtc[i][2];
                        c.p[i].pos -= len_half;
                        if (applypbc)
                            c.geo.boundary(c.p[i].pos);
                        if (c.geo.collision(c.p[i].pos, 0))
                            throw std::runtime_error("particle-container collision");
                    }
                    return true;
                }
            } else
                throw std::runtime_error("xtcfile<->container particle mismatch");
        } else
            throw std::runtime_error("xtc file cannot be read");
        return false; // end of file or not opened
    }

    /**
     * This will take an arbitrary particle vector and add it
     * to an xtc file. If the file is already open, coordinates will
     * be added, while a new file is created if not.
     * Coordinates are shifted and converted to nanometers.
     * Box dimensions for the frame must be manually
     * set by the `setbox()` function before calling this.
     */
    template <class Titer1, class Titer2 /** particle vector iterator */>
    bool save(const std::string &file, Titer1 begin, Titer2 end) {
        if (begin != end) {
            if (xd == nullptr)
                xd = xdrfile_open(&file[0], "w");
            if (xd != nullptr) {
                rvec *x = new rvec[ranges::distance(begin, end)];
                size_t N = 0;
                for (auto j = begin; j != end; ++j) {
                    x[N][0] = j->pos.x() * 0.1 + xdbox[0][0] * 0.5; // AA->nm
                    x[N][1] = j->pos.y() * 0.1 + xdbox[1][1] * 0.5; // move inside sim. box
                    x[N][2] = j->pos.z() * 0.1 + xdbox[2][2] * 0.5; //
                    N++;
                }
                write_xtc(xd, N, step_xtc++, time_xtc++, xdbox, x, prec_xtc);
                delete[] x;
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
    void setLength(double);
    void setLength(const Point &);
};

std::vector<int> fastaToAtomIds(const std::string &); //!< Convert FASTA sequence to atom id sequence

/**
 * @brief Create particle vector from FASTA sequence with equally spaced atoms
 *
 * Particle positions is generated as a random walk
 */
std::vector<Particle> fastaToParticles(const std::string &, double = 7, const Point & = {0, 0, 0});

/**
 * @brief Load structure file into particle vector
 */
bool loadStructure(const std::string &, std::vector<Particle> &, bool, bool = true);

/**
 * @brief Create (compressed) output stream
 * @param filename Output filename
 * @param openmode Mode for opening file
 * @param compression Set to true for zlib compression
 * @return unique pointer to output stream
 */
std::unique_ptr<std::ostream> makeOutputStream(const std::string &, std::ios_base::openmode, bool);

/**
 * @brief Create (compressed) input stream from file
 * @param filename Input filename
 * @param openmode Mode for opening file
 * @return unique pointer to input stream
 * @note zlib compression is auto-detected
 */
std::unique_ptr<std::istream> makeInputStream(const std::string &, std::ios_base::openmode);

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
    FormatSpaceTrajectory(std::ostream &ostream) {
        if (ostream)
            output_archive = std::make_unique<cereal::BinaryOutputArchive>(ostream);
    }
    FormatSpaceTrajectory(std::istream &istream) {
        if (istream)
            input_archive = std::make_unique<cereal::BinaryInputArchive>(istream);
    }
    template <class Tspace> void load(Tspace &) {
      assert(input_archive != nullptr);
    } //!< Load single frame from stream
    template <class Tspace> void save(const Tspace &) {
      assert(output_archive != nullptr);
    } //!< Save single frame from stream
};

} // namespace Faunus
