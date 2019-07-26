#pragma once

#include "core.h"
#include "particle.h"
#include <range/v3/distance.hpp>

namespace Faunus {

#ifndef __cplusplus
#define __cplusplus
#endif
#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"

namespace IO {

/* @brief Read lines from file into vector */
bool readFile(const std::string &file, std::vector<std::string> &v);

/**
 * @brief Write std::string to file
 * @param file Filename
 * @param s String to write
 * @param mode `std::ios_base::out` (new, default) or `std::ios_base::app` (append)
 */
bool writeFile(const std::string &file, const std::string &s, std::ios_base::openmode mode = std::ios_base::out);

/**
 * @brief Strip lines matching a pattern
 * @param v vector of std::string
 * @param pat Pattern to search for
 */
void strip(std::vector<std::string> &v, const std::string &pat);

} // namespace IO

/**
 * @brief Read/write AAM file format
 * @todo Make "p" vector private.
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

    static std::string p2s(const Particle &a, int i) {
        std::ostringstream o;
        o.precision(5);
        auto &atom = Faunus::atoms.at(a.id);
        o << atom.name << " " << i + 1 << " " << a.pos.transpose() << " " << a.charge << " " << atom.mw << " "
          << atom.sigma / 2 << std::endl;
        return o.str();
    }

    // convert string line to particle
    static Particle &s2p(const std::string &s, Particle &a);

  public:
    typedef std::vector<Particle> Tpvec;
    static bool load(const std::string &file, Tpvec &target, bool _keepcharges = true);

    static bool save(const std::string &file, const Tpvec &pv);
};

/**
 * @brief PQR format
 * @date December 2007
 *
 * Saves particles as a PQR file. This format is very similar
 * to PDB but also contains charges and radii
 */
class FormatPQR {
  private:
    // Write box dimensions (standard PDB format)
    static std::string writeCryst1(const Point &len, const Point &angle = {90, 90, 90});

  public:
    typedef std::vector<Particle> Tpvec;

    /**
     * @brief Simple loader for PQR files
     *
     * This will load coordinates from a PQR file into
     * a particle vector. Atoms are identified by the
     * `ATOM` and `HETATM` keywords and the PQR atom
     * name is used to identify the particle from `AtomMap`.
     * If the unit cell dimension, given by `CRYST1`, is
     * found a non-zero vector is returned.
     *
     * @param file PQR file to load
     * @param p Target particle vector to *append* to
     * @return Vector with unit cell dimensions. Zero length if not found.
     * @note This is a lazy and unsafe loader that doesn't properly
     *       respect the PQR fixed column format. Has been tested to
     *       work with files from VMD, pdb2pqr, and Faunus.
     */
    static Point load(const std::string &file, Tpvec &p, bool keepcharges);

    /**
     * @brief Read trajectory. Each frame must be separated by the "END" keyword.
     * @param file File name
     * @param v Vector of particle vectors
     */
    static void load(const std::string &file, std::vector<Tpvec> &v);

    /**
     * @param file Filename
     * @param p Particle vector
     * @param len Unit cell dimensions (optional)
     * @param n Number of atoms in each residue (default: 1e9)
     */
    static bool save(const std::string &file, const Tpvec &p, Point len = Point(0, 0, 0), int n = 1e9);
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
    typedef std::vector<Particle> Tpvec;
    static bool save(const std::string &file, const Tpvec &p, const Point &box = {0, 0, 0});

    /*
     * @brief Load XYZ file with atom positions
     *
     * Loads coordinates from XYZ file into particle vector. Remaining atom properties
     * are taken from the defined atom list; a warning is issued of the atom name
     * is unknown.
     *
     * @param file Filename
     * @param p Destination particle vector
     * @param append True means append to `p` (default). If `false`,`p` is first cleared.
     */
    static bool load(const std::string &file, Tpvec &p, bool append = true);
};

/**
 * @brief MXYZ format
 * @date June 2013
 *
 * Saves particles as a modifiedXYZ file. This format has number of particles at the first line
 * comment on second line, which we use to have a box information, and this is followed by positions,
 * direction and patch direction on each line
 *
 * @todo Under construction
 */
class FormatMXYZ {
  private:
    static std::string p2s(const Particle &a, int);

    static Particle &s2p(const std::string &s, Particle &a);

  public:
    typedef std::vector<Particle> Tpvec;
    static bool save(const std::string &file, const Tpvec &p, const Point &len, int time);

  public:
    static bool load(const std::string &file, Tpvec &p, Point &len);
};

/**
 * @brief Gromacs GRO format
 * @date December 2007
 * @todo Non cubic dimensions
 */
class FormatGRO {
  private:
    std::vector<std::string> v;

    void s2p(const std::string &s, Particle &dst);

  public:
    typedef std::vector<Particle> Tpvec;
    double len; //!< Box side length (cubic so far)

    /**
     * @brief Load GRO file into particle vector
     * @param file Filename
     * @param p Destination particle vector
     */
    bool load(const std::string &file, Tpvec &p);

    bool save(const std::string &file, Tpvec &p, std::string mode = "");

    template <class Tspace> bool static save(const std::string &file, Tspace &spc) {
        typedef typename Tspace::Tparticle Tparticle;
        int nres = 1, natom = 1;
        char buf[79];
        Point halflen = spc.geo.getLength() * 0.5;
        std::ostringstream o;
        o << "Generated by Faunus -- http://faunus.sourceforge.net\n" << spc.p.size() << "\n";
        for (auto &g : spc.groups) {
            for (auto &i : g) {
                Point a = (i.pos + halflen) / 10; // angstron->nm
                std::string &name = atoms.at(i.id).name;
                sprintf(buf, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n", nres, name.c_str(), name.c_str(), natom++, a.x(), a.y(),
                        a.z());
                o << buf;
            }
            nres++;
        }
        o << spc.geo.getLength().transpose() / 10 << "\n";
        return IO::writeFile(file, o.str());
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
    template <class Tspace> bool loadnextframe(Tspace &c, bool setbox = true, bool applypbc = false) {
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
    bool open(std::string s);
    void close();

    FormatXTC(double len);
    ~FormatXTC();

    void setbox(double x, double y, double z);
    void setbox(double len);
    void setbox(const Point &p);
};

/**
 * @brief Convert FASTA sequence to atom id sequence
 * @param fasta FASTA sequence, capital letters.
 * @return vector of verified and existing atom id's
 */
inline auto fastaToAtomIds(const std::string &fasta) {
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
    names.reserve(fasta.size());

    for (auto c : fasta) {     // loop over letters
        auto it = map.find(c); // is it in map?
        if (it == map.end())
            throw std::runtime_error("Invalid FASTA letter '" + std::string(1, c) + "'");
        else
            names.push_back(it->second);
    }
    return Faunus::names2ids(atoms, names);
}

/**
 * @brief Create particle vector from FASTA sequence with equally spaced atoms
 *
 * Particle positions is generated as a random walk
 */
std::vector<Particle> fastaToParticles(const std::string &fasta, double spacing = 7, const Point &origin = {0, 0, 0});

/**
 * @brief Load structure file into particle vector
 * @param file filename to load (aam, pqr, xyz, ...)
 * @param dst destination particle vector
 * @param append if true, expand dst vector
 * @param keepcharges if true, ignore AtomData charges
 * @return true if successfully loaded
 */
bool loadStructure(const std::string &file, std::vector<Particle> &dst, bool append, bool keepcharges = true);

} // namespace Faunus
