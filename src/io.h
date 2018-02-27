#pragma once

#include "core.h"
#include "geometry.h"

namespace Faunus {

//    namespace Xtc {
#ifndef __cplusplus
#define __cplusplus
#endif
#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"

//    }

    namespace IO {

        /* @brief Read lines from file into vector */
        inline bool readFile(const std::string &file, std::vector<std::string> &v) {
            std::ifstream f( file );
            if (f) {
                std::string s;
                while (getline(f,s))
                    v.push_back(s);
                f.close();
                return true;
            }
            std::cerr << "# WARNING! FILE " << file << " NOT READ!\n";
            return false;
        }

        /**
         * @brief Write std::string to file
         * @param file Filename
         * @param s String to write
         * @param mode `std::ios_base::out` (new, default) or `std::ios_base::app` (append)
         */
        inline bool writeFile(const std::string &file, const std::string &s,
                std::ios_base::openmode mode=std::ios_base::out) {
            std::ofstream f(file, mode);
            if (f) {
                f << s;
                return true;
            }
            return false;
        }

        /**
         * @brief Strip lines matching a pattern
         * @param v vector of std::string
         * @param pat Pattern to search for
         */
        inline void strip(std::vector<std::string> &v, const std::string &pat) {
            auto iter=v.begin();
            while (iter!=v.end())
                if ((*iter).find(pat) != std::string::npos)
                    v.erase(iter);
                else ++iter;
        }

    }//namespace

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

            template<class Tparticle>
                static std::string p2s(const Tparticle &a, int i) {
                    std::ostringstream o;
                    o.precision(5);
                    double radius = atoms<Tparticle>.at(a.id).sigma/2;
                    o << atoms<Tparticle>.at(a.id).name << " " << i+1 << " "
                        << a.pos.transpose() << " "
                        << a.charge << " " << a.mw << " " << radius << endl;
                    return o.str();
                }

            template<class Tparticle>
                static Tparticle& s2p(const std::string &s, Tparticle &a) {
                    std::stringstream o;
                    std::string name, num;
                    double radius;
                    o << s;
                    o >> name;
                    o >> num >> a.pos.x() >> a.pos.y() >> a.pos.z() >> a.charge >> a.mw >> radius;
                    if (a.id==0)
                        std::cerr << "Warning: Atom name " << name << " is not in the atom list.\n";
                    return a;
                }

        public:
            template<class Tpvec>
                static bool load(const std::string &file, Tpvec &target) {
                    std::vector<std::string> v;
                    target.clear();
                    if (IO::readFile(file,v)==true) {
                        IO::strip(v,"#");
                        unsigned int n=atoi(v[0].c_str());
                        target.resize(n);
                        for (unsigned int i=1; i<=n; i++)
                            s2p(v.at(i), target.at(i-1));
                        return true;
                    }
                    return false;
                }

            template<class Tpvec>
                static bool save(const std::string &file, const Tpvec &pv) {
                    std::ostringstream o;
                    o << pv.size() << endl;
                    for (size_t i=0; i<pv.size(); i++)
                        o << p2s(pv[i], i);
                    return IO::writeFile(file, o.str());
                }
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
            template<class Tvec>
                static std::string writeCryst1(const Tvec &len, Tvec angle=Tvec(90,90,90)) {
                    char buf[500];
                    sprintf(buf, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
                            len.x(),len.y(),len.z(),angle.x(),angle.y(),angle.z());
                    return std::string(buf);
                }
        public:

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
            template<class Tparticle, class Talloc>
                static Point load(const std::string &file, std::vector<Tparticle,Talloc> &p) {
                    Point len(0,0,0);
                    std::ifstream in(file);
                    if (in) {
                        Tparticle a;
                        int iatom, ires;
                        std::string line,key,aname,rname;
                        while (std::getline(in, line)) {
                            std::stringstream o(line);
                            while (o >> key)
                                if (key=="ATOM" || key=="HETATM") {
                                    double radius;
                                    o >> iatom >> aname;
                                    auto it = findName( atoms<Tparticle>, aname );
                                    if (it==atoms<Tparticle>.end())
                                        throw std::runtime_error("PQR load error: unknown atom name '" + aname + "'.");
                                    a = it->p;
                                    o >> rname >> ires
                                        >> a.pos.x() >> a.pos.y() >> a.pos.z() >> a.charge >> radius;
                                    p.push_back(a);
                                } else if (key=="CRYST1")
                                    o >> len.x() >> len.y() >> len.z();
                        }
                    }
                    return len;
                }

            /**
             * @brief Read trajectory. Each frame must be separated by the "END" keyword.
             * @param file File name
             * @param v Vector of particle vectors
             */
            template<class Tparticle, class Talloc>
                static void load( const std::string &file,  std::vector <std::vector<Tparticle, Talloc>> &v )
                {
                    std::ifstream in(file);
                    if ( in )
                    {
                        Tparticle a;
                        std::vector<Tparticle, Talloc> p;
                        int iatom, ires;
                        std::string line, key, aname, rname;
                        while ( std::getline(in, line))
                        {
                            std::stringstream o(line);
                            while ( o >> key )
                                if ( key == "ATOM" || key == "HETATM" )
                                {
                                    double radius;
                                    o >> iatom >> aname;
                                    a = atoms<Tparticle>[aname];
                                    o >> rname >> ires >> a.x() >> a.y() >> a.z() >> a.charge >> radius;
                                    p.push_back(a);
                                    if ( a.id == 0 )
                                        std::cerr << "Warning: Atom name " << aname << " is not in the atom list.\n";
                                }
                                else if ( key == "END" )
                                {
                                    v.push_back(p);
                                    p.clear();
                                    p.reserve(v.back().size());
                                }
                        }
                    }
                    else
                        throw std::runtime_error("Error loading PQR trajectory " + file);
                }

            /**
             * @param file Filename
             * @param p Particle vector
             * @param len Unit cell dimensions (optional)
             * @param n Number of atoms in each residue (default: 1e9)
             */
            template<class Tpvec, class Tvec=Point>
                static bool save(const std::string &file, const Tpvec &p, Tvec len=Point(0,0,0), int n=1e9) {
                    typedef typename Tpvec::value_type Tparticle;
                    int nres=1, natom=1;
                    char buf[500];
                    std::ostringstream o;
                    if (len.norm()>1e-6)
                        o << writeCryst1(len);
                    for (auto &p_i : p) {
                        auto& prop = atoms<Tparticle>.at(p_i.id);
                        std::string name=prop.name;
                        double radius = prop.sigma/2;
                        sprintf(buf, "ATOM  %5d %-4s %-4s%5d    %8.3f %8.3f %8.3f %.3f %.3f\n",
                                natom++, name.c_str(), name.c_str(), nres,
                                (p_i.pos+len/2).x(), (p_i.pos+len/2).y(), (p_i.pos+len/2).z(), p_i.charge, radius); // move particles inside the sim. box
                        o << buf;
                        if ( atoms<Tparticle>.at(p_i.id).name=="CTR" )
                            nres++;
                        else if (natom % n == 0)
                            nres++;
                    }
                    o << "END\n";
                    return IO::writeFile(file, o.str());
                }
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
        template<class Tpvec>
            static bool save(const std::string &file, const Tpvec &p, const Point &box={0,0,0}) {
                typedef typename Tpvec::value_type Tparticle;
                std::ostringstream o;
                o << p.size() << "\n" << box.transpose() << "\n";
                for (auto &i : p)
                    o << atoms<Tparticle>.at(i.id).name << " " << i.pos.transpose() << "\n";
                return IO::writeFile(file, o.str());
            }

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
        template<class Tparticle, class Talloc>
            static bool load(const std::string &file, std::vector<Tparticle,Talloc> &p, bool append=true) {
                std::ifstream f( file );
                if (f)
                {
                    if (append==false)
                        p.clear();
                    Tparticle a;
                    std::string comment, name;
                    size_t n;
                    f >> n;
                    p.reserve( p.size() + n );
                    std::getline(f, comment); // ">>" token doesn't gobble new line
                    std::getline(f, comment); // read comment line
                    for (size_t i=0; i<n; i++) {
                        f >> name;
                        auto it = findName( atoms<Tparticle>, name );
                        if (it==atoms<Tparticle>.end())
                            throw std::runtime_error("XYZ load error: unknown atom name '" + name + "'.");
                        a = it->p;
                        f >> a.pos.x() >> a.pos.y() >> a.pos.z();
                        p.push_back(a);
                    }
                    if (!p.empty())
                        return true;
                }
                return false;
            }
    };

    /**
     * @brief MXYZ format
     * @date June 2013
     *
     * Saves particles as a modifiedXYZ file. This format has number of particles at the first line
     * comment on second line, which we use to have a box information, and this is followed by positions,
     * direction and patch direction on each line
     */
    class FormatMXYZ {
        private:
            template<class Tparticle>
                static std::string p2s(const Tparticle &a, int i) {
                    std::ostringstream o;
                    o.precision(5);
                    o << a.transpose() << " " << a.dir.transpose() << " "
                        << a.patchdir.transpose() << "\n";
                    return o.str();
                }

            template<class Tparticle>
                static Tparticle& s2p(const std::string &s, Tparticle &a) {
                    std::stringstream o;
                    o << s;
                    o >> a.x() >> a.y() >> a.z() >> a.dir.x() >> a.dir.y() >> a.dir.z()
                        >> a.patchdir.x() >> a.patchdir.y() >> a.patchdir.z();
                    a.init();
                    return a;
                }
        public:
            template<class Tpvec, class Tvec=Point>
                static bool save(const std::string &file, const Tpvec &p, const Point &len, int time) {
                    std::ostringstream o;
                    o << p.size() << "\n"
                        << "sweep " << time << "; box " << len.transpose() << "\n";
                    for (size_t i=0; i<p.size(); i++)
                        o << p2s(p[i], i);
                    return IO::writeFile(file, o.str(), std::ios_base::app);
                }

        public:
            template<class Tpvec, class Tpoint>
                static bool load(const std::string &file, Tpvec &p, Tpoint &len) {
                    std::stringstream o;
                    std::vector<std::string> v;
                    if (IO::readFile(file,v)==true) {
                        IO::strip(v,"#");
                        size_t n=atoi(v[0].c_str());
                        if (p.size() != n)
                            std::cerr << "# mxyz load error: number of particles in xyz file " << n
                                << " does not match input file (" << p.size() << ")!" << endl;
                        o << v[1].erase(0, v[1].find_last_of("x")+1);
                        o >> len.x() >> len.y() >> len.z();
                        for (size_t i=2; i<n+2; i++)
                            s2p(v.at(i), p.at(i-2));
                        return true;
                    }
                    return false;
                }
    };

    /**
     * @brief Gromacs GRO format
     * @date December 2007
     * @todo Non cubic dimensions
     */
    class FormatGRO {
        private:
            std::vector<std::string> v;

            template<class Tparticle>
                inline void s2p(const std::string &s, Tparticle &dst) {
                    std::stringstream o;
                    std::string name;
                    o << s.substr(10,5) << s.substr(20,8) << s.substr(28,8) << s.substr(36,8);
                    o >> name >> dst.x >> dst.y >> dst.z;
                    dst=atoms<Tparticle>.at(name);
                    return 10*dst; //nm->angstrom
                }

        public:
            double len;            //!< Box side length (cubic so far)

            /**
             * @brief Load GRO file into particle vector
             * @param file Filename
             * @param p Destination particle vector
             */
            template<class Tparticle, class Talloc>
                bool load(const std::string &file, std::vector<Tparticle,Talloc> &p) {
                    p.clear();
                    v.resize(0);
                    if (IO::readFile(file,v)==true) {
                        int last=atoi(v[1].c_str())+1;
                        for (int i=2; i<=last; i++) {
                            Tparticle a;
                            s2p( v[i], a );
                            p.push_back( a );
                        }
                        return true;
                    }
                    return false;
                }

            template<class Tparticle, class Talloc>
                bool save(const std::string &file, std::vector<Tparticle,Talloc> &p, std::string mode="") {
                    int nres=1, natom=1;
                    char buf[79];
                    double halflen=len/20; // halflen in nm
                    std::ostringstream o;
                    o << "Generated by Faunus -- http://faunus.sourceforge.net"
                        << std::endl << p.size() << std::endl;
                    for (auto &pi : p) {
                        std::string name=atoms<Tparticle>.at(pi.id).name;
                        sprintf(buf, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
                                nres,name.c_str(),name.c_str(),natom++,
                                pi.x()/10+halflen, pi.y()/10+halflen, pi.z()/10+halflen );
                        o << buf;
                        if ( atoms<Tparticle>[pi.id].name=="CTR" )
                            nres++;
                    }
                    if (len>0)
                        o << len/10 << " " << len/10 << " " << len/10 << std::endl; // box side length in nm
                    if ( mode=="append" )
                        return IO::writeFile(file, o.str(), std::ios_base::app);
                    else
                        return IO::writeFile(file, o.str(), std::ios_base::out);
                }

            template<class Tspace>
                bool static save(const std::string &file, Tspace &spc) {
                    typedef typename Tspace::Tparticle Tparticle;
                    int nres=1, natom=1;
                    char buf[79];
                    Point halflen=spc.geo.getLength() * 0.5;
                    std::ostringstream o;
                    o << "Generated by Faunus -- http://faunus.sourceforge.net\n" << spc.p.size() << "\n";
                    for (auto &g : spc.groups) {
                        for (auto &i : g) {
                            Point a = (i.pos + halflen)/10; // angstron->nm
                            std::string &name = atoms<Tparticle>.at(i.id).name;
                            sprintf(buf, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
                                    nres,name.c_str(),name.c_str(),natom++, a.x(), a.y(), a.z() );
                            o << buf;
                        }
                        nres++;
                    }
                    o << spc.geo.getLength().transpose()/10 << "\n";
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
            XDRFILE *xd=NULL;        //!< file handle
            matrix xdbox;       //!< box dimensions
            rvec *x_xtc;        //!< vector of particle coordinates
            float time_xtc, prec_xtc=1000, prec_tol=1e-3;
            int natoms_xtc, step_xtc;
        public:
            inline int getNumAtoms() { return natoms_xtc; }

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
            template<class Tspace>
                bool loadnextframe(Tspace &c, bool setbox=true, bool applypbc=false) {
                    if (xd!=NULL)
                    {
                        if (natoms_xtc==(int)c.p.size())
                        { 
                            int rc = read_xtc(xd, natoms_xtc, &step_xtc, &time_xtc, xdbox, x_xtc, &prec_xtc);
                            if (rc==0)
                            {
                                Geometry::Cuboid* geo = dynamic_cast<Geometry::Cuboid*>(&c.geo);
                                if (geo==nullptr)
                                    throw std::runtime_error("Cuboid-like geometry required");
                                Point len_half = 0.5*geo->getLength();
                                if (setbox)
                                    geo->setLength( Point( 10.0*xdbox[0][0], 10.0*xdbox[1][1], 10.0*xdbox[2][2] ) );
                                for (size_t i=0; i<c.p.size(); i++) {
                                    c.p[i].pos.x() = 10.0*x_xtc[i][0];
                                    c.p[i].pos.y() = 10.0*x_xtc[i][1];
                                    c.p[i].pos.z() = 10.0*x_xtc[i][2];
                                    c.p[i].pos -= len_half;
                                    if (applypbc)
                                        geo->boundary( c.p[i].pos );
                                    if ( geo->collision(c.p[i].pos, 0) )
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
             * set by the `ioxtc::setbox()` function before calling this.
             */
            template<class Tgroup, class Tparticle, class Talloc>
                bool save(const std::string &file, const std::vector<Tparticle,Talloc> &p, Tgroup g) {
                    if (!g.empty()) {
                        if ( xd == NULL )
                            xd = xdrfile_open(&file[0], "w");
                        if ( xd != NULL ) {
                            rvec *x = new rvec[ g.size() ];
                            unsigned int i=0;
                            for ( auto &j : g ) {
                                x[i][0] = j.pos.x()*0.1 + xdbox[0][0]*0.5; // AA->nm
                                x[i][1] = j.pos.y()*0.1 + xdbox[1][1]*0.5; // move inside sim. box
                                x[i][2] = j.pos.z()*0.1 + xdbox[2][2]*0.5; //
                                i++;
                            }
                            write_xtc( xd, g.size(), step_xtc++, time_xtc++, xdbox, x, prec_xtc );
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
            inline bool open(std::string s) {
                if (xd!=NULL)
                    close();
                xd = xdrfile_open(&s[0], "r");
                if (xd!=NULL) {
                    int rc = read_xtc_natoms(&s[0], &natoms_xtc); // get number of atoms
                    if (rc==exdrOK) {
                        x_xtc = new rvec [natoms_xtc]; // resize coordinate array
                        return true;
                    }
                } else
                    std::cerr << "# ioxtc error: xtc file could not be opened." << endl;
                return false;
            }

            inline void close() {
                xdrfile_close(xd);
                xd=NULL;
                delete[] x_xtc;
            }

            FormatXTC(double len) {
                prec_xtc = 1000.;
                time_xtc=step_xtc=0;
                setbox(len);
                xd=NULL;
                x_xtc=NULL;
            }

            ~FormatXTC()
            {
                close();
            }

            inline void setbox(double x, double y, double z) {
                assert(x>0 && y>0 && z>0);
                for (int i=0; i<3; i++)
                    for (int j=0; j<3; j++)
                        xdbox[i][j]=0;
                xdbox[0][0]=0.1*x; // corners of the
                xdbox[1][1]=0.1*y; // rectangular box
                xdbox[2][2]=0.1*z; // in nanometers!
            }

            inline void setbox(double len) { setbox(len,len,len); }

            inline void setbox(const Point &p) { setbox(p.x(), p.y(), p.z()); }

    };

    template<class Tpvec, class Enable = void>
        struct loadStructure {
            bool operator()(const std::string &file, Tpvec &dst, bool append) {
                if (append==false)
                    dst.clear();
                std::string suffix = file.substr(file.find_last_of(".") + 1);
                if ( suffix == "xyz" )
                    FormatXYZ::load(file, dst);
                if ( !dst.empty() ) return true;
                return false;
            }
        }; //!< XYZ file into given particle vector (fallback if charges not available)

    template<class Tpvec>
        struct loadStructure<Tpvec, std::enable_if_t<std::is_base_of<Charge,typename Tpvec::value_type>::value>> {
            bool operator()(const std::string &file, Tpvec &dst, bool append)
            {
                if (append==false)
                    dst.clear();
                std::string suffix = file.substr(file.find_last_of(".") + 1);
                if ( suffix == "aam" )
                    FormatAAM::load(file, dst);
                if ( suffix == "pqr" )
                    FormatPQR::load(file, dst);
                if ( suffix == "xyz" )
                    FormatXYZ::load(file, dst);
                if ( !dst.empty() ) return true;
                return false;
            }
        }; //!< Load AAM/PQR/XYZ file into given particle vector


}//namespace
