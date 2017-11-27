#pragma once

#include "core.h"

namespace Faunus {

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
            std::cout << "# WARNING! FILE " << file << " NOT READ!\n";
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
            std::ofstream f(file.c_str(), mode);
            cout << "Writing to file '" << file << "'. ";
            if (f) {
                f << s;
                f.close();
                cout << "OK!\n";
                return true;
            }
            cout << "FAILED!\n";
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
                o << atoms<Tparticle>[a.id].name << " " << i+1 << " "
                << a.pos.transpose() << " "
                << a.charge << " " << a.mw << " " << a.radius << endl;
                return o.str();
            }

            template<class Tparticle>
            static Tparticle& s2p(const std::string &s, Tparticle &a) {
                std::stringstream o;
                std::string name, num;
                o << s;
                o >> name;
                //a = atoms<Tparticle>[name];
                o >> num >> a.pos.x() >> a.pos.y() >> a.pos.z() >> a.charge >> a.mw >> a.radius;
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
                                    o >> iatom >> aname;
                                    a = findName( atoms<Tparticle>, aname)->p;
                                    o >> rname >> ires
                                        >> a.pos.x() >> a.pos.y() >> a.pos.z() >> a.charge >> a.radius;
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
                                    o >> iatom >> aname;
                                    a = atoms<Tparticle>[aname];
                                    o >> rname >> ires >> a.x() >> a.y() >> a.z() >> a.charge >> a.radius;
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
                        std::string name=atoms<Tparticle>[p_i.id].name;
                        sprintf(buf, "ATOM  %5d %-4s %-4s%5d    %8.3f %8.3f %8.3f %.3f %.3f\n",
                                natom++, name.c_str(), name.c_str(), nres,
                                (p_i.pos+len/2).x(), (p_i.pos+len/2).y(), (p_i.pos+len/2).z(), p_i.charge, p_i.radius); // move particles inside the sim. box
                        o << buf;
                        if ( atoms<Tparticle>[p_i.id].name=="CTR" )
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
     * the first line; comment on second line followed by positions of all
     * particles xyz position on each line
     */
    struct FormatXYZ {
        template<class Tpvec, class Tvec=Point>
            static bool save(const std::string &file, const Tpvec &p, Tvec len=Tvec(0,0,0)) {
                typedef typename Tpvec::value_type Tparticle;
                std::ostringstream o;
                o << p.size() << "\nGenerated by Faunus\n";
                for (auto &i : p)
                    o << atoms<Tparticle>[i.id].name << " " << Tvec(i).transpose() << "\n";
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
                            throw std::runtime_error("FormatXYZ: unknown atom name '" + name + "'.");
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
                    dst=atoms<Tparticle>[name];
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
                        std::string name=atoms<Tparticle>[pi.id].name;
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
                bool save(const std::string &file, Tspace &spc) {
                    typedef typename Tspace::Tparticle Tparticle;
                    assert( std::fabs(len*len*len-spc.geo->getVolume())<1e-6
                            && "Did you forget to pass box size to FormatGRO?" );
                    std::string name;
                    int nres=1, natom=1;
                    char buf[79];
                    double halflen=len/2;
                    std::ostringstream o;
                    o << "Generated by Faunus -- http://faunus.sourceforge.net"
                        << std::endl << spc.p.size() << std::endl;
                    for (auto g : spc.groupList()) {
                        for (auto i : *g) {
                            name=atoms<Tparticle>[ spc.p[i].id ].name;
                            sprintf(buf, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
                                    nres,name.c_str(),name.c_str(),natom++,
                                    spc.p[i].x()/10+halflen, spc.p[i].y()/10+halflen, spc.p[i].z()/10+halflen );
                            o << buf;
                        }
                        nres++;
                    }
                    if (len>0)
                        o << len << " " << len << " " << len << std::endl;
                    return IO::writeFile(file, o.str());
                }
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
