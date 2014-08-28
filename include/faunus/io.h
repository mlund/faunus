#ifndef FAU_IO_H
#define FAU_IO_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/geometry.h>
#include <faunus/group.h>
#include <faunus/space.h>

#ifndef __cplusplus
#define __cplusplus
#endif
#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif
#include <xdrfile/xdrfile_trr.h>
#include <xdrfile/xdrfile_xtc.h>

#endif
#include <xdrfile/xdrfile_trr.h>
#include <xdrfile/xdrfile_xtc.h>


namespace Faunus {

  namespace IO {

    /* @brief Read lines from file into vector */
    inline bool readFile(const std::string &file, std::vector<string> &v) {
      std::ifstream f(file.c_str() );
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
     * @brief Write string to file
     * @param file Filename
     * @param s String to write
     * @param mode `std::ios_base::out` (new, default) or `std::ios_base::app` (append)
     */
    inline bool writeFile(string file, string s,
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
     * @param v vector of string
     * @param pat Pattern to search for
     */
    inline void strip(std::vector<string> &v, const string &pat) {
      auto iter=v.begin();
      while (iter!=v.end())
        if ((*iter).find(pat)!=string::npos)
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
          o << atom[a.id].name << " " << i+1 << " " << a.x() << " " << a.y() <<" "<< a.z() << " "
            << a.charge << " " << a.mw << " " << a.radius << endl;
          return o.str();
        }

      template<class Tparticle>
        static Tparticle& s2p(const std::string &s, Tparticle &a) {
          std::stringstream o;
          string name, num;
          o << s;
          o >> name;
          a = atom[name];
          o >> num >> a.x() >> a.y() >> a.z() >> a.charge >> a.mw >> a.radius;
          return a;
        }

    public:
      template<class Tpvec>
        static bool load(const string &file, Tpvec &target) {
          vector<string> v;
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
        static bool save(const string &file, const Tpvec &pv) {
          std::ostringstream o;
          o << pv.size() << endl;
          for (size_t i=0; i<pv.size(); i++)
            o << p2s(pv[i], i);
          return IO::writeFile(file, o.str());
        }
  };

  /**
   *@brief Dipolar particles in VRML format
   */

  class DipoleWRL {
    public:
      template<class Tspace>
        void saveDipoleWRL(string filename, Tspace &spc, Group &sol) {
          auto L = spc.geo.len.x();
          auto L2 = L/2;

          std::ofstream f;
          f.open(filename);
          f << "#VRML V2.0 utf8\n" << endl;
          f << "Viewpoint {description" << "\"View 1\"" << "position" << "  "  << L*2 
            << "  " << "0.00    0.00" << endl
            << "orientation 0 1 0 0 }"<< endl
            << " Shape {appearance Appearance { material Material { " << endl
            << "diffuseColor 0.9 0.9 0.9 transparency 1.0 }}" << endl
            << "geometry Box { size"<< "  " << L << "  " << L << "  " << L << "}}" << endl
            << "Shape {appearance Appearance {material Material {" << endl
            << "emissiveColor 0.8 1.0 0.6 }}" << endl
            << "geometry IndexedLineSet {coord Coordinate {" << endl 
            << "point [" << endl

            << -L2 << "  " <<   L2  << "  " <<   L2 << endl 
            <<  L2 << "  " <<   L2  << "  " <<   L2 << endl 
            <<  L2 << "  " <<   L2  << "  " <<  -L2 << endl 
            << -L2 << "  " <<   L2  << "  " <<  -L2 << endl 
            << -L2 << "  " <<  -L2  << "  " <<   L2 << endl 
            <<  L2 << "  " <<  -L2  << "  " <<   L2 << endl 
            <<  L2 << "  " <<  -L2  << "  " <<  -L2 << endl 
            << -L2 << "  " <<  -L2  << "  " <<  -L2 << endl 
            << "]}" << endl 
            << "coordIndex [" << endl 
            << "0, 1, 2, 3, 0, -1," << endl 
            << "4, 5, 6, 7, 4, -1," << endl 
            << "0, 4, -1," << endl 
            << "1, 5, -1," << endl 
            << "2, 6, -1," << endl 
            << "3, 7, -1," << endl 
            << "]}}" << endl;



          for (auto i : sol) { 
            auto  x = Point(spc.p[i]).transpose().x(); 
            auto  y = Point(spc.p[i]).transpose().y(); 
            auto  z = Point(spc.p[i]).transpose().z(); 
            auto dipx = spc.p[i].mu.transpose().x();
            auto dipy = spc.p[i].mu.transpose().y();
            auto dipz = spc.p[i].mu.transpose().z();
            auto r = spc.p[i].radius;
            f << "Transform { translation" << "  " <<  x << "  " << y << "  "<< z << endl
              << " children [ DEF par_0 Shape{ appearance Appearance { material" << endl
              << " Material {diffuseColor       0.00      1.00      1.00 transparency 0.2 }}" << endl 
              << "  geometry Sphere {radius" <<"  " << r <<"}}]}" << endl;

            auto size = sqrt((dipx * dipx) + (dipy*dipy) + (dipz*dipz)); 
            auto cosT = dipy/size;
            auto angle = acos(cosT);

            auto Xdipx = (-1.0) * dipx;

            f << "Transform { translation" << "  " <<  x << "  " << y << "  "<< z << endl
              <<" rotation" << "  "<< dipz << "  " << " 0" << "  "<< Xdipx << " " <<  angle << endl 

              <<  " children[DEF tip_0 Shape{appearance Appearance{material Material{" << endl 

              <<  "diffuseColor       1.00      0.00      0.00 }}  geometry Cone { bottomRadius"<<"  "<< 0.3*r<<"  "<< "height" << "  "<< 0.6*r << "}}]}" << endl;
          }
          f.close();
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
        static string writeCryst1(const Tvec &len, Tvec angle=Tvec(90,90,90)) {
          char buf[100];
          sprintf(buf, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n", 
              len.x(),len.y(),len.z(),angle.x(),angle.y(),angle.z());
          return string(buf);
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
      template<class Tpvec>
        static Point load(const std::string &file, Tpvec &p) {
          Point len(0,0,0);
          std::ifstream in(file);
          if (in) {
            typename Tpvec::value_type a;
            int iatom, ires;
            std::string line,key,aname,rname;
            while (std::getline(in, line)) {
              std::stringstream o(line);
              while (o >> key)
                if (key=="ATOM" || key=="HETATM") {
                  o >> iatom >> aname;
                  a=atom[aname];
                  o >> rname >> ires
                    >> a.x() >> a.y() >> a.z() >> a.charge >> a.radius; 
                  p.push_back(a);
                } else if (key=="CRYST1")
                  o >> len.x() >> len.y() >> len.z();
            }
          }
          return len;
        }

      /**
       * @param file Filename
       * @param p Particle vector
       * @param len Unit cell dimensions (optional)
       * @param n Number of atoms in each residue (default: 1e20)
       */
      template<class Tpvec, class Tvec=Point>
        static bool save(const string &file, const Tpvec &p, Tvec len=Point(0,0,0), unsigned int n=1e9) {
          unsigned int nres=1, natom=1;
          char buf[100];
          std::ostringstream o;
          if (len.norm()>1e-6)
            o << writeCryst1(len);
          for (auto &p_i : p) {
            string name=atom[p_i.id].name;
            sprintf(buf, "ATOM  %5d %-4s %-4s%5d    %8.3f %8.3f %8.3f %.3f %.3f\n",
                natom++, name.c_str(), name.c_str(), nres,
                (p_i+len/2).x(), (p_i+len/2).y(), (p_i+len/2).z(), p_i.charge, p_i.radius); // move particles inside the sim. box
            o << buf;
            if ( atom[p_i.id].name=="CTR" )
              nres++;
            else if (natom % n == 0)
              nres++;
          }
          o << "END\n";
          return IO::writeFile(file, o.str());
        }
      /*
         sprintf(sd,"%-6s%5s %4s%c%-4s%c%4s%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s\n",
         recordname, indexbuf, atomname, altlocchar, resnamebuf, chain[0], 
         residbuf, insertion[0], x, y, z, occ, beta, segnamebuf, elementsymbol);
         */

      template<class Tgeo, class Tparticle>
        static bool save(const string &file, Space<Tgeo,Tparticle> &spc) {
          Group* oldg=nullptr;
          unsigned int nres=0, natom=1, cnt=0;
          char buf[100];
          string name, resname;
          std::ostringstream o;
          o << writeCryst1(spc.geo.len);
          for (auto i : spc.p) {
            name=atom[i.id].name;
            Group* g = spc.findGroup(cnt++);
            if (g!=oldg) {
              oldg=g;
              nres++;
            }
            resname = (g->name.empty() ? name : g->name);
            resname.resize(3,'x'); 
            i+=spc.geo.len/2;
            sprintf(buf, "ATOM  %5d %-4s %-4s%5d    %8.3f %8.3f %8.3f %.3f %.3f\n",
                natom++, name.c_str(), resname.c_str(), nres, i.x(), i.y(), i.z(),
                i.charge, i.radius);
            o << buf;
          }
          o << "END\n";
          assert(cnt==spc.p.size());
          return IO::writeFile(file, o.str());
        }

  };

  /**
   * @brief XYZ format
   * @date June 2013
   *
   * Saves particles as a XYZ file. This format has number of particles at the first line
   * comment on second line followed by positions of all particles xyz position on each line
   */
  class FormatXYZ {
    public:
      template<class Tpvec, class Tvec=Point>
        static bool save(const string &file, const Tpvec &p, Tvec len=Tvec(0,0,0)) {
          std::ostringstream o;
          o << p.size() << "\nGenerated by Faunus\n";
          for (auto &i : p)
            o << atom[i.id].name << " " << Tvec(i).transpose() << "\n";
          return IO::writeFile(file, o.str());
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
          o << a.x() << " " << a.y() <<" "<< a.z() << " " << a.dir.x() << " " << a.dir.y() <<" "<< a.dir.z() << " "
            << a.patchdir.x() << " " << a.patchdir.y() <<" "<< a.patchdir.z() << " " << std::endl;
          return o.str();
        }

      template<class Tparticle>
        static Tparticle& s2p(const std::string &s, Tparticle &a) {
          std::stringstream o;
          o << s;
          o >> a.x() >> a.y() >> a.z() >> a.dir.x() >> a.dir.y() >> a.dir.z()
            >> a.patchdir.x() >> a.patchdir.y() >> a.patchdir.z();
          return a;
        }
    public:
      template<class Tpvec, class Tvec=Point>
        static bool save(const string &file, const Tpvec &p, const Point &len, const unsigned int time) {
          char buf[200];

          std::ostringstream o;
          o << p.size() << "\n";
          sprintf(buf, "sweep %d; box %f %f %f \n", time, len.x(),len.y(),len.z());
          o << buf;
          for (size_t i=0; i< p.size(); i++)
            o << p2s(p[i], i);
          return IO::writeFile(file, o.str(), std::ios_base::app);
        }

    public:
      template<class Tpvec>
        static bool load(const string &file, Tpvec &p, Point &len) {
          std::stringstream o;
          vector<string> v;

          if (IO::readFile(file,v)==true) {
            IO::strip(v,"#");
            unsigned int n=atoi(v[0].c_str());
            //target.resize(n);
            if (p.size() != n)
              std::cerr << "# mxyz load error: number of particles in xyz file " << n
                << " does not match input file (" << p.size() << ")!" << endl;
            o << v[1].erase(0, v[1].find_last_of("x")+1);
            o >> len.x() >> len.y() >> len.z();
            for (unsigned int i=2; i<n+2; i++)
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
      std::vector<string> v;

      inline particle s2p(const string &s) {
        std::stringstream o;
        string name;
        double x,y,z;
        o << s.substr(10,5) << s.substr(20,8) << s.substr(28,8) << s.substr(36,8);
        o >> name >> x >> y >> z;
        particle p;
        p=atom[name]; 
        return 10*p; //nm->angstrom
      }

    public:
      double len;            //!< Box side length (cubic so far)
      p_vec p;

      inline bool load(const string &file) {
        p.clear();
        v.resize(0);
        if (IO::readFile(file,v)==true) {
          int last=atoi(v[1].c_str())+1;
          for (int i=2; i<=last; i++)
            p.push_back( s2p(v[i]) );
          return true;
        }
        return false;
      }

      template<class Tparticle, class Talloc>
        bool save(const string &file, std::vector<Tparticle,Talloc> &p, string mode="") {
          int nres=1, natom=1;
          char buf[79];
          double halflen=len/20; // halflen in nm
          std::ostringstream o;
          o << "Generated by Faunus -- http://faunus.sourceforge.net"
            << std::endl << p.size() << std::endl;
          for (auto &pi : p) {
            string name=atom[pi.id].name;
            sprintf(buf, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
                nres,name.c_str(),name.c_str(),natom++,
                pi.x()/10+halflen, pi.y()/10+halflen, pi.z()/10+halflen );
            o << buf;
            if ( atom[pi.id].name=="CTR" )
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
        bool save(string file, Tspace &spc) {
          assert( abs(len*len*len-spc.geo->getVolume())<1e-6
              && "Did you forget to pass box size to FormatGRO?" );
          string name;
          int nres=1, natom=1;
          char buf[79];
          double halflen=len/2;
          std::ostringstream o;
          o << "Generated by Faunus -- http://faunus.sourceforge.net"
            << std::endl << spc.p.size() << std::endl;
          for (auto g : spc.groupList()) {
            for (auto i : *g) {
              name=atom[ spc.p[i].id ].name;
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
      XDRFILE *xd;        //!< file handle
      matrix xdbox;       //!< box dimensions
      rvec *x_xtc;        //!< vector of particle coordinates
      float time_xtc, prec_xtc;
      int natoms_xtc, step_xtc;
    public:
      std::vector<Group*> g;          //!< List of PBC groups to be saved as whole

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
        bool loadnextframe(Tspace &c) {
          if (xd!=NULL) {
            if (natoms_xtc==(int)c.p.size()) { 
              int rc = read_xtc(xd, natoms_xtc, &step_xtc, &time_xtc, xdbox, x_xtc, &prec_xtc);
              if (rc==0) {
                Geometry::Cuboid* geo = dynamic_cast<Geometry::Cuboid*>(c.geo);
                assert(geo!=nullptr && "Geometry must to derived from Cuboid");
                geo->setlen( Point( xdbox[0][0], xdbox[1][1], xdbox[2][2] ) );
                for (size_t i=0; i<c.p.size(); i++) {
                  c.p[i].x() = x_xtc[i][0];
                  c.p[i].y() = x_xtc[i][1];
                  c.p[i].z() = x_xtc[i][2];
                  c.p[i] = c.p[i]*10 - geo->len_half;
                  c.trial[i] = Point(c.p[i]);
                  if ( geo->collision(c.p[i]) ) {
                    std::cerr << "# ioxtc load error: particle-container collision!" << endl;
                    return false;
                  }
                } 
                return true;
              }
            } else
              std::cerr << "# ioxtc load error: xtcfile-container particle mismatch!" << endl;
          } else
            std::cerr << "# ioxtc load error: xtc file not available for reading!" << endl;
          return false;
        }

      /**
       * This will take an arbitrary particle vector and add it
       * to an xtc file. If the file is already open, coordinates will
       * be added, while a new file is created if not.
       * Coordinates are shiftet and converted to nanometers.
       * Box dimensions for the frame must be manually
       * set by the `ioxtc::setbox()` function before calling this.
       */
      template<class Tpoint, class Talloc>
        bool save(const string &file, const std::vector<Tpoint,Talloc> &p) {
          if (xd==NULL)
            xd=xdrfile_open(&file[0], "w");
          if (xd!=NULL) {
            rvec *x = new rvec [p.size()];
            unsigned int i=0;
            for (auto &pi : p) {
              x[i][0] = pi.x()*0.1 + xdbox[0][0]*0.5; // AA->nm
              x[i][1] = pi.y()*0.1 + xdbox[1][1]*0.5; // move inside sim. box
              x[i][2] = pi.z()*0.1 + xdbox[2][2]*0.5; //
              i++;
            }
            write_xtc(xd,p.size(),step_xtc++,time_xtc++,xdbox,x,prec_xtc);
            delete[] x;
            return true;
          }
          return false;
        }

      /*
       * @brief Save a frame to trj file (PBC)
       *
       * Save all particles in Cuboid to xtc file. Molecules added to the `ioxtc::g`
       * vector will be made whole (periodic boundaries are temporarily undone). Box
       * dimensions are taken from the Cuboid class and the particles are shifted so
       * that origin is in the corner of the box (Gromacs practice)
       *
       * @param file Name of the output xtc file
       * @param c Cuboid container from which particles and box dimensions are read.
       *
       * @warning This is broken!
       */
      template<class T1, class T2>
        bool save(const string &file, Space<T1,T2> &c) {
          Geometry::Cuboid* geo = dynamic_cast<Geometry::Cuboid*>(&c.geo);
          assert(geo!=nullptr && "Only Cuboid geometries classes allowed.");
          if (geo==nullptr)
            return false;
          setbox(geo->len.x(), geo->len.y(), geo->len.z());
          auto p=c.p;
          for (auto gi : g) {
            gi->translate( c, -gi->cm );  // b.trial is moved to origo -> whole!
            for (auto j : *gi)
              p[j] = c.trial[j] + gi->cm; // move back to cm without periodicity
            gi->undo(c);                  // restore to original PBC location
          }
          for (auto &pi : p)
            pi+=geo->len_half;            // gromacs origo is in the corner of the box
          return save(file, p);           // while in Cuboid we use the middle
        }

      template<class Tpvec>
        bool save(const string &file, Tpvec &p, std::vector<Group> &g) {
          p_vec t;
          for (auto &gi : g)
            for (auto j : gi)
              t.push_back( p[j] );
          return save(file, t);
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

  class FormatTopology {
    private:
      int rescnt;

      template<class Tspace>
        string writeAtomTypes(const Tspace &spc) {
          char w=8;
          std::ostringstream o;
          o.precision(6);
          o << endl << "[ atomtypes ]" << endl;
          std::set<particle::Tid> ids; // found particle id's in Space
          for (auto &p : spc.p)
            ids.insert( p.id );
          for (auto i : ids) {
            o << std::right << setw(w) << atom[i].name << setw(w) << atom[i].mw 
              << setw(w) << atom[i].charge << setw(w) << "A" << setw(w) << 2*atom[i].radius/10
              << setw(w) << atom[i].eps << endl;
          }
          return o.str();
        }

      template<class Tspace>
        string writeMoleculeType(const Group &g, const Tspace &spc) {
          if (g.empty() || g.isAtomic())
            return "";
          rescnt++;
          std::map<particle::Tid, unsigned short> idcnt;
          char w=10;
          std::ostringstream o;
          o.precision(6);
          o << endl << "[ moleculetype ]" << endl << g.name << "      1" << endl;
          o << endl << "[ atoms ]" << endl;
          o << ";   nr     type     resnr   residue   atom         cgnr      charge    mass" << endl;
          for (auto i : g) {
            auto id=spc.p[i].id;
            idcnt[id]++;
            o << setw(w-5) << i-g.front()+1 << setw(w) << atom[id].name << setw(w) << rescnt
              << setw(w) << g.name << setw(5) << atom[id].name
              << std::left << setw(3) << idcnt[id] << std::right << setw(w) << (i%15)
              << setw(w) << spc.p[i].charge << setw(w) << spc.p[i].mw << endl;
          }
          return o.str();
        }

    public:
      FormatTopology() : rescnt(0) {}

      template<class Tspace>
        bool save(const string &topfile, Tspace &spc) {
          std::set<string> done;
          std::ofstream f(topfile.c_str());
          if (f) {
            f << "[ defaults ]\n"
              << ";nbfunc     comb-rule      gen-pairs     fudgeLJ      fudgeQQ\n"
              << "   1             2            yes         1.0000      1.0000\n";
            f << writeAtomTypes(spc);
            for (auto g : spc.groupList() ) {
              if ( done.find(g->name)==done.end() ) {
                f << writeMoleculeType(*g, spc);
                done.insert(g->name);
              }
            }
            f << "[ system ]\n" << "something" << endl;
            f << "[ molecules ]\n" << "prot " << spc.groupList().size() << endl;
            return true;
          }
          return false;
        }
  };

  /** @brief Generates TCL script for adding bonds in VMD */
  template<class Tbondlist>
    std::string VMDBonds(const Tbondlist &bondlist) {
      std::ostringstream o;
      for (auto &b : bondlist)
        o << "topo addbond " << b.first.first << " " << b.first.second << "\n";
      return o.str();
    }  

  /** @brief Trajectory of charges per particle
   *  @author Chris Evers
   *  @date May 2011, Lund
   *
   *  Saves a trajectory of the charges for all particles in a particle vector
   */
  class FormatQtraj {
    private:
      bool append;
    public:
      FormatQtraj() : append(false) {};

      /** @brief Save groups */
      template<class Tpvec>
        bool save(const string &file, const Tpvec &p, vector<Group> &g) {
          decltype(p) t;
          for (size_t i=0; i<g.size(); i++)
            for (auto j : g[i])
              t.push_back( p[j] );
          return save(file, t);
        }

      /** @brief Save a frame to trj file */
      template<class Tpvec>
        bool save(const string &file, const Tpvec &p) {
          std::ostringstream o;
          o.precision(6);
          for (size_t i=0; i<p.size(); i++)
            o << p[i].charge << " ";
          o << endl;
          if ( append==true )
            return IO::writeFile(file, o.str(), std::ios_base::app);
          else
            append=true;
          return IO::writeFile(file, o.str(), std::ios_base::out);
        }
  };

  /*
   * @brief Add bonded peptide from fasta sequence
   * @param spc Space
   * @param bonded Bonded energy class, i.e. `Energy::Bonded`
   * @param pairpot Bond potential, i.e. `Potential::Harmonic`
   * @param fasta Fasta sequence (captital letters, no spacing)
   * @return Group with peptide -- remember to enroll in space using `Space::enroll`
   * @note Untested
   */
  template<class Tspace, class Tbonded, class Tpairpot>
    Group addFastaSequence(Tspace &spc, Tbonded &bonded, const Tpairpot &pairpot, const string &fasta) {
      std::map<char,string> map;
      map['A']="ALA";
      map['R']="ARG";
      map['N']="ASN";
      map['D']="ASP";
      map['C']="CYS";
      map['E']="GLU";
      map['Q']="GLN";
      map['G']="GLY";
      map['H']="HIS";
      map['I']="ILE";
      map['L']="LEU";
      map['K']="LYS";
      map['M']="MET";
      map['F']="PHE";
      map['P']="PRO";
      map['S']="SER";
      map['T']="THR";
      map['W']="TRP";
      map['Y']="TYR";
      map['V']="VAL";

      Group g;
      typename Tspace::ParticleVector p;
      typename Tspace::ParticleVector::value_type a;

      // interpret fasta sequence
      for (auto c : fasta) {
        if (map.find(c)!=map.end() ) {
          a=atom[ map[c] ];
          p.push_back(a);
        } else std::cerr << "Unknown character in fasta sequence!";
      }

      if (!p.empty()) {
        for (auto &i : p)
          spc.geo.randompos(i);        // random positions
        g = spc.insert(p);             // add to space
        g.name = fasta;
        cout << "FASTA group: " << g.info() << "\n";
        for (int i=g.front(); i<g.back(); i++)
          bonded.add(i, i+1, pairpot); // add bonds
      }
      return g;
    }

}//namespace
#endif
