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
   *   name, number, x, y, z, charge number, weight, radius
   * - Positions and radii should be in angstroms
   * - Currently, data in the number field is ignored.
   * - No particular spacing is required.
   *
   *     # information
   *     # more information
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
            << a.charge << " " << a.mw << " " << a.radius << std::endl;
          return o.str();
        }

      template<class Tparticle>
        static Tparticle& s2p(const std::string &s, Tparticle &a) {
          std::stringstream o;
          string name, num;
          o << s;
          o >> name >> num >> a.x() >> a.y() >> a.z() >> a.charge >> a.mw >> a.radius;
          a.id = atom[name].id;
          a.hydrophobic = atom[a.id].hydrophobic;
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
   * @brief PQR format
   * @date December 2007
   *
   * Saves particles as a PQR file. This format is very simular
   * to PDB but also contains charges and radii of the proteins.
   */
  class FormatPQR {
    public:
      template<class Tpvec>
        static bool save(const string &file, const Tpvec &p) {
          int nres=1, natom=1;
          char buf[100];
          std::ostringstream o;
          // index, atom->name, atom->resname, atom->resid,x, y, z, atom->charge, atom->radius
          for (auto &p_i : p) {
            string name=atom[p_i.id].name;
            sprintf(buf, "ATOM  %5d %-4s %-4s%5d    %8.3f %8.3f %8.3f %.3f %.3f\n",
                natom++, name.c_str(), name.c_str(), nres,
                p_i.x(), p_i.y(), p_i.z(), p_i.charge, p_i.radius );
            o << buf;
            if ( atom[p_i.id].name=="CTR" ) nres++;
          }
          return IO::writeFile(file, o.str());
        }
  };

  /*!
   * \brief Gromacs GRO format
   * \date December 2007
   * \author Mikael Lund
   * \todo Non cubic dimensions
   */
  class FormatGRO {
    private:
      std::vector<string> v;
      particle s2p(string &);
    public:
      double len;            //!< Box side length (cubic so far)
      p_vec p;
      bool load(string);

      template<class Tparticle, class Talloc>
        bool save(const string &file, std::vector<Tparticle,Talloc> &p) {
          int nres=1, natom=1;
          char buf[79];
          double halflen=len/2;
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
            o << len << " " << len << " " << len << std::endl;
          return IO::writeFile(file, o.str());
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

  /*! \brief GROMACS xtc compressed trajectory file format
   *  \author Mikael Lund
   *  \date June 2007-2011, Prague / Malmo
   *
   *  Saves simulation frames to a Gromacs xtc trajectory file including
   *  box information if applicable. Molecules with periodic boundaries
   *  can be saved as "whole" by adding their groups to the public g-vector.
   */
  class FormatXTC {
    private:
      //p_vec p; //!< internal particle vector for temporary data
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
       * measure we do a container collision check to see if all particles are within the Cuboid
       * boundaries.
       *
       * \note The container particle vector *must* match the number of particles in the xtc file. If not
       *       an error message will be issued and the function will abort.
       * \note You may want to transfer the new box size to the pair potential if periodic boundaries are used.
       */
      template<class Tspace>
        bool loadnextframe(Tspace &c) {
          if (xd!=NULL) {
            if (natoms_xtc==(int)c.p.size()) { 
              int rc = read_xtc(xd, natoms_xtc, &step_xtc, &time_xtc, xdbox, x_xtc, &prec_xtc);
              if (rc==0) {
                static double ten=10;
                Geometry::Cuboid* geo = dynamic_cast<Geometry::Cuboid*>(c.geo);
                Point l( xdbox[0][0], xdbox[1][1], xdbox[2][2] );
                geo->setlen(l*ten);
                for (size_t i=0; i<c.p.size(); i++) {
                  c.p[i].x() = x_xtc[i][0]*ten - geo->len_half.x();     // store pos. in container.                 
                  c.p[i].y() = x_xtc[i][1]*ten - geo->len_half.y();
                  c.p[i].z() = x_xtc[i][2]*ten - geo->len_half.z();
                  c.trial[i].x()=c.p[i].x();
                  c.trial[i].y()=c.p[i].y();
                  c.trial[i].z()=c.p[i].z();
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

      /*!
       * This will take an arbitrary particle vector and add it
       * to an xtc file. No shifting is done - only modification is conversion
       * from aangstom to nanometers. The box dimensions for the frame must be manually
       * set by the ioxtc::setbox() function before calling this.
       */
      template<class Tpoint, class Talloc>
        bool save(const string &file, const std::vector<Tpoint,Talloc> &p) {
          if (xd==NULL)
            xd=xdrfile_open(&file[0], "w");
          if (xd!=NULL) {
            rvec *x = new rvec [p.size()];
            int i=0;
            for (auto &pi : p) {
              x[i][0] = (pi.x() ) * 0.1;      // AA->nm
              x[i][1] = (pi.y() ) * 0.1;
              x[i][2] = (pi.z() ) * 0.1;
              i++;
            }
            write_xtc(xd,p.size(),step_xtc++,time_xtc++,xdbox,x,prec_xtc);
            delete[] x;
            return true;
          }
          return false;
        }

      /*!
       * @brief Save a frame to trj file (PBC)
       *
       * Save all particles in Cuboid to xtc file. Molecules added to the ioxtc::g
       * vector will be made whole (periodic boundaries are temporarily undone). Box
       * dimensions are taken from the Cuboid class and the particles are shifted so
       * that origin is in the corner of the box (Gromacs practice)
       *
       * \param file Name of the output xtc file
       * \param c Cuboid container from which particles and box dimensions are read.
       */
      template<class T1, class T2>
        bool save(const string &file, const Space<T1,T2> &c) {
          Geometry::Cuboid* geo = dynamic_cast<Geometry::Cuboid*>(c.geo);
          assert(geo!=nullptr && "Only Cuboid geometries classes allowed.");
          if (geo==nullptr)
            return false;
          setbox(geo->len.x(), geo->len.y(), geo->len.z());
          std::vector<Point> p(c.p.size());
          for (auto gi : g) {
            gi->translate( c, -gi->cm );             // b.trial is moved to origo -> whole!
            for (auto j : *gi)
              p[j] = c.trial[j] + gi->cm;            // move back to cm without periodicity
            gi->undo(c);                             // restore to original PBC location
          }
          for (auto &pi : p)
            pi+=geo->len_half;                       // gromacs origo is in the corner of the box
          return save(file, p);                      // while in Cuboid we use the middle
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
      bool open(string s) {
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

      void close() {
        xdrfile_close(xd);
        xd=NULL;
        delete[] x_xtc;
      }

      FormatXTC(float len) {
        prec_xtc = 1000.;
        time_xtc=step_xtc=0;
        setbox(len);
        xd=NULL;
        x_xtc=NULL;
      }

      void setbox(double x, double y, double z) {
        assert(x>0 && y>0 && z>0);
        for (int i=0; i<3; i++)
          for (int j=0; j<3; j++)
            xdbox[i][j]=0;
        xdbox[0][0]=0.1*x; // corners of the
        xdbox[1][1]=0.1*y; // rectangular box
        xdbox[2][2]=0.1*z; // in nanometers!
      }

      void setbox(float len) { setbox(len,len,len); }

      void setbox(const Point &p) { setbox(p.x(), p.y(), p.z()); }

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
      FormatTopology();

      template<class Tspace>
        bool save(string topfile, Tspace &spc) {
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

  /*!
   * \brief File IO for faste protein sequences
   */
  template<class Tbonded>
    class FormatFastaSequence {
      private:
        std::map<char,string> map; //!< Map one letter code (char) to three letter code (string)
        Tbonded bond;
      public:
        FormatFastaSequence(double harmonic_k, double harmonic_req) : bond(harmonic_k, harmonic_req) {
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
        }

        p_vec interpret(string seq) {
          p_vec p;
          particle a;
          for (auto c : seq) {
            if (map.find(c)!=map.end() ) {
              a=atom[ map[c] ];
              p.push_back(a);
            }
          }
          return p;
        }

        /*!
         * Inserts at end of particle vector
         */
        template<class Tspace>
          Group insert(string fasta, Tspace &spc, Tbonded &b) {
            p_vec p = interpret(fasta);
            Group g( p.size() );
            if (p.size()>0) {
              for (auto &a : p)
                if ( spc.insert(a) )
                  g.resize( g.size()+1 );
              for (int i=g.front(); i<g.back(); i++)
                b.add(i, i+1, bond );
            }
            return g;
          }

        template<class Tspace>
          Group include(string file, Tspace &spc, Tbonded &b) {
            Group g;
            return g;
          }

    };
}//namespace
#endif
