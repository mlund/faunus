#ifndef FAU_IO_H
#define FAU_IO_H
#include <faunus/common.h>
#include <faunus/bonded.h>

#ifndef __cplusplus
#define __cplusplus
#endif
#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif
#include "xdrfile/xdrfile_trr.h"
#include "xdrfile/xdrfile_xtc.h"

namespace Faunus {
  
  class Group;
  class Molecular;

  /*! \brief Basic file I/O routines
   *  \author Mikael Lund
   */
  class io {
    private:
      string s;
    public:
      bool readfile(string, vector<string>&);       //!< Read entire file to vector
      bool writefile(string, string, std::ios_base::openmode=std::ios_base::out); //!< Save string to file
      void strip(vector<string> &, string="#");     //!< Remove lines containing pattern
      void splash(string);                          //!< Splashes the content of a file
  };

  /*! \brief Read/write AAM file format
   *  \author Mikael Lund
   */
  class ioaam {
    private:
      string p2s(particle &, int);
      particle s2p(string &);
      io fio;
    public:
      ioaam();
      p_vec p;
      bool load(string);
      bool save(string, p_vec&);
  };

  /*!
   * \brief PQR format
   * \date December 2007
   * \author Mikael Lund
   */
  class iopqr {
    private:
      io fio;
    public:
      iopqr();
      p_vec p;                   //!< Placeholder for loaded data
      bool save(string, p_vec&); //!< Save with particle charge
  };

  /*!
   * \brief Gromacs GRO format
   * \date December 2007
   * \author Mikael Lund
   * \todo Non cubic dimensions
   */
  class iogro {
    private:
      vector<string> v;
      particle s2p(string &);
      io fio;
    public:
      double len;            //!< Box side length (cubic so far)
      p_vec p;
      bool load(string);
      bool save(string, p_vec&);
  };

  /*! \brief GROMACS xtc compressed trajectory file format
   *  \author Mikael Lund
   *  \date June 2007-2011, Prague / Malmo
   *
   *  Saves simulation frames to a Gromacs xtc trajectory file including
   *  box information if applicable. Molecules with periodic boundaries
   *  can be saved as "whole" by adding their groups to the public g-vector.
   */
  class ioxtc {
    private:
      p_vec p; //!< internal particle vector for temporary data
      XDRFILE *xd;        //!< file handle
      matrix xdbox;       //!< box dimensions
      rvec *x_xtc;        //!< vector of particle coordinates
      float time_xtc, prec_xtc;
      int natoms_xtc, step_xtc;
    public:
      vector<Molecular*> g;                    //!< List of PBC groups to be saved as whole
      ioxtc(float);                            //!< Constructor that sets an initially cubic box
      bool open(string);                       //!< Open xtc file for reading
      bool loadnextframe(Space&);              //!< Load a single frame into cuboid
      bool save(string, const p_vec&);         //!< Save a frame to trj file.
      bool save(string, Space&);               //!< Save a frame to trj file (PBC)
      bool save(string, p_vec&, vector<Group>&);//!< Save groups
      void setbox(float);                      //!< Set box size to be saved in frame (cubic)
      void setbox(double,double,double);       //!< Set box size to be saved in frame
      void close();                            //!< Close trj file
  };

  /*! \brief Trajectory of charges per particle
   *  \author Chris Evers
   *  \date May 2011, Lund
   *
   *  Saves a trajectory of the charges for all particles in a particle vector
   */
  class ioqtraj {
    private:
      bool append;
      p_vec load(string);
    public:
      ioqtraj();
      bool save(string, p_vec&);   //!< Save a frame to trj file.
      bool save(string, p_vec&, vector<Group> &); //!< Save groups
  };
  
  class xyfile {
    private:
      io fio;
      std::ofstream f;
      unsigned int cnt;
    public:
      xyfile(string);
      void add(double, double);
      void close();
  };

  class FastaSequence {
    private:
      std::map<char,string> map; //!< Map one letter code (char) to three letter code (string)
      Energy::HarmonicBond bond;
    public:
      p_vec interpret(string);
      FastaSequence(double=0.76, double=4.9);
      Group insert(string, Space&, Energy::ParticleBonds&);
      Group include(string, Space&, Energy::ParticleBonds&);
  };
}//namespace
#endif
