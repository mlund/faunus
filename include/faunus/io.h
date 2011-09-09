#ifndef FAU_IO_H
#define FAU_IO_H
#include <faunus/common.h>
#include <faunus/point.h>

#ifndef __cplusplus
#define __cplusplus
#endif
#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif
#include "xdrfile/xdrfile_trr.h"
#include "xdrfile/xdrfile_xtc.h"

namespace Faunus {
  
  class group;
  class molecular;
  class space;

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
      vector<particle> p;
      bool load(string);
      bool save(string, vector<particle>&);
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
      vector<particle> p;                   //!< Placeholder for loaded data
      bool save(string, vector<particle> &); //!< Save with particle charge
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
      vector<particle> p;
      bool load(string);
      bool save(string, vector<particle>&);
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
      vector<particle> p; //!< internal particle vector for temporary data
      XDRFILE *xd;        //!< file handle
      matrix xdbox;       //!< box dimensions
      rvec *x_xtc;        //!< vector of particle coordinates
      float time_xtc, prec_xtc;
      int natoms_xtc, step_xtc;
    public:
      vector<molecular*> g;                    //!< List of PBC groups to be saved as whole
      ioxtc(float);                            //!< Constructor that sets an initially cubic box
      bool open(string);                       //!< Open xtc file for reading
      bool loadnextframe(space&);              //!< Load a single frame into cuboid
      bool save(string, const vector<particle> &);//!< Save a frame to trj file.
      bool save(string, space&);               //!< Save a frame to trj file (PBC)
      bool save(string, vector<particle> &, vector<group> &); //!< Save groups
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
      vector<particle> load(string);
    public:
      ioqtraj();
      bool save(string, vector<particle> &);   //!< Save a frame to trj file.
      bool save(string, vector<particle> &, vector<group> &); //!< Save groups
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
}//namespace
#endif
