#ifndef FAU_IO_H
#define FAU_IO_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/energy.h>

#ifndef __cplusplus
#define __cplusplus
#endif
#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif
#include "xdrfile/xdrfile_trr.h"
#include "xdrfile/xdrfile_xtc.h"

#endif

namespace Faunus {
  
  class Group;
  class GroupMolecular;

  /*
   * \brief Basic file I/O routines
   * \author Mikael Lund
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

  /*!
   * \brief Read/write AAM file format
   * \author Mikael Lund
   *
   * The AAM format is a simple format for loading particle positions, charges, radii and
   * molecular weights. The structure is as follows:
   * \li Lines beginning with # are ignored and can be placed anywhere
   * \li The first non-# line gives the number of particles
   * \li Every subsequent line gives atom information in the format: name, number, x, y, z, charge number, weight, radius
   * \li Positions and radii should be in angstroms
   * \li Currently, data in the number field is ignored.
   * \li No particular spacing is required.
   *
   * \code
   * # information
   * # more information
   * 2
   * Na    1     10.234 5.4454 -2.345  +1    22.0   1.7
   * Cl    2    5.011     1.054  20.02   -1   35.0   2.0
   * \endcode
   */
  class FormatAAM {
    private:
      string p2s(particle &, int);
      particle s2p(string &);
      io fio;
    public:
      FormatAAM();
      p_vec p;
      bool load(string);
      bool save(string, p_vec&);
  };

  /*!
   * \brief PQR format
   * \date December 2007
   * \author Mikael Lund
   */
  class FormatPQR {
    private:
      io fio;
    public:
      FormatPQR();
      p_vec p;                   //!< Placeholder for loaded data
      bool save(string, p_vec&); //!< Save with particle charge
  };

  /*!
   * \brief Gromacs GRO format
   * \date December 2007
   * \author Mikael Lund
   * \todo Non cubic dimensions
   */
  class FormatGRO {
    private:
      vector<string> v;
      particle s2p(string &);
      io fio;
    public:
      double len;            //!< Box side length (cubic so far)
      p_vec p;
      bool load(string);
      bool save(string, p_vec&);
      bool save(string, Space&);
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
      p_vec p; //!< internal particle vector for temporary data
      XDRFILE *xd;        //!< file handle
      matrix xdbox;       //!< box dimensions
      rvec *x_xtc;        //!< vector of particle coordinates
      float time_xtc, prec_xtc;
      int natoms_xtc, step_xtc;
    public:
      vector<GroupMolecular*> g;                    //!< List of PBC groups to be saved as whole
      FormatXTC(float);                            //!< Constructor that sets an initially cubic box
      bool open(string);                       //!< Open xtc file for reading
      bool loadnextframe(Space&);              //!< Load a single frame into cuboid
      bool save(string, const p_vec&);         //!< Save a frame to trj file.
      bool save(string, Space&);               //!< Save a frame to trj file (PBC)
      bool save(string, p_vec&, vector<Group>&);//!< Save groups
      void setbox(float);                      //!< Set box length - cubic
      void setbox(double,double,double);       //!< Set box length - xyz
      void setbox(const Point&);               //!< Set box length - xyz from vector
      void close();                            //!< Close trj file
  };

  class FormatTopology {
    private:
      int rescnt;
      string writeAtomTypes(const Space&);
      string writeMoleculeType(const Group&, const Space&);
    public:
      FormatTopology();
      bool save(string, Space&); //!< Generate topology from Space
  };

  /*! \brief Trajectory of charges per particle
   *  \author Chris Evers
   *  \date May 2011, Lund
   *
   *  Saves a trajectory of the charges for all particles in a particle vector
   */
  class FormatQtraj {
    private:
      bool append;
      p_vec load(string);
    public:
      FormatQtraj();
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

  /*!
   * \brief File IO for faste protein sequences
   */
  class FormatFastaSequence {
    private:
      std::map<char,string> map; //!< Map one letter code (char) to three letter code (string)
      Potential::Harmonic bond;
    public:
      p_vec interpret(string);
      FormatFastaSequence(double=0.76, double=4.9);
      Group insert(string, Space&, Energy::Bonded&);
      Group include(string, Space&, Energy::Bonded&);
  };
}//namespace
#endif
