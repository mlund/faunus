#ifndef _io_h
#define _io_h
#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "species.h"
#include "group.h"
#include "container.h"

#ifdef GROMACS
#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif
#include "xtcio.h"
#endif

//--------------------------------
//!\brief Basic file I/O routines
//!\author Mikael Lund
class io {
  private:
    string s;
  public:
    bool readfile(string, vector<string>&);      //!< Read entire file to vector
    bool writefile(string, string);              //!< Save string to file
    void strip(vector<string> &, string="#");    //!< Remove lines containing pattern
};

//----------------------------------------
/*!\brief General class for particle I/O.
 * \author Mikael Lund
 * \todo Unfinished!
 *
 * The purpose of this class is to provide general read/write
 * routines for particle I/O. This can be used to load/save
 * structures, coordinate space, pdb files, povray etc.
 */
class iopart : private io {
  friend class ioaam;
  friend class iopov;
  friend class ioxyz;
  friend class ioxtz;
  private:
  species *spcPtr;
  vector<string> v; 
  vector<particle> p;
  public:
  iopart(species &spc) {spcPtr=&spc; }
  virtual vector<particle> load(string)=0;            //!< Load from disk
  virtual bool save(string, vector<particle>&)=0;     //!< Save to disk
};

//------------------------------------
/*! \brief Write XYZ structure files
 *  \author Mikael Lund
 */
class ioxyz : public iopart {
  private:
    vector<particle> load(string) {};
  public:
    ioxyz(species &);
    bool save(string, vector<particle>&);
};


//------------------------------------
/*! \brief Read/write AAM file format
 *  \author Mikael Lund
 *  \todo Unfinished!
 */
class ioaam : public iopart {
  private:
    string p2s(particle &, int=0); 
    particle s2p(string &); 
  public:
    ioaam(species &);
    vector<particle> load(string);
    bool save(string, vector<particle>&);
};

//-----------------------------------------------
/*! \brief Persistence of Vision Raytracer output
 *  \author Mikael Lund
 *  \todo Unfinished!
 */
class iopov : public iopart {
  private:
    ostringstream o;
    string p2s(particle &, int=0);
    void header();
    vector<particle> load(string) {}; 
  public:
//    iopov(species &);
    iopov(container &, species &);
    void clear();                       //!< Clear output buffer
    void box(float);                    //!< Add cubic box
    void cell(float);                   //!< Add spherical cell
    void light(float);                  //!< Add light source
    void connect(point&, point&, float);//!< Connect two points w. a cylinder
    void camera();                      //!< Specify camera location and viewpoint
    bool save(string, vector<particle>&);
};

#ifdef GROMACS
//-----------------------------------------------
/*! \brief GROMACS xtc compressed trajectory fileformat
 *  \author Mikael Lund
 *  \date June 2007, Prague
 *  \todo Filename ignored, should be changed. Static box length.
 *  \note Distances are stored in nanometers.
 *
 *  This class is used for output of configurations
 *  to a GROMACS xtc file, originally designed for MD trajectories.
 *  The MC configurations in the xtc file can subsequently be used
 *  in a number of other programs VMD, for example, as well as analysed
 *  using a range of tools as part of GROMACS -- distribution
 *  functions etc.
 *  It can also be used to store lenghty simulations as commonly
 *  done in MD.
 */
class ioxtc : public iopart {
  private:
    vector<particle> load(string) {}
    rvec x[1000];
    int xd;
    float box[3][3], time, step;
  public:
    ioxtc(container &);
    bool save(string, vector<particle> &);
    void close();
};
#endif

#endif
