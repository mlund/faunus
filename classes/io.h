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
    iopov(species &);
    void clear();                    //!< Clear output buffer
    void box(float);                 //!< Add cubic box
    void cell(float);                //!< Add spherical cell
    void light(float);               //!< Add light source
    void camera();                   //!< Specify camera location and viewpoint
    bool save(string, vector<particle>&);
};

#endif
