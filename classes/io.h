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
   bool readfile(string, vector<string> &);     //!< Read entire file to vector
   bool writefile(string, string &);            //!< Save string to file
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
    vector<string> v;                           //!< All lines go here
    vector<particle> p;
    virtual unsigned char first()=0;            //!< Line w. first particle
    virtual unsigned char last()=0;             //!< Line w. last particle
    virtual particle s2p(string &)=0;           //!< String -> particle
    virtual string p2s(particle &,int=0)=0;     //!< Particle -> string
    virtual string header(vector<particle>&)=0; //!< Return header string
  public:
    iopart(species &);
    vector<particle> load(string);              //!< Load from disk
    bool save(vector<particle> &, string );     //!< Save to disk
};

//------------------------------------
/*! \brief Read/write AAM file format
 *  \author Mikael Lund
 *  \todo Unfinished!
 */
class ioaam : public iopart {
  private:
    unsigned char first();
    unsigned char last();
    string header(vector<particle>&);
    string p2s(particle &, int=0); 
    particle s2p(string &); 
  public:
    ioaam(species &);
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
    vector<particle> load(string); //!< Disable load routine
  public:
    iopov(species &);
    void box(float);                 //!< Add cubic box
    void cell(float);                //!< Add spherical cell
    void light(point);               //!< Add light source
    void camera(point, point);       //!< Specify camera location and viewpoint
    string header(vector<particle>&, int i=0);
};

#endif
