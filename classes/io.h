#ifndef _io_h
#define _io_h
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "species.h"
#include "group.h"

/*! \brief File I/O for structures, coordinates etc.
 *  \author Mikael Lund
 *  \todo Could be made more elegant using an "loadfile" class inherited
 *        by classes for the AAM format, PDB format etc.
 */
class io {
 private:
  struct filefmt {
    double x,y,z,mw,charge;
    string atomname, aminoname;
    unsigned int atomnr, aminonr;
  };

 public:
  vector<particle> loadaam(species &, string);  //!< Load structure in AAM format
  bool saveaam(species &, string, vector<particle> &, group &); //!< Save structure in AAM format
};

/*!
 * \brief General class for particle I/O.
 * \author Mikael Lund
 * \todo Unfinished!
 * \warning This class will eventually replace \link io \endlink
 *
 * The purpose of this class is to provide general read/write
 * routines for particle I/O. This can be used to load/save
 * structures, coordinate space, pdb files, povray etc.
 */
class iofile {
  friend class ioaam;
  private:
    species *spcPtr;
    vector<string> v;
    virtual particle s2p(string &s) {};
    virtual string p2s(particle &p) {};
    bool readfile(string);                      //!< Read entire file
  public:
    iofile(species &);
    vector<particle> load(string);              //!< Load from disk
    bool save(vector<particle> &, string );     //!< Save to disk
};

/*! \brief Read/write AAM file format
 *  \author Mikael Lund
 *  \toto Unfinished!
 */
class ioaam : public iofile {
  private:
    string p2s(particle &); //!< string->particle
    particle s2p(string &); //!< particle->string
  public:
    ioaam(species &);
};

/*! \brief Persistence of Vision Raytracer output
 *  \author Mikael Lund
 *  \todo Unfinished
 */
class iopov : public iofile {
  private:
    string p2s(particle &); //!< string->particle
    particle s2p(string &); //!< particle->string
    vector<particle> load(string); //!< Disable load routine
  public:
    iopov(species &);
};

#endif
