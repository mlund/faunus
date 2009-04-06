#ifndef FAU_IO_H
#define FAU_IO_H
#include "faunus/common.h"
#include "faunus/group.h"
#include "faunus/container.h"
#include "faunus/titrate.h"


#ifdef GROMACS
#ifndef __cplusplus
#define __cplusplus
#endif
#include "xdrfile/xdrfile_trr.h"
#include "xdrfile/xdrfile_xtc.h"
#endif

//#include "xdrfile/xdrfile_trr.h"
//#endif GROMACS

namespace Faunus {
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
      void splash(string);                       //!< Splashes the content of a file
  };

  class particleIO {
    protected:
      string warning();
      string format;
    public:
      virtual bool save(container &, string)=0;
      virtual bool load(container &, string)=0;
      virtual bool load(vector<particle> &, string)=0;
      virtual string info() {
        return "";
      }
  };

  class io_aam : public particleIO {
    bool load(container &, string);
    bool load(vector<particle> &, string);
    bool sabe(container &, string);
  };

  //----------------------------------------
  /*!\brief General class for particle I/O.
   * \author Mikael Lund
   * \todo Unfinished!
   * \todo move transformations s2p p2s to iopart
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
    friend class iogro;
    friend class iopqr;
    private:
    atoms *atomPtr;
    vector<string> v; 
    vector<particle> p;
    public:
    iopart(atoms &a) {
      atomPtr=&a;
    }
    virtual vector<particle> load(string)=0;            //!< Load from disk
    virtual bool save(string, vector<particle>&)=0;     //!< Save to disk
  };

  //------------------------------------
  /*! \brief Write XYZ structure files, only 
   *  \brief intended for particles.p IO
   *  \author Mikael Lund
   */
  class ioxyz : public iopart {
    friend class particles;

    private:
    particle s2p(string &);
    particles *sys;

    public:
    ioxyz(atoms &);
    bool save(string, vector<particle>&);
    vector<particle> load(string);
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
      ioaam(atoms &);
      vector<particle> load(string);
      void load(container&,inputfile&,vector<macromolecule>&);//!< Read proteins from disk
      void loadlattice(container&,inputfile&,vector<macromolecule>&);//!< Read proteins from disk on to a lattice
      bool load(container &, string); //!< Reads a configuration from disk
      bool save(string, vector<particle>&);
  };

  //-----------------------------------------------
  /*! \brief Persistence of Vision Raytracer output
   *  \author Mikael Lund
   *  \todo Unfinished!
   */
  class iopov : public iopart {
    private:
      std::ostringstream o;
      string p2s(particle &, int=0);
      void header();
      vector<particle> load(string);
    public:
      iopov(container &);
      //iopov(container &, atoms &);
      void clear();                       //!< Clear output buffer
      void box(float);                    //!< Add cubic box
      void cell(float);                   //!< Add spherical cell
      void light(float);                  //!< Add light source
      void connect(point&, point&, float);//!< Connect two points w. a cylinder
      void camera();                      //!< Specify camera location and viewpoint
      bool save(string, vector<particle>&);
  };

  /*!
   * \brief PQR format
   * \date December 2007
   * \author Mikael Lund
   */
  class iopqr : public iopart {
    private:
      string p2s(particle &, int=0) { return string(); }
      void header() {}
      vector<particle> load(string) { return vector<particle>(); }
    public:
      iopqr(atoms &);
      bool save(string, vector<particle> &);            //!< Save with particle charge
      bool save(string, vector<particle> &, titrate &); //!< Save with average charges
  };


  /*!
   * \brief Gromacs GRO format
   * \date December 2007
   * \author Mikael Lund
   */
  class iogro : public iopart {
    private:
      string p2s(particle &, int=0) { return string(); }
      void header() {};
      vector<particle> load(string) { return vector<particle>(); }
      float len;
    public:
      iogro(atoms &, inputfile &);
      bool save(string, vector<particle> &);
      bool save(string, box &);
  };

#ifdef GROMACS
  //-----------------------------------------------
  /*! \brief GROMACS xtc compressed trajectory fileformat
   *  \author Mikael Lund
   *  \date June 2007, Prague
   *  \todo Filename ignored, should be changed. Static box length.
   *        The XTC format is now included in OpenBabel2!
   *  \note Distances are stored in nanometers.
   *
   *  This class is used for output of configurations
   *  to a GROMACS xtc file, originally designed for MD trajectories albeit
   *  with no forces included.
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
      rvec x[3300];
      XDFILE xd;
      float box[3][3], time, step;
    public:
      ioxtc(container::container &, float);
      bool save(string, vector<particle> &);
      void setbox(float);
      void close();
  };
#endif
};//namespace
#endif
