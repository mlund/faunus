#ifndef FAU_IO_H
#define FAU_IO_H
#include "faunus/common.h"
#include "faunus/titrate.h"

#ifndef __cplusplus
#define __cplusplus
#endif
#ifndef CPLUSPLUS
#define CPLUSPLUS
#endif
#include "xdrfile/xdrfile_trr.h"
#include "xdrfile/xdrfile_xtc.h"

namespace Faunus {
  
  class container;
  class inputfile;
  class box;
  class group;
  class macromolecule;
  class particle;
  class particles;

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
    vector<string> v; 
    vector<particle> p;
    public:
    virtual vector<particle> load(string)=0;            //!< Load from disk
    virtual bool save(string, vector<particle>&)=0;     //!< Save to disk
  };

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
    ioxyz();
    bool save(string, vector<particle>&);
    vector<particle> load(string);
  };

  /*! \brief Read/write AAM file format
   *  \author Mikael Lund
   */
  class ioaam : public iopart {
    private:
      string p2s(particle &, int=0); 
      particle s2p(string &); 
    public:
      ioaam();
      vector<particle> load(string);
      void load(container&,inputfile&,vector<macromolecule>&);//!< Read proteins from disk
      void loadlattice(container&,inputfile&,vector<macromolecule>&);//!< Read proteins from disk on to a lattice
      bool load(container &, string); //!< Reads a configuration from disk
      bool save(string, vector<particle>&);
  };

  /*! \brief Persistence of Vision Raytracer output
   *  \author Mikael Lund
   */
  class iopov : public iopart {
    private:
      std::ostringstream o;
      string p2s(particle &, int=0);
      void header();
      vector<particle> load(string);
    public:
      iopov(container &);
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
      iopqr();
      bool save(string, vector<particle> &);            //!< Save with particle charge
      bool save(string, vector<particle> &, titrate &); //!< Save with average charges
      bool save(string, vector<particle> &, vector<group> &); //!< Save groups
  };

  /*!
   * \brief Gromacs GRO format
   * \date December 2007
   * \author Mikael Lund
   */
  class iogro : public iopart {
    private:
      particle s2p(string &);
      string p2s(particle &, int=0) { return string(); }
      void header() {};
      float len;
    public:
      iogro(inputfile &);
      bool save(string, vector<particle> &);
      bool save(string, box &);
      vector<particle> load(string);
  };

  /*! \brief GROMACS xtc compressed trajectory fileformat
   *  \author Mikael Lund
   *  \date June 2007, Prague
   *
   *  Saves simulation frames to a Gromacs xtc trajectory file including
   *  box information if applicable. Molecules with periodic boundaries
   *  can be saves as "whole" by adding their groups to the public g-vector
   *  when saving with save(string,box).
   */
  class ioxtc : public iopart {
    private:
      vector<particle> p;
      vector<particle> load(string);
      XDRFILE *xd;
      matrix xdbox;
      float time, step, prec_xtc;
    public:
      vector<group*> g;                        //!< List of PBC groups to be saved as whole
      ioxtc(float);
      bool OpenTrajectory(string);             //!< Not finished!
      bool LoadFrame(int, vector<particle> &); //!< Not finished!
      bool save(string, vector<particle> &);   //!< Save a frame to trj file.
      bool save(string, box &);                //!< Save a frame to trj file (PBC)
      bool save(string, vector<particle> &, vector<group> &); //!< Save groups
      void setbox(float);                      //!< Set box size to be saved in frame (cubic)
      void setbox(double,double,double);       //!< Set box size to be saved in frame
      void close();                            //!< Close trj file
  };
  
  class xyfile {
    private:
      std::ofstream f;
      unsigned int cnt;
    public:
      xyfile(string);
      void add(double, double);
      void close();
  };

};//namespace
#endif
