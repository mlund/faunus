#ifndef FAU_GROUP_H
#define FAU_GROUP_H

#include "faunus/common.h"
#include "faunus/average.h"
#include "faunus/particles.h"
#include "faunus/container.h"
#include "faunus/pot_spring.h"
#include "faunus/inputfile.h"

namespace Faunus {
  /*! \brief Groups set of particles, such as molecules, salt etc.
   *  \author mikaek lund
   *  \date Lund 2004
   *
   *  A group defines the range of a a protein, the mobile ions,
   *  a chain or whatever. If is used extensively in MC movements,
   *  energy calculations and many, many other places.
   */
  class group {
    protected:
      slump slp;
      string title;
    public:
      particle cm, cm_trial; 
      short int beg,end;    ///< Define range in particle vector. [beg;end]
      string name;          ///< Informative (arbitrary) name
      group(int=0);         ///< Constructor, initialize data.

      void set(short int,short int);        ///< Set particle range, "beg" and "end".
      short int size();                     ///< Number of particles in group
      virtual short int random();           ///< Picks a random particle within this group
      bool find(unsigned int) const;        ///< Check if particle is part of the group
      virtual double charge(vector<particle> &);//!< Calculate total charge
      point masscenter(vector<particle> &); //!< Calculate center-of-mass
      point masscenter(container &);        //!< Calc. center-of-mass
      virtual string info();                //!< Print information
      void operator+=(group);
      group operator+(group);

      void move(container &, point);                      //!< Translate group
      bool overlap(container &);                          //!< Test overlap w all particles
      void undo(particles &);
      void accept(particles &);                           //!< Accept a move
      void add(container &, vector<particle>, bool=false);//!< Add a particle vector
      void add(container &, particle::type, short);       //!< Add particles w. collision check
      short int count(vector<particle> &, particle::type);//!< Count number of specific particles
      virtual unsigned short displace(container&,double); //!< Displace random particle
      virtual void isobaricmove(container &, double){};   //!< Pressure scaling
      virtual unsigned short nummolecules();              //!< Number of molecules
      unsigned short numhydrophobic(vector<particle> &);  //!< Number of hydrophobic particles
  };

  /*!
   * \brief Group for handling salt particles
   */
  class salt : public group {
    private:
      short nanion, ncation;
    public:
      salt(particle::type=particle::NA,
          particle::type=particle::CL);
      particle::type anion, cation; //!< Anion and cation types
      double muex;                  //!< Excess chemical potential
      void add(container &, inputfile &);                 //!< Add salt as specified in config file
      string info(container &);     //!< Show info
      virtual void isobaricmove(container &, double);
  };

  class macromolecule : public group {
    public:
      macromolecule();
      point mu;            //!< Dipole moment
      average<float> Q;    //!< Total charge. Updated with space::charge()
      average<float> Q2;   //!< Total charge squared.
      average<float> dip;  //!< Dipole moment scalar.
      average<float> dip2; //!< Dipole moment scalar squared.

      string info();                           //!< Show info
      string info(container &);                //!< Show more info!
      void center(container &);                //!< Center group in origo (0,0,0)
      double charge(vector<particle> &);       //!< Calculate total charge
      double radius(vector<particle> &);       //!< Calculate radius
      double dipole(vector<particle> &);       //!< Calculate dipole moment
      void zmove(container &, double);         //!< Move in z-direction, only
      void rotate(container &, double);        //!< Rotate around a point
      void rotate(container &, point, double);
      void rotate(container &, point, point, double); //!< Rotate around arbitrary point
      using group::add;
      void add(container &, inputfile &);      //!< Add according to inputfile
      void operator=(group);                   //!< Copy from group
      virtual void isobaricmove(container &,double);//!< Displace CM with scale difference
      virtual unsigned short nummolecules();
  };

  /*! \brief Group container for an array of molecules
   *  \author mikaek lund
   *  \date Asljunga 2008
   *
   * This group derivative is for handling an atom collection of
   * molecules in the particle vector. A [] operator is
   * implemented for convenient access to individual molecules.
   * Note that the group::random() function is redefined to point to a
   * random molecule instead of an atom.\n
   * Example:\n
   * \code
   * int i = spc.random(); // pick random molecule
   * group w = spc[i];     // ..and return it as a group
   * w.beg; // --> first atom of i'th molecule
   * \endcode
   */
  class molecules : public group {
    private:
      group sel;                //!< A temporary group class
    public:
      molecules(unsigned short);//!< Constructor. Specify number of atoms in each molecule.
      unsigned short numatom;   //!< Number of atoms in each molecule
      short random();           //!< Pick a random molecule (NOT atom)
      string info();            //!< Show information
      group operator[](unsigned short); //!< Access n'th molecule
  };

  /*! \brief Group container for SPC/E water (and similar three-point water models)
   *  \author mikaek lund
   *  \date Asljunga 2008
   */
  class spc : public molecules {
    public:
      spc();
      string info();
      float dp_trans;           //!< Translational displacement
      float dp_rot;             //!< Rotational displacement
      average<float> dip;       //!< Average dipole moment
      average<float> dip2;      //!< Average squared dipole moment
  };

  /*! \brief Freely jointed chain with harmonic spring potentials
   *  \author mikaek lund
   */
  class chain : public group, private pot_spring {
    public:
      chain(container &, int, particle::type, double );
      chain(container &, int, particle::type, double , point &);
      double k;                   //!< Spring constant
      double req;                 //!< Equilibrium distance between monomers
      bool graftpoint;       //!< Is chain grafted? (default == false)
      point *GP;             //!< Pointer to graft point
      void add(container &, int, particle::type);
      void addgrafted(container &, int, particle::type, point &);
      double monomerenergy(vector<particle> &, short);    //!< Spring energy of a monomer
      double internalenergy(vector<particle> &);          //!< Internal spring energy
      //!< Spring potential
      inline double quadratic(point &p1, point &p2) {
        double r=p1.dist(p2)-req;
        return k*r*r;
      }
  };

  /*
  // A brush lying in the xy-plane

  class brush : public group, private pot_spring {
  private:
  vector<chain>;
  public:
  brush(int, int)
  };
  */

  class planarsurface : public group {
    public:
      point plane;                //!< Plane definition
      point offset;               //!< Plane translation
      planarsurface();
      void project(point &);      //!< Project point onto plane
  };

  /*! \brief This class is used for zwitter ions present on a surface.
   *  \author mikaek lund
   *
   *  One end of the zwitter ion (even particles) can move only on the
   *  surface while the other ends (odd particles) are allowed to move
   *  in any direction. The two ends are connected by a spring.
   */
  class zwittermembrane : public planarsurface, private pot_spring {
    public:
      short mate(short);                                  //!< Find Zwitter ion partners (two "mates")
      void add(container&, particle, particle, short=1);  //!< Add Zwitter ion(s)
      unsigned short displace(container&, double);        //!< Displace random particle
      double selfenergy(particles &);                     //!< Return internal spring energy
  };
}//namespace
#endif

