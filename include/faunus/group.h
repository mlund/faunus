#ifndef FAU_GROUP_H
#define FAU_GROUP_H

#include "faunus/common.h"
#include "faunus/average.h"
#include "faunus/particles.h"
#include "faunus/container.h"
#include "faunus/pot_spring.h"
#include "faunus/inputfile.h"
#include "faunus/iobabel.h"

namespace Faunus {
  /*! \brief Groups set of particles, such as molecules, salt etc.
   *  \author Mikael Lund
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
      bool find(unsigned int) const;        ///< Check if particle is part of the group
      short int size() const;               ///< Number of particles in group
      virtual double charge(const vector<particle> &);//!< Calculate total charge
      virtual short int random();           ///< Picks a random particle within this group
      point masscenter(const vector<particle> &); //!< Calculate center-of-mass
      point masscenter(const container &);  //!< Calc. center-of-mass
      virtual string info();                //!< Print information
      bool operator==(const group&) const;
      group& operator+=(const group&);
      const group operator+(const group&) const;

      virtual void move(container &, point);              //!< Translate group
      void invert(vector<particle> &, point &);           //!< Invert a group
      bool overlap(container &);                          //!< Test overlap w all particles
      virtual void undo(particles &);
      virtual void accept(particles &);                   //!< Accept a move
      void add(container &, vector<particle>, bool=false);//!< Add a particle vector
      void add(container &, unsigned char, short);        //!< Add particles w. collision check
      bool swap(container &, group &);                    //!< Swap location of two groups
      short int count(vector<particle> &, unsigned char); //!< Count number of specific particles
      virtual unsigned short displace(container&, point); //!< Displace random particle
      virtual void isobaricmove(container &, double){};   //!< Pressure scaling
      virtual unsigned short nummolecules();              //!< Number of molecules
      unsigned short numhydrophobic(vector<particle> &);  //!< Number of hydrophobic particles
      bool swap(container &, int);                        //!< Move group to a new position
  };

  /*!
   * \brief Group for handling salt particles
   */
  class salt : public group {
    private:
      short nanion, ncation;
    public:
      using group::info;
      salt(unsigned char=0, unsigned char=0);
      unsigned char anion, cation;  //!< Anion and cation types
      double muex;                  //!< Excess chemical potential
      void add(container &, inputfile &);                 //!< Add salt as specified in config file
      string info(container &);     //!< Show info
      virtual void isobaricmove(container &, double);
  };

  class macromolecule : public group {
    public:
      macromolecule();
      point mu;            //!< Dipole moment
      float conc;
      average<double> Q;    //!< Total charge. Updated with space::charge()
      average<double> Q2;   //!< Total charge squared.
      average<double> dip;  //!< Dipole moment scalar.
      average<double> dip2; //!< Dipole moment scalar squared.

      string info();                              //!< Show info
      string info(container &);                   //!< Show more info!
      void center(container &);                   //!< Center group in origo (0,0,0)
      double charge(const vector<particle> &);    //!< Calculate total charge
      double getcharge(const vector<particle> &); //!< Calculate total charge
      double radius(vector<particle> &);          //!< Calculate radius
      double gradius(vector<particle> &);         //!< Calculate radius of gyration
      double vradius(vector<particle> &);         //!< Volume based protein radius
      double dipole(vector<particle> &);          //!< Calculate dipole moment
      void zmove(container &, double);            //!< Move in z-direction, only
      virtual void rotate(container &, double, double=0); //!< Rotate around a point
      void rotate(container &, point, double, double);
      void rotate(container &, point, point, double); //!< Rotate around arbitrary point
      void transrot(container &, double, double); 
      using group::add;
      void add(container &, inputfile &);         //!< Add according to inputfile
      macromolecule& operator=(const group&);     //!< Copy from group
      virtual void isobaricmove(container &,double);  //!< Displace CM with scale difference
      virtual unsigned short nummolecules();
  };

  /*! \brief Group container for an array of molecules
   *  \author Mikael Lund
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
  class molecules : public macromolecule {
    private:
      group sel;                //!< A temporary group class
    public:
      molecules(unsigned short);//!< Constructor. Specify number of atoms in each molecule.
      unsigned short numatom;   //!< Number of atoms in each molecule
      short random();           //!< Pick a random molecule (NOT atom)
      string info();            //!< Show information
      group operator[](unsigned short); //!< Access n'th molecule
      void add(container &, vector<particle> &, short=1);
      vector<int> pick(int );
  };

  /*! \brief Class for polymer molecules
   *  \author Mikael Lund
   *  \date Lund 2009
   */
  class polymer : public macromolecule {
    public:
      vector< vector<unsigned short> > nb; //!< Neighbor list
      polymer();
      vector<unsigned short> neighbors(unsigned short) const;
      string info();
      bool areneighbors(unsigned short, unsigned short) const;
#ifdef BABEL
      bool babeladd( container &, inputfile & ); //!< Load molecule from disk using OpenBabel
#endif
  };

  /*! \brief Class for phospholipid-membrane
   *  \author Bjoern Persson
   *  \date Lund 2009
   *
   */
  class popscmembrane : public group{
    public:
      double scratio;          //Ratio of pops in percent
      double headarea;         //Headgroup area
      vector<polymer> pops;    //Pops vector
      vector<polymer> popc;    //Popc vector
      popscmembrane();
      void load(inputfile &, slit &);  //Scannes for "scratio", "headarea"
      string info();
  };
/*! \brief Class for porphyrin dendrimer
 *
 */
  class glu3 :public macromolecule {
    public:
      polymer chains;  //Glutamic chains
      macromolecule core;
      glu3(container &, inputfile &);
      string info();
  };

#ifdef HYPERSPHERE
  /*! \brief Hypersphere groups
   *  \author Martin Trulsson
   *  \date Lund, 2009
   */
  class hypermolecule : public macromolecule {
    public:
      void add(container &, vector<particle>, bool=false); //!< Add particle vector to group and container
      void move(container &, point);                       //!< Translate group
      void rotate(container &, double, double=0);          //!< Rotate group
  };
#endif
}//namespace
#endif

