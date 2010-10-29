#ifndef FAU_GROUP_H
#define FAU_GROUP_H

#include "faunus/common.h"
#include "faunus/average.h"
#include "faunus/particles.h"
#include "faunus/point.h"

namespace Faunus {
  class container;
  class slit;
  class inputfile;

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
      int beg,end;    //!< Define range in particle vector. [beg;end]
      string name;          //!< Informative name. Avoid spaces.
      group(int=0);         //!< Constructor, initialize data.

      void set(int, int);                        //!< Set particle range, "beg" and "end".
      bool find(unsigned int) const;             //!< Check if particle is part of the group
      int size() const;                          //!< Number of particles in group
      double charge(const vector<particle> &);   //!< Calculate total charge
      virtual int random();                      //!< Picks a random particle within this group
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
      void add(container &, vector<particle>, bool=false);//!< Add a particle vector at random position
      void add(container &, unsigned char, int);          //!< Add particles w. collision check
      bool swap(container &, group &);                    //!< Swap location of two groups
      int count(vector<particle> &, unsigned char);       //!< Count number of specific particles
      virtual int displace(container&, point);            //!< Displace random particle
      virtual void isobaricmove(container &, double){};   //!< Pressure scaling
      virtual int nummolecules();                         //!< Number of molecules
      unsigned short numhydrophobic(vector<particle> &);  //!< Number of hydrophobic particles
      bool swap(container &, int);                        //!< Move group to a new position
      bool saveCharges(string filename, vector<particle> &p); //!< Save all charges in group to disk 
      bool loadCharges(string filename, vector<particle> &p); //!< Load all charges in group from disk 
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
      average<double> rg2;  //!< Radius of gyration squared
      average<double> ree2; //!< End-to-end distance squared

      string info();                              //!< Show info
      string info(container &);                   //!< Show more info!
      void center(container &);                   //!< Center group in origo (0,0,0)
      double charge(const vector<particle> &);    //!< Calculate total charge
      double getcharge(const vector<particle> &); //!< Calculate total charge
      double radius(vector<particle> &);          //!< Calculate radius
      double gradius(vector<particle> &);         //!< Calculate radius of gyration
      double sqmassgradius(container &);          //!< Calculate mass weighted squared radius of gyration
      double sqend2enddistance(container &);      //!< Calculate squared end-to-end distance
      double vradius(vector<particle> &);         //!< Volume based protein radius
      double dipole(vector<particle> &);          //!< Calculate dipole moment
      double dipole(const container &);           //!< Calculate dipole moment
      void zmove(container &, double);            //!< Move in z-direction, only
      virtual void rotate(container &, double, double=0); //!< Rotate around a point
      void rotate(container &, point, double, double);
      void rotate(container &, point, point, double); //!< Rotate around arbitrary point
      void transrot(container &, double, double); 
      using group::add;
      void add(container &, inputfile &);         //!< Add according to inputfile
      macromolecule& operator=(const group&);     //!< Copy from group
      virtual void isobaricmove(container &,double);  //!< Displace CM with scale difference
      virtual int nummolecules();
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
      group sel;                 //!< A temporary group class
    public:
      molecules(unsigned short); //!< Constructor. Specify number of atoms in each molecule.
      unsigned short numatom;    //!< Number of atoms in each molecule
      int random();              //!< Pick a random molecule (NOT atom)
      string info();             //!< Show information
      group operator[](int);     //!< Access n'th molecule
      void add(container &, vector<particle> &, short=1);
      vector<int> pick(int );
  };

  /*! \brief Class for polymer molecules
   *  \author Mikael Lund
   *  \date Lund 2009
   */
  class polymer : public macromolecule {
    public:
      vector< vector<int> > nb; //!< Neighbor list
      polymer();
      vector<int> neighbors(int) const;
      string info();
      bool areneighbors(int, int) const;
      string getVMDBondScript();                 //!< Print TCL script for VMD to create bonds
#ifdef BABEL
      bool babeladd( container &, inputfile & ); //!< Load molecule from disk using OpenBabel
#endif
  };

  /*!
   * \brief Class for phospholipid-membrane
   * \author Bjoern Persson
   * \date Lund 2009
   */
  class popscmembrane : public group {
    protected:
      double scratio;               //!< Ratio of pops in percent
      double headarea;              //!< Headgroup area (i.e. density)

    public:
      vector<polymer> pops;         //!< Pops polymers
      vector<polymer> popc;         //!< Popc polymers
      popscmembrane();
      void load(inputfile&, slit&); //!< Scanns input object for "scratio", "headarea"
      string info();                //!< Get info string
      string info(slit &);          //!< Get expanded info string
      string getVMDBondScript();    //!< Print TCL script that tells VMD to draw bonds
  };
  
  /*!
   * \brief Class for porphyrin dendrimer
   */
#ifdef BABEL
  class glu3 :public macromolecule {
    public:
      polymer chains;  //Glutamic chains
      macromolecule core;
      glu3(container &, inputfile &);
      string info();
  };
#endif

#ifdef HYPERSPHERE
  /*!
   *  \brief Hypersphere groups
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

