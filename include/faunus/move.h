#ifndef FAUNUS_MCMOVE_H
#define FAUNUS_MCMOVE_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/average.h>

namespace Faunus {
  //class average;
  class unittest;
  class Energybase;

  namespace Move {

    /*!
     * \brief Optimize Monte Carlo displacement parameters for optimal mean square displament
     * \date Lund, 2011
     * \author Mikael Lund
     */
    class displacement_optimizer {
      private:
        double* dpPtr;   //!< Pointer to displacement parameter
        //average<double>* msqPtr; //!< Pointer to mean-square displacement
        double dp_min;
        double dp_max;

      public:
        //void set_target(double&, average&);
    };

    /*!
     * \brief Base class for Monte Carlo moves
     * \author Mikael Lund
     * \date Lund, 2007-2011
     *
     * The is a base class and derived classes MUST implement a number of
     * functions. Importanty, it is important to have an argumentfree move()
     * function as this can be attached to groups so that these have the
     * information to move automatically.
     */

    class Movebase {
      protected:
        Energy::Energybase* pot;         //!< Pointer to energy functions
        Space* spc;
        string title;                    //!< title of move
        string cite;                     //!< litterature reference, url, DOI etc.
        string prefix;                   //!< inputmap prefix
        unsigned long int cnt;           //!< total number of trial moves
        unsigned long int cnt_accepted;  //!< number of accepted moves
        char w;                          //!< info string text width
        double dusum;                    //!< Sum of all energy changes made by this move
        const double infty;              //!< Large value to represent infinity

        virtual void trialMove()=0;      //!< Do a trial move
        virtual void acceptMove()=0;     //!< Accept move, store new coordinates etc.
        virtual void rejectMove()=0;     //!< Reject move, revert to old coordinates etc.
        virtual double energyChange()=0; //!< Returns energy change of trialMove
        bool run() const;                //!< Runfraction test
        bool metropolis(const double &) const; //!< Metropolis criteria

      public:
        Movebase(Energy::Energybase&, Space&, string);             //!< Constructor
        virtual ~Movebase() {};
        virtual double move()=0;     //!< Attempt a move and return energy change
        double runfraction;          //!< Fraction of times calling move() should result in an actual move
        string info();               //!< Returns information string
        virtual double totalEnergy();//!< Total energy relevant for energy drift tracking
        //void unittest(unittest&);  //!< Perform unit test
    };

    /*
    class translate : public Movebase {
      protected:
        void trialMove();
        void acceptMove();
        void rejectmove();
        double energychange();
      public:
        translate(string="translate_", Energybase&, space&);
        unsigned int group;
        point dp;                   //!< Displacement vector
        double move();
        string info();
    };
    */

    /*!
     * \brief Translation of single particles or single particles in a group
     * \author Mikael Lund
     * \date Lund, 2011
     *
     * To move a single particle, specify its position in the space particle vector
     * in iparticle. For randomly moving all particles in a group (typically salt),
     * point igroup to the appropriate group in the space class g vector.
     */
    class ParticleTranslation : public Movebase {
      private:
        typedef std::map<short, average<double> > map_type;
        map_type accmap; //!< Single particle acceptance map
        map_type sqrmap; //!< Single particle mean square displacement map
        Group* igroup;   //!< Group pointer in which particles are moved randomly (NULL if none, default)
        int iparticle;   //!< Select single particle to move (-1 if none, default)
      protected:
        void trialMove();
        void acceptMove();
        void rejectMove();
        double energyChange();
      public:
        ParticleTranslation(InputMap&, Energy::Energybase&, Space&, string="mv_particle");
        void setGroup(Group&); //!< Select group in which to randomly pick particles from
        void setParticle(int); //!< Select single particle in p_vec to move
        double move();         //!< Move selected particle once or n times in selected group of length n
        point dir;             //!< Displacement directions (default: x=y=z=1)
        string info();
        double totalEnergy();  //!< Total energy for drift tracking
    };

    class MoleculeTranslation : public Movebase {
    };

    class MoleculeRotation : public Movebase {
    };


  }//namespace
}//namespace
#endif
