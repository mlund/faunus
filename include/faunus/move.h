#ifndef FAUNUS_MCMOVE_H
#define FAUNUS_MCMOVE_H

#include <faunus/common.h>
#include <faunus/point.h>

namespace Faunus {
  //class average;
  class unittest;
  class energybase;

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

    class movebase {
      protected:
        void pad(std::ostringstream&);   //!< Text padding for info() function
        Energy::energybase* pot;                 //!< Pointer to energy functions
        space* spc;
        string title;                    //!< title of move
        string cite;                     //!< litterature reference, url, DOI etc.
        string prefix;                   //!< inputfile prefix
        unsigned long int cnt;           //!< total number of trial moves
        unsigned long int cnt_accepted;  //!< number of accepted moves
        char iw;                         //!< width of first column of info string
        double dusum;                    //!< Sum of all energy changes made by this move

        virtual void trialmove()=0;      //!< Do a trial move
        virtual void acceptmove()=0;     //!< Accept move, store new coordinates etc.
        virtual void rejectmove()=0;     //!< Reject move, revert to old coordinates etc.
        virtual double energychange()=0; //!< Returns energy change of trialmove
        bool run() const;                //!< Runfraction test
        bool metropolis(const double &) const;

      public:
        movebase(Energy::energybase&, space&, string);             //!< Constructor
        virtual double move()=0;    //!< Attempt a move and return energy change
        double runfraction;         //!< Fraction of times calling move() should result in an actual move
        string info();              //!< Returns information string
        //void unittest(unittest&);   //!< Perform unit test
    };

    /*
    class translate : public movebase {
      protected:
        void trialmove();
        void acceptmove();
        void rejectmove();
        double energychange();
      public:
        translate(string="translate_", energybase&, space&);
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
    class translate_particle : public movebase {
      protected:
        void trialmove();
        void acceptmove();
        void rejectmove();
        double energychange();
      public:
        int igroup;     //!< Group in which particles are moved randomly (-1 if no group, default)
        int iparticle;  //!< Select single particle to move (-1 if none, default)
        point dir;      //!< Displacement directions (default: x=y=z=1)
        translate_particle(inputfile &in, Energy::energybase&, space&, string="mv_particle");
        double move();  //!< Move iparticle once or n times in igroup of length n
        string info();
    };

  }//namespace
}//namespace
#endif
