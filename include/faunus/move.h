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

    class movebase {
      protected:
        void pad(std::ostringstream&);   //!< Text padding for info() function
        energybase* potPtr;              //!< Pointer to energy functions
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

      public:
        movebase(string);             //!< Constructor
        virtual double move()=0;    //!< Attempt a move and return energy change
        double runfraction;         //!< Fraction of times calling move() should result in an actual move
        string info();              //!< Returns information string
        //void unittest(unittest&);   //!< Perform unit test
    };

    class translate : public movebase {
      protected:
        void trialmove();
        void acceptmove();
        void rejectmove();
        double energychange();
      public:
        point dp;                   //!< Displacement vector
        translate(std::string="translate_");
        double move();
        string info();
    };

  }//namespace
}//namespace
#endif
