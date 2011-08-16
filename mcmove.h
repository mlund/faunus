#ifndef FAUNUS_MCMOVE_H
#define FAUNUS_MCMOVE_H

#include <faunus/common.h>

namespace Faunus {
  class average;
  class unitcheck;

  /*!
   * \brief Optimize Monte Carlo displacement parameters for optimal mean square displament
   * \date Lund, 2011
   * \author Mikael Lund
   */
  class displacement_optimizer {
    private:
      double* dpPtr;   //!< Pointer to displacement parameter
      average<double>* msqPtr; //!< Pointer to mean-square displacement
      double dp_min;
      double dp_max;

    public:
      set_target(double&, average&);
  };

  class mcmove {
    protected:
      string title;                    //!< title of move
      string cite;                     //!< litterature reference, url, DOI etc.
      unsigned long int cnt;           //!< total number of trial moves
      unsigned long int cnt_accepted;  //!< number of accepted moves
      char iw;                         //!< width of first column of info string
      double dusum;                    //!< Sum of all energy changes made by this move

      virtual void trialmove()=0;      //!< Do a trial move
      virtual double energychange()=0; //!< Returns energy change of trialmove
      virtual void acceptmove()=0;     //!< Accept move, store new coordinates etc.
      virtual void rejectmove()=0;     //!< Reject move, revert to old coordinates etc.

    public:
      mcmove(string);             //!< Constructor
      virtual double move()=0;    //!< Attempt a move and return energy change
      double runfraction;         //!< Fraction of times calling move() should result in an actual move
      string info();              //!< Returns information string
      void unittest(unittest &);  //!< Perform unit test
  };

  mcmove::mcmove(string pfx) {
    cnt=cnt_accepted=0;
    dusum=0;
    iw=15;
  }

  void mcmove::unittest(unittest&) {

  }

  string mcmove::info() {
    std::ostringstream p, o;
    p << "#   " << setw(iw);
    o << "# " << title << endl
      << p << "More information:" << "  " << cite << endl
      << p << "Runfraction" << "= " << runfraction << endl;
    if (cnt>0) {
      o << p << "Number of trials" << "= " << cnt << endl
        << p << "Acceptance" << "= " << double(cnt_accepted)/cnt << endl
        << p << "Total energy change" << "= " << dusum << " kT" << endl;
    }
    return o.str();
  }

}//namespace
#endif
