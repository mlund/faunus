#include "faunus/moves/base.h"

namespace Faunus {
   
  markovmove::markovmove(ensemble &e, container &c, energybase &i) : ens(&e), con(&c), pot(&i) {
    du=utot=dp=deltadp=0;
    cnt=naccept=0;
    runfraction=1;
    dp_opt=false;
    dp_min=0;
    dp_max=10;
    dp_width=(dp_max-dp_min)/10;
    dp_cnt=0;
    dp_N=1000;
    dp_dist.init(dp_width, dp_min, dp_max);
  }
  
  /*!
   * \param in Inputfile object to read parameters from. Each keyword is prefixed by the string "prefix".
   *
   * The following keywords are searched for:
   * - prefix_dp (displacement parameter)
   * - prefix_dpopt (switch dp optimization on or off)
   * - prefix_dpmin (minimum dp to sample)
   * - prefix_dpmax (maximum dp to sample)
   * - prefix_dpwidth (distance between dp points)
   * - prefix_dpsamples (number of dp parameters to sample)
   */
  void markovmove::getInput(inputfile &in) {
    runfraction = in.getflt(prefix+"runfrac",runfraction);
    dp     = in.getflt(prefix+"dp", dp);
    dp_opt = in.getboo(prefix+"dpopt", dp_opt);
    if (dp_opt==true) {
      dp_min   = in.getflt(prefix+"dpmin", dp_min);
      dp_max   = in.getflt(prefix+"dpmax", dp_max);
      dp_width = in.getflt(prefix+"dpwidth", dp_width);
      dp_N     = in.getint(prefix+"dpsamples", dp_N);
      dp_dist.init(dp_width, dp_min, dp_max);
      if (dp<dp_min || dp>dp_max)
        std::cerr << "Warning: " << prefix << "dp is out of range!" << endl;
    }
  }

  double markovmove::newdp() {
    double x=dp;
    if (dp_opt==true && dp_N>0) {
      x=slp.random_one()*(dp_max-dp_min) + dp_min;
      if (x<1e-5)
        x+=dp_width;
    }
    return x;
  }
  
  double markovmove::optimaldp() {
    double x=dp;
    if (dp_opt==true) {
      x=dp_dist.x_atmax_y();
      if (x<1e-5)
        x+=dp_width;
    }
    return x;
  }
  
  /*!
   * This virtual function should be called from all move functions. The base class version does the following:
   *  - Sets the energy changr to zero
   *  - Increases the counter for the move (cnt)
   *  - Performs displacement parameter optimization (if dp_opt is true)
   * \todo Call from all move functions!
   */
  double markovmove::move() {
    cnt++;
    du=0;
    if (dp_opt==true) {
      if (dp_N<=0)                  // do we still need to optimize?
        dp=optimaldp();             // no! set dp to optimum value
      else {
        dp_cnt++;                   // yes! increase counter for current dp
        if (dp_cnt>500.) {          // should we still sample this particular dp?
          dp_cnt=0;                 //   no! reset counter
          dp_dist(dp)+=dpsqr.avg(); //   store L^2 in table
          dpsqr.reset();            //   and reset average
          dp=newdp();               //   set new displacement parameter
          dp_N--;                   //   decrease total counter
        }
      }
    }
    return du;
  }

  string markovmove::info() {
    std::ostringstream o;
    o << endl
      << "# " << name << ":" << endl;
    if (cnt>0) {
      o.precision(2);
      o << "#   Acceptance                = " << accepted()*100 << " \%" << endl;
      o.precision(4);
      o << "#   Number of trials          = " << cnt << endl
        << "#   Runfraction               = " << runfraction*100 << " \%" << endl
        << "#   Energy change (kT)        = " << utot << " " << utot/cnt << " "
                                              << utot/(accepted()*cnt) << endl;
      if (dp!=0)
        o << "#   Displacement parameter    = " << dp << endl;
      if (cite.empty()==false)
        o << "#   More information:           " << cite << endl;
      if (dp!=0)
        if (dpsqr.cnt>0)
          o << "#   Average displacement      = " << sqrt(dpsqr.avg()) << endl
            << "#   Mean square displacement  = " << dpsqr.avg() << endl;
    }
    if (dp_opt==true) {
      o << "#   Displacement Optimization:" << endl
        << "#      Completed                      = " << ((dp_N>0) ? "no" : "yes") << endl
        << "#      Min/max displacement parameter = " << dp_min << " " << dp_max << " " << dp_width << endl
        << "#      Optimal displacement parameter = " << optimaldp()
        << " (L^2=" << dp_dist(optimaldp()).avg() << ")" << endl;
      if (cnt>dp_N) {
        o << "#      DP vs. L^2 distribution:" << endl;
        for (double x=dp_min; x<=dp_max; x+=3*dp_width)
          o << "#      " << std::setw(4) << x << " " << dp_dist(x).avg() << endl;
      }
    }
    return o.str();
  }

  bool markovmove::run(float p) {
    return (p>slp.random_one()) ? true : false;
  }

  float markovmove::accepted() {
    return naccept/float(cnt);
  }

  /*!
   * This function will adjust the displacement parameter in a way
   * that the acceptance ration lies within a certain tange. Useful
   * for equilibration runs -- do not use it in production runs.
   * \param max Maximum percentage of accepted moves
   * \warning This violates the detailed balance criteria.
   * \param min Minimum percentage of accepted moves
   * \author Mikael Lund
   * \todo Specify a maxmimum dp
   */
  void markovmove::adjust_dp(float min, float max) {
    float a=accepted()*100.;
    if (a>max) dp+=deltadp;
    if (a<min) dp-=deltadp;
    if (dp<=0) dp=deltadp;
  }

  void markovmove::check(checkValue &test) {
    string s = name;
    std::remove(s.begin(), s.end(), ' ');
    test.check(s + "Accepted", accepted() );
  }

}//namespace
