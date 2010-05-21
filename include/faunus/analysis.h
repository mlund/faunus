#ifndef FAU_ANALYSIS_H
#define FAU_ANALYSIS_H

//#include "faunus/histogram.h"
#include "faunus/average.h"
#include "faunus/container.h"
#include "faunus/slump.h"
#include "faunus/io.h"
#include "faunus/xytable.h"
#include "faunus/hardsphere.h"
#include "faunus/energy/base.h"

namespace Faunus {
  /*!
   * Base class for analysis functions
   * \author Mikael Lund
   * \date Prague, 2007
   */
  class analysis {
    protected:
      slump slp;
      bool runtest();             //!< True if we should run, false if not.
    public:
      float runfraction;          //!< Fraction of times analysis should be run
      virtual string info()=0;    //!< Information/results
      analysis() { runfraction=1; }
  };

  /*!
   * This class can encompass several distributions
   * of average values (y) as a function of some variable (x)
   * Each distribution is automatically identified by specifying
   * an arbitrary name while adding to it.
   * This will typically be used to average a certain property as
   * a function of, say, an intermolecular separation. It is assumed
   * that the x range is identical for all distributions.
   *
   * \author Mikael Lund
   * \date December 2007
   */
  class distributions : public analysis {
    private:
      float xmax, xmin; // maximum/minimum encountered x value
      float xmax_set, xmin_set, dx;
      io fio;
      vector<string> s;
      vector< xytable<float,average<float> > > d;
      unsigned short find(string);
    public:
      distributions(float=0.5, float=0, float=0); //!< Constructor
      bool add(string, float, float);  //!< Add value to distribution
      bool write(string);  //!< Write distributions to a file
      string info();       //!< Write distributions to a string
  };

  class gfactor : public analysis {
    private:
      average<float> g;
      double N,V,mu2,MU2;
      macromolecule mol;
      double gamma;
      double gscale;
  //    histogram h;
    public:
      gfactor(container &, molecules &) ;//: h(float(0.1), float(0.), float(40.));
      double add(vector<particle> &, molecules &);
      string info();
  };

  class systemenergy : public analysis {
    private:
      double u0;
      vector<double> confu;        //!< Vector to track system energy in time
      io fio;
    public:
      double cur;
      double sum;                 //!< Initial energy + all changes
      average<double> uavg;
      systemenergy(double);
      void update(double);        //!< Specify current system energy and recalc averages
      void track();               //!< Add a time element to confu
      void operator+=(double);    //!< Add system energy change
      string info();              //!< Info
      void write(); 
      string confuout();
      double drift();             //!< Measured energy drift
      void check(checkValue &);   //!< Output testing
  };

  /*!
   * \brief Calculates free energy pathway
   * \author Mikael Lund
   * \note Can't really remember what this is good for...ML
   *
   * Starting from some point this class will generate eight
   * other points around it (cubic) and calculate the excess chemical
   * potential in these eight points. The one with the lowest value
   * will be the center of a point that will be used as the center
   * for eight new points etc. This will go on until the path reaches a
   * target point as specified in the constructor.
   */
  template<class T_pairpot> class widompath : public analysis {
    private:
      std::ostringstream ostr;
      bool finished;
      unsigned short cnt;
      float r;
      particle goal;
      vector<float> trajmu;
      vector<point> cube, trajpos;
      vector< average<float> > cubemu;
      void setcube(point &);
    public:
      particle p;                         //!< Particle to insert
      widompath(particle &, particle &);  //!< Constructor
      void povpath(iopov&);               //!< Plot trajectory in povray
      void update(container&, interaction<T_pairpot>&, iopov&);
      string info();                      //!< Information
  };

  //! \param beg Particle to insert, including starting position
  //! \param end Destination particle - Pathway ends here!
  template<class T_pairpot> widompath<T_pairpot>::widompath(particle &beg, particle &end) {
    r=5.;
    runfraction=0.2;
    finished=false;
    cnt=100;
    p=beg;
    goal=end;
    setcube(p);
  }
  template<class T_pairpot> void widompath<T_pairpot>::update(
      container &con, interaction<T_pairpot> &pot, iopov &pov) {
    if (runtest()==false || finished==true)
      return;
    unsigned char i,imax;
    if (cnt>0) {
      cnt--;
      for (i=0; i<cube.size(); i++) {
        p=cube[i]; 
        cubemu[i] += (con.collision(p)==false) ?
          exp(-pot.energy(con.p, p)) : exp(-100);
      }
    } else {
      imax=0;
      for (i=1; i<cube.size(); i++)
        if (cubemu[i].avg()>cubemu[imax].avg())
          imax=i;
      p=cube[imax];
      trajpos.push_back(p); 
      trajmu.push_back( -log(cubemu[imax].avg()) );
      //cout << p << " " << -log(cubemu[imax].avg()) << endl;
      setcube(p);
      cnt=100;
      if (p.dist(goal)<4*p.radius) {
        finished=true;
        //povpath(pov);
      }
    }
  }
  template<class T_pairpot> void widompath<T_pairpot>::povpath(iopov &pov) {
    for (unsigned short i=0; i<trajpos.size()-1; i++)
      pov.connect(trajpos[i], trajpos[i+1], 0.5);
    //pov.connect(trajpos[trajpos.size()-1], goal, 0.5);
  }
  template<class T_pairpot> void widompath<T_pairpot>::setcube(point &p) {
    cube.resize(8);
    cubemu.resize( cube.size() );
    for (unsigned char i=0; i<cube.size(); i++) {
      cubemu[i].reset();
      cube[i].x = p.x + r*slp.random_half();
      cube[i].y = p.y + r*slp.random_half();
      cube[i].z = p.z + r*slp.random_half();
    } 
  }

  template<class T_pairpot> string widompath<T_pairpot>::info() {
    std::ostringstream o;
    o << "# WIDOMPATH ANALYSIS:" << endl
      << "#   Ghost particle (z,r) = " << p << " " << p.charge << " " << p.radius << endl
      << "#   Target position      = " << goal << endl
      << "#   Points in trajectory = " << trajpos.size() << endl
      << "#   Target reached       = ";
    if (finished==true) o << "Yes";
    else o << "NO!";
    o << endl;
    for (int i=0; i<trajpos.size(); i++)
      o << i << " " << trajmu[i] << " " << trajpos[i] << endl;
    return o.str();
  }

  /*!
   * This class analyses angular correlations between
   * any two macromolecules.
   * \author Mikael Lund
   * \date ADFA, 2008
   */
  class angularcorr : public analysis {
    private:
      point m1,m2;
    public:
      angularcorr();
      void update(macromolecule &, macromolecule &, distributions &);
  };

  /*!
   * This will analyse the avg. concentration of a particle
   * type or group in a spherical shell around a central particle
   * (binding site)
   *
   * \author Mikael Lund
   * \date ADFA, 2008
   */
  class twostatebinding : public analysis {
    private:
      float r2, vol;
    public:
      average<float> conc; //!< Local or "bound" concentration
      twostatebinding(double);
      void update(container &, point &, group &);
      void update(container &, point &, vector<macromolecule> &, int=0);
      void update(container &c, point &, unsigned char);
      string info();
      string info(double);
  };

  class aggregation :public analysis {
    private:
      vector<macromolecule*> g;      // Vector of pointers to all macro.mol. in analysis
      vector<average<double> > RG2;  // Average radius of gyration square of N-particle cluster
      container *con;                // Pointer to the container
      vector<int> dist;              // Vector to store the histogram
      vector<macromolecule*> agg;    // Vector of pointers to molecules in aggregate
      vector<macromolecule*> unagg;  // --  //  -- to those not yet identified in a aggregate
      double CNT;                    // Counter
      hardsphere coll;               // Hardsphere booleans to define aggregates
      double sep;                    // Cluster definition, every mol. sep. by less than this
    public:                          // (in between any 'atom') is considered part of an aggregate
      aggregation(container &, vector<macromolecule> &i, double);
      void count();                  // Update the averages and histogram
      void write(string);            // Print result to 'string'
      string info() {                // Info
        std::ostringstream o;
        o << endl
          << "# AGGREGATION COUNT"<< endl
          << "#     Cluster definition = "<<sep<<" AA"<<endl
          << "#     Analysis performed "<<CNT<<" times."<<endl;
        return o.str();
      }
  };

  /* Calculates the virial pressure by calculating the forces between
   * the particles. The force is calculated by taking the derivative
   * of the pair-potential and does therefore *not* work with
   * hardsphere potentials.
   *
   * \author Mikael Lund
   * \date Lund, 2008
   */
  class virial : public analysis {
    private:
      double conc;
    public:
      double dr; //!< r-step when taking the derivative of the pair potential.
      virial(container &);
      virial(container &, vector<macromolecule> &);
      average<double> pex; //!< Excess pressure
      void sample(container &, energybase &);
      void sample(container &, energybase &, vector<macromolecule> &);
      void check(checkValue &t); //!< Output checking
      string info();
  };

  class pointpotential : public analysis {
    private:
      struct data {
        point p;
        string name;
        average<double> phi, expphi;
      };
    public:
      vector<data> list;
      void add(point, string);
      void sample(container &, energybase &);
      string info();
  };

  /*
   * \author Bjorn Persson?
   * \todo Document this class
   */
  class diskoverlap : public analysis {
    private:
      int s1, s2, s3, i, j, k;
      vector< average < double > > size;
      vector< average < double > > asym;
      vector<double> sscale;
      vector<double> ascale;
      vector<double> scnt;
      vector<double> acnt;
      double cnt;
      point origin, dummy;
    public:
      diskoverlap(vector<point> &);
      void check(vector<point> &);
      void blockavg();
      string info();
  };
}//namespace
#endif
