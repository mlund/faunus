#ifndef FAU_ANALYSIS_H
#define FAU_ANALYSIS_H

#include "faunus/histogram.h"
#include "faunus/average.h"
#include "faunus/slump.h"
#include "faunus/io.h"
#include "faunus/xytable.h"
#include "faunus/hardsphere.h"
#include "faunus/energy/base.h"

namespace Faunus {
  class container;
  class macromolecule;
  class point;
  class particle;
  class group;
  class histogram;

  /*!
   * \brief Base class for analysis functions
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
   * \brief Matrix of averages along some coordinate
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
      bool cntwrite(string); //!< Write counter distributions to a file
      string info();       //!< Write distributions to a string
      string cntinfo();    //!< Write the counter of each averaged value in the distribution to a string
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


  /*!
   * \brief Tracks total system energy and drifts
   */
  class systemenergy : public analysis {
    private:
      double u0;
      vector<double> confu;        //!< Vector to track system energy in time
      io fio;
    public:
      double cur;
      double sum;                 //!< Initial energy + all changes
      average<double> uavg;
      systemenergy();
      systemenergy(double);
      void initialize(double);    //!< Initialize all data
      void update(double);        //!< Specify current system energy and recalc averages
      void track();               //!< Add a time element to confu
      systemenergy & operator+=(double);    //!< Add system energy change
      string info();              //!< Info
      void write(); 
      string confuout();
      double drift();             //!< Measured energy drift
      void check(checkValue &);   //!< Output testing
  };

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
      vector<macromolecule*> g;        // Vector of pointers to all macro.mol. in analysis
      vector<int> dist;                // Vector to store number of counts of each aggregate size
      vector<average<double> > RG2;    // Average radius of gyration square of N-particle cluster
      vector<average<double> > contact;// Average number of contacts per protein in each cluster size
      vector<histogram> intdist;       // Internal distribution function on N-particle cluster
      histogram      avgintdist;       // Internal distribution averaged over all aggregates
      container *con;                  // Pointer to the container
      vector<macromolecule*> agg;      // Vector of pointers to molecules in aggregate
      vector<macromolecule*> unagg;    // --  //  -- to those not included or yet not found
      vector<macromolecule*> remaining;// --  //  -- to those remaining to be found in each iteration
      double CNT;                      // Counter
      hardsphere coll;                 // Hardsphere booleans to define aggregates
      double sep;                      // Cluster definition, every mol. sep. by less than this
    public:                            // (in between any 'atom') is considered part of an aggregate
      aggregation(container &, vector<macromolecule> &i, double);
      aggregation(container &, vector<polymer> &i, double);
      void count();                    // Update the averages and histogram
      void write(string);              // Print result to 'string'
      string info() {                  // Info
        std::ostringstream o;
        o << endl
          << "# AGGREGATION COUNT"<< endl
          << "#     Cluster definition = "<<sep<<" AA"<<endl
          << "#     Analysis performed "<<CNT<<" times."<<endl;
        return o.str();
      }
  };

  /*!
   * Calculates the virial pressure by calculating the forces between
   * the particles. The force is calculated by taking the derivative
   * of the pair-potential and does therefore *not* work with
   * hardsphere potentials.
   *
   * \brief Virial pressure analysis
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
      vector<data> list;
    public:
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
	

  /*!
   * \brief Analyse pairing of mobile particles
   * \author Mikael Lund
   * \date Asljunga 2010
   *
   * This class will look for pairs of designated species and analyse if
   * they are within a certain threshold. The number of free ions and
   * pairs are averaged so as to evaluate the thermodynamic association constant.
   */
  class pairing : public analysis {
  private:
    average<double> cpair;  //!< concentration of pairs
    average<double> c1;     //!< concentration of free type 1
    average<double> c2;     //!< concentration of free type 2
    double r2;              //!< threshold squared
    char pid1, pid2;        //!< pairs to look for
  public:
    pairing(inputfile &);            //!< Constructor - input via inputfile class
    pairing(string, string, double); //!< Constructor - hard coded input
    void sample(container &);        //!< Look for pairs and free particles
    string info();                   //!< Information string
  };

  /*!
   * \brief Osmotic pressure in the cell model.
   * \author Mikael Lund
   * \date Lund, 2010
   *
   * This class will analyse the concentration of mobile ions at the
   * boundary of the spherical simulation container. This density is
   * directly related to the osmotic coefficient.
   */
  class osmoticpressure : public analysis {
    private:
      unsigned int cnt;
      double width;                       //!< Width of the cell boundary
      vector<unsigned int> hist;          //!< Mobile concentration profile
      cell* cPtr;
      double getConc(double);             //!< Get concentration (in M) at radial distance from cell center
    public:
      osmoticpressure(cell &);            //!< Constructor
      average<double> rhoid;              //!< Average salt concentration at boundary and in bulk
      void sample(group &);               //!< Sample concentration profile
      string info();                      //!< Information string
      void check(checkValue &);           //!< Unit testing
      bool write(string);                 //!< Write concentration profile to disk
  };
}//namespace
#endif
