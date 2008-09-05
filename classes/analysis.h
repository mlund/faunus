#ifndef FAU_ANALYSIS_H
#define FAU_ANALYSIS_H

#include "average.h"
#include "container.h"
#include "potentials.h"
#include "slump.h"
#include "io.h"
#include "xytable.h"
#include "hardsphere.h"

namespace Faunus {

/*!
 * Base class for analysis functions
 * \author Mikael Lund
 * \date Prague, 2007
 */
class analysis {
  protected:
    slump slp;
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

class systemenergy : public analysis {
  private:
    double u0,sum,cur;
    average<double> uavg, u2avg;
    vector<double> confu;  //!< Vector to track system energy in time
    io fio;
  public:
    systemenergy(double);
    void update(double);        //!< Specify current system energy and recalc averages
    void track();               //!< Add a time element to confu
    void operator+=(double);    //!< Add system energy change
    string info() {
      write();                  //!< Print dynamics of system energy 
      ostringstream o;
      o << endl << "# SYSTEM ENERGY (kT):" << endl;
      if (uavg.sum>0)
        o << "#   Averages <U> <U^2> = " << uavg.avg() << " " << u2avg.avg() << endl
          << "#   sqrt(<U^2>-<U>^2)  = " << sqrt(u2avg.avg()-uavg.avg()*uavg.avg()) << endl;
      o << "#   Initial energy     = " << u0 << endl
        << "#   Initial + changes  = " << sum << endl
        << "#   Current energy     = " << cur << endl
        << "#   Absolute drift     = " << abs(cur-sum) << endl;
      return o.str();
    }
    string confuout() {
      int j=confu.size();
      ostringstream o;
      o << endl << "# SYSTEM ENERGY (kT):"<< endl;
      for (int i=0;i<j;i++)
        o << i+1 << " " << confu[i] << endl;
      return o.str();
    }
    void write(); 
};

/*! \brief Widom method for excess chemical potentials
 *  \author Mikael Lund
 *  \todo Expand with a corrected one-particle insertion (Woodward+Svensson's charge scaling)
 *
 *  This class will insert "ghost" particle(s) so as to
 *  calculate the (mean) excess chemical potential.
 */
template<class T_pairpot>
class widom : public analysis {
  private:
    unsigned long long int cnt;
    double u;
    particle a,b;
    container *con;
    average<double> expsum; 
    interaction<T_pairpot> *pot;
  public:
    widom(container &c, interaction<T_pairpot> &i,
        particle::type t1, particle::type t2) {
      cnt=0;
      con=&c;
      pot=&i;
      a=con->get(t1);
      b=con->get(t2);
    }
    string info();                              //!< Print results of analysis
    double muex() { return -log(expsum.avg()); }//!< Excess chemical potential
    double gamma() { return exp(muex()); }      //!< Activity coefficient
    void insert(unsigned short=100);            //!< Widom insertions
};
//! Insert a salt pair and evaluate the excess chemical potential
//! \param n Number of insertions
template<class T>
void widom<T>::insert(unsigned short n) {
  while (n>0) {
    n--;
    cnt++;
    con->randompos(a);
    con->randompos(b);
    u=pot->energy(con->p, a) +
      pot->energy(con->p, b) +
      pot->pair.pairpot(a,b)*pot->pair.f;
    expsum+=exp(-u);
  }
}
template<class T>
string widom<T>::info() {
  ostringstream o;
  o << endl
    << "# WIDOM PARTICLE INSERTION ANALYSIS:" << endl
    << "#   Number of insertions = " << cnt << endl
    << "#   Ion pair charges     = " << a.charge << ", " << b.charge << endl
    << "#   Excess chemical pot. = " << muex()  << endl
    << "#   Mean activity coeff. = " << gamma() << endl;
  return o.str();
}

/*!
 * \brief Calculates free energy pathway
 * \author Mikael Lund
 *
 * Starting from some point this class will generate eight
 * other points around it (cubic) and calculate the excess chemical
 * potential in these eight points. The one with the lowest value
 * will be the center of a point that will be used as the center
 * for eight new points etc. This will go on until the path reaches a
 * target point as specified in the constructor.
 */
template<class T_pairpot>
class widompath : public analysis {
  private:
    ostringstream ostr;
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
template<class T_pairpot>
widompath<T_pairpot>::widompath(particle &beg, particle &end) {
  r=5.;
  runfraction=0.2;
  finished=false;
  cnt=100;
  p=beg;
  goal=end;
  setcube(p);
}
template<class T_pairpot>
void widompath<T_pairpot>::update(
    container &con, interaction<T_pairpot> &pot, iopov &pov) {
  if (slp.runtest(runfraction)==false || finished==true)
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
template<class T_pairpot>
void widompath<T_pairpot>::povpath(iopov &pov) {
  for (unsigned short i=0; i<trajpos.size()-1; i++)
    pov.connect(trajpos[i], trajpos[i+1], 0.5);
  //pov.connect(trajpos[trajpos.size()-1], goal, 0.5);
}
template<class T_pairpot>
void widompath<T_pairpot>::setcube(point &p) {
  cube.resize(8);
  cubemu.resize( cube.size() );
  for (unsigned char i=0; i<cube.size(); i++) {
    cubemu[i].reset();
    cube[i].x = p.x + r*slp.random_half();
    cube[i].y = p.y + r*slp.random_half();
    cube[i].z = p.z + r*slp.random_half();
  } 
  /*}
  cube[0].x+=r;
  cube[0].y+=r;
  cube[0].z+=r;
  cube[1]=-cube[0]+p;
  cube[2].x-=r;
  cube[2].y+=r;
  cube[2].z+=r;
  cube[3]=-cube[2]+p;
  cube[4].x+=r;
  cube[4].y-=r;
  cube[4].z+=r;
  cube[5]=-cube[4]+p;
  cube[6].x+=r;
  cube[6].y+=r;
  cube[6].z-=r;
  cube[7]=-cube[6]+p;
  cout << p << endl << cube[0] << endl << cube[1] << endl << cube[2]
    << endl << cube[3] << endl << cube[4] << endl << cube[5] << endl
    << cube[6] << endl << cube[7] << endl;
    */
}
template<class T_pairpot>
string widompath<T_pairpot>::info() {
  ostringstream o;
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
    twostatebinding(double radius) {
      r2=radius*radius;
      vol=4./3.*acos(-1)*pow(radius,3);
    }
    void update(container &con, point &site, group &g) {
      unsigned short i,cnt=0;
      for (i=g.beg; i<=g.end; i++)
        if (con.sqdist(site,con.p[i])<r2) {
          cnt=1;
          break;
        }
      conc+=float(cnt)/vol;
    }
    void update(container &con, point &site, vector<macromolecule> &g, int first=0) {
      unsigned short i,j,cnt=0;
      for (j=first; j<g.size(); j++)
        for (i=g[j].beg; i<=g[j].end; i++)
          if (con.sqdist(site, con.p[i])<r2) {
            cnt+=1;
            break;
          }
      conc+=float(cnt)/vol;
    }
    void update(container &con, point &site, particle::type id) {
      unsigned short i,n=con.p.size(),cnt=0;
      for (i=0; i<n; i++)
        if (con.p[i].id==id)
          if (con.sqdist(site,con.p[i])<r2) cnt++;
      conc+=float(cnt)/vol;
    }
    string info() { return string(0); }
    string info(double bulkconc) {
      ostringstream o;
      o << endl << "# TWOSTATE BINDING ANALYSIS:" << endl
        << "#   More information:  J. Phys. Chem. 1995, 99, 10412" << endl
        << "#   Site radius      = " << sqrt(r2) << endl
        << "#   Avg. site conc.  = " << conc.avg() << endl;
      if (bulkconc>0)
        o << "#   Site excess      = " << conc.avg()/bulkconc
          << " (" << -log(conc.avg()/bulkconc) << " kT)" << endl;
      return o.str();
    }
};

class aggregation :public analysis {
  private:
    vector<macromolecule*> g;
    vector<average<double> > RG2;
    container *con;
    vector<int> dist;  // Vector to store the histogram
    vector<macromolecule*> agg;
    vector<macromolecule*> unagg;
    double CNT;
    hardsphere coll;
    double sep;
  public:
    aggregation(container &, vector<macromolecule> &i, double);
    void count();
    void write(string);
    string info() {
      ostringstream o;
      o << endl
        << "# AGGREGATION COUNT"<< endl<<endl;
    }
};
}//namespace
#endif
