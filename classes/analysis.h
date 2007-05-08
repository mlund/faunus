#ifndef _ANALYSIS_H
#define _ANALYSIS_H

#include "average.h"
#include "container.h"
#include "potentials.h"
#include "io.h"

class analysis {
  public:
    virtual string info()=0;
};

class systemenergy : public analysis {
  private:
    double u0,sum,cur;
    average<double> uavg, u2avg;
  public:
    systemenergy(double);
    void setcurrent(double);    //!< Specify current system energy
    void operator+=(double);    //!< Add system energy change
    string info() {
      uavg+=sum;
      u2avg+=sum*sum;
      ostringstream o;
      o << "# System energy (kT):" << endl
        << "#   Averages <U> <U^2> = " << uavg.avg() << " " << u2avg.avg() << endl
        << "#   Initial energy     = " << u0 << endl
        << "#   Initial + changes  = " << sum << endl
        << "#   Current energy     = " << cur << endl
        << "#   Absolute drift     = " << abs(cur-sum) << endl;
      return o.str();
    }
};

/*! \brief Widom method for excess chemical potentials
 *  \author Mikael Lund
 *  \todo Expand to one-particle insertion (Woodward+Svensson)
 *
 *  This class will insert "ghost" particles so as to
 *  calculate the mean excess chemical potential.
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
  o << "# Widom Analysis:" << endl
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
    point goal;
    vector<float> trajmu;
    vector<point> cube, trajpos;
    vector< average<float> > cubemu;
    void update(container&, interaction<T_pairpot>&);
    void setcube(point &);
  public:
    particle p;                       //!< Particle to insert
    widompath();
    void povpath(iopov&);             //!< Creates a povray trajectory
};

template<class T_pairpot>
widompath<T_pairpot>::widompath() {
  finished=false;
  cnt=100;
  p.charge=+1;
  p.radius=2.0;
  cube.resize(8);
  cubemu.resize( cube.size() );
}

template<class T_pairpot>
void widompath<T_pairpot>::update(
    container &con, interaction<T_pairpot> &pot) {
  unsigned char i,imax;
  if (finished==false)
    return;
  if (cnt>0) {
    cnt--;
    for (i=0; i<cube.size(); i++) {
      p=cube[i]; 
      if (con.collision(p)==false)
        cubemu[i]+=exp(-pot.energy(con.p, p));
      else cubemu[i]+=1e3;
    }
  } else {
    imax=0;
    for (i=1; i<cube.size(); i++)
      if (cubemu[i].avg()>cubemu[imax].avg())
        imax=i;
    p=cube[imax];
    trajpos.push_back(p); 
    trajmu.push_back( -log(cubemu[imax].avg()) );
    setcube();
    cnt=100;
    if (p.dist(goal)<p.radius) {
    }
  }
}

template<class T_pairpot>
void widompath<T_pairpot>::povpath(iopov &pov) {
  for (unsigned short i=0; i<trajpos.size()-1; i++)
    pov.connect(trajpos[i], trajpos[i+1], 0.5);
  pov.connect(trajpos[trajpos.size()-1], goal, 0.5);
}

template<class T_pairpot>
void widompath<T_pairpot>::setcube(point &p) {
  for (unsigned char i=0; i<cube.size(); i++) {
    cube[i]=p;
    cubemu[i].reset();
  }
  cube[0].x+=r;
  cube[0].y+=r;
  cube[0].z+=r;
  cube[1]=-cube[0];
  cube[2].x-=r;
  cube[2].y+=r;
  cube[2].z+=r;
  cube[3]=-cube[2];
  cube[4].x+=r;
  cube[4].y-=r;
  cube[4].z+=r;
  cube[5]=-cube[4];
  cube[6].x+=r;
  cube[6].y+=r;
  cube[6].z-=r;
  cube[7]=-cube[6];
}
#endif
