#ifndef FAU_POT_TABLE_H
#define FAU_POT_TABLE_H

#include "faunus/potentials/base.h"
#include "faunus/xytable.h"
#include "faunus/potentials/pot_hccoulomb.h"
#include "faunus/potentials/pot_hsminimage.h"

namespace Faunus {
  /*!
   * \brief Load pair potentials from disk
   * \author Mikael Lund
   * \date Prague, 2008
   *
   * Class to load potential of mean force (PMF) data from file(s). If
   * the distance between the particles is outside the data range a 
   * simple Coulomb/Hardsphere potential will be applied.
   *
   * The PMF-file must contain a list of space-separated
   * distances and energies. No text is tolerated! The
   * data need not be equidistant and will be averaged
   * to fit the current resolution of 0.1 Aangstrom.
   *
   * If the data resolution is less than 0.1 AA or if
   * data is not found (for example from 0 to contact)
   * a default value of +30 kT will be assigned.
   *
   * \note Check your loaded potentials!
   *
   */
  class pot_table : public pot_hscoulomb {
    protected:
      xytable<double, double> pmf[particle::LAST][particle::LAST];
      double xres; //!< Distance resolution
      double nan;  //!< Energy of no data found.
    public:
      string pmfdir; //!< Directory containing PMF data

      pot_table(const inputfile &in) : pot_hscoulomb(in) {
        xres=0.1;
        nan=30;
        name+="/Empirical data potential";
        pmfdir=in.getstr("pmfdir","./");
        for (int i=particle::FIRST; i<particle::LAST; i++)
          for (int j=particle::FIRST; j<particle::LAST; j++)
            pmf[i][j].xres=xres;
      }

      double pairpot (const particle &p1, const particle &p2) {
        unsigned short i=p1.id,j=p2.id;
        if (i>j)
          std::swap(i,j);
        if (pmf[i][j].xmax()<0.01)
          return pot_hscoulomb::pairpot(p1,p2);
        double r2=p1.sqdist(p2), xmax=pmf[i][j].xmax()-0.5;
        return ( r2 > xmax*xmax ) ?
          pot_hscoulomb::pairpot(p1,p2) :  // use Coulomb pot. outside data
          pmf[i][j]( sqrt(r2) );           // ...else use table.
      }

      bool loadpmf(const species &spc, particle::type t1, particle::type t2) {
        double x,y;
        xytable< double, Faunus::average<double> > avgpmf(xres);
        string file = pmfdir + "/" + spc.d[t1].name + "-" + spc.d[t2].name + ".dat";
        if (t1>t2)
          std::swap(t1,t2);
        std::ifstream fs( file.c_str() );
        if (fs) {
          while (!fs.eof()) {
            fs >> x >> y;
            avgpmf(x)+=y;
          }
          fs.close();
          for (x=0; x<avgpmf.xmax(); x+=avgpmf.xres) {
            if (avgpmf(x).cnt==0)
              avgpmf(x)+=nan;
            pmf[t1][t2](x)=avgpmf(x).avg()/f;
          }
          if (pmf[t1][t2].xmax()>0)
            return true;
        }
        return false;
      }

      void loadpmf(const container &c) {
        vector<particle::type> id; // vector of particle id's found in container.
        for (unsigned short i=0; i<c.p.size(); i++) {
          vector<particle::type>::iterator iter = std::find(id.begin(), id.end(), c.p[i].id);
          if (iter==id.end())
            id.push_back(c.p[i].id);
        }
        for (unsigned short i=0; i<id.size(); i++)
          for (unsigned short j=i; j<id.size(); j++)
            if( loadpmf(c,id[i],id[j])==false )
              loadpmf(c,id[j],id[i]);
      }

      virtual string info() {
        std::ostringstream o;
        o << pot_hscoulomb::info()
          << "#   PMF directory     = " << pmfdir << endl;
        return o.str();
      }

      /*! Show info + list of loaded pmf data */
      virtual string info(species &spc) {
        std::ostringstream o;
        o << info()
          << "#   PMF Info: (a,b,resolution,r_max,file)" << std::endl; 
        for (unsigned short i=particle::FIRST; i<particle::LAST; i++)
          for (unsigned short j=particle::FIRST; j<particle::LAST; j++)
            if (pmf[i][j].xmax()>0) 
              o << "#     " << spc.d[i].name << " " << spc.d[j].name << " "
                << pmf[i][j].xres << " "
                << pmf[i][j].xmax() << std::endl;
        o << std::endl; 
        return o.str();
      }
  };

 /*!
   * \brief As Faunus::pot_table but with periodic boundaries.
   * \author Mikael Lund
   * \date Prague, 2008
   * \todo info() output could be more elegant.
   */
  class pot_tableMI : public pot_table {
    private:
      pot_hsminimage coulomb;
    public:
      pot_tableMI(const inputfile &in) : pot_table(in), coulomb(in) {};
      inline double pairpot(const particle &p1, const particle &p2) {
        unsigned short i=p1.id,j=p2.id;
        if (i>j)
          std::swap(i,j);
        if (pmf[i][j].xmax()<0.01)
          return coulomb.pairpot(p1,p2);
        double r2=p1.sqdist(p2), xmax=pmf[i][j].xmax()-0.5;
        return ( r2 > xmax*xmax ) ?
          coulomb.pairpot(p1,p2) :   // use Coulomb pot. outside data
          pmf[i][j]( sqrt(r2) );     // ...else use table.
      }
      string info() { return pot_table::info() + coulomb.info(); }
      string info(species &spc) { return pot_table::info(spc); }
  };
}//namespace
#endif
