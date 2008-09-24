#ifndef FAU_POT_DATAPMF_H
#define FAU_POT_DATAPMF_H

#include "faunus/xydata.h"
#include "faunus/potentials/base.h"
#include "faunus/potentials/pot_coulomb.h"

namespace Faunus {
  /*!
   * \brief Load pair potentials from disk
   * \author Mikael Lund
   * \date Prague, 2007
   *
   * Class to load potential of mean force (PMF) data from file(s). If
   * the distance between the particles is outside the data range a 
   * simple Coulomb + LJ potential will be applied.
   */
  class pot_datapmf : public pot_coulomb {
    private:
      xydata<double> pmfd[particle::LAST][particle::LAST];
    public:
      string pmfdir; //!< Directory containing PMF data
      pot_datapmf(const inputfile &in) : pot_coulomb(in) {
        name+="/Empirical data potential";
        pmfdir=in.getstr("pmfdir","./");
      }
      double pairpot (const particle &p1, const particle &p2) {
        unsigned short i=p1.id,j=p2.id;
        if (i>j)
          std::swap(i,j);
        if (pmfd[i][j].xmax<0.01)
          return pot_coulomb::pairpot(p1,p2);
        double r2=p1.sqdist(p2);
        return ( r2 > pmfd[i][j].xmax * pmfd[i][j].xmax ) ? 
          pot_coulomb::pairpot(p1,p2) :  // use Coulomb pot. outside data 
          pmfd[i][j].x2y(sqrt(r2));      // ...else use table. 
      }

      /*!\param spc Species class.
       * \param pmfdir Directory in which to search for PMF's
       * \param type Particle name to search for. I.e. "NA" or "CL"
       */
      void loadpmf(const species &spc, particle::type id) {
        string n1,n2;
        unsigned short i,j=id;
        for (i=particle::FIRST; i<particle::LAST; i++) {
          n1=spc.d[i].name+"-"+spc.d[j].name;
          n2=spc.d[j].name+"-"+spc.d[i].name;
          if (loadpmf( spc, n1 + ".dat")==false)
            loadpmf( spc, n2 + ".dat");
        }
      }

      /* Search through particle vector and attempt to load pmf's for all possible
       * particle combinations.
       */
      void loadpmf(const container &c) {
        vector<particle::type> id;
        for (unsigned short i=0; i<c.p.size(); i++) {
          vector<particle::type>::iterator iter = std::find(id.begin(), id.end(), c.p[i].id);
          if (iter==id.end())
            id.push_back(c.p[i].id);
        }
        for (unsigned short i=0; i<id.size(); i++)
          for (unsigned short j=i; j<id.size(); j++)
            loadpmf(c, c.d[id[i]].name+"-"+c.d[id[j]].name+".dat"); 
      }

      /*! Load PMF(s) from a file. File format: Each set starts
       * with "#$ type1 type2 length". Several sets can be present
       * in the same file.
       */
      bool loadpmf(const species &spc, string filename) {
        filename=pmfdir+"/"+filename;
        string s,a_str,b_str;
        int a,b,len;
        vector<double> x,y;
        std::ifstream fh(filename.c_str());
        if (fh) {
          //scan file, word-by-word
          while (!fh.eof()) {
            fh >> s;
            if (s.find("#$")==0) {
              s.clear();
              fh >> a_str >> b_str >> len;
              x.resize(len);
              y.resize(len);
              a = spc.id(a_str);
              b = spc.id(b_str);
              if (a>b)
                std::swap(a,b);
              for (int i=0; i<len; i++) {
                fh >> x[i] >> y[i];
                y[i]=y[i]/f; // unit: kT/f
              };
              if (pmfd[a][b].xmax==0) {
                pmfd[a][b].add( x, y );
                pmfd[a][b].comment=filename;
              };
            };
          };
          fh.close();
          return true;
        }
        return false;
      }

      string info() {
        std::ostringstream o;
        o << pot_lj::info()
          << "#   Bjerrum length    = " << f << endl
          << "#   PMF directory     = " << pmfdir << endl;
        return o.str();
      }

      /*! Show info + list of loaded pmf data */
      string info(species &spc) {
        std::ostringstream o;
        o << info()
          << "#   PMF Info: (a,b,resolution,r_max,file)" << std::endl; 
        for (unsigned short i=particle::FIRST; i<particle::LAST; i++)
          for (unsigned short j=particle::FIRST; j<particle::LAST; j++)
            if (pmfd[i][j].xmax>0) 
              o << "#     " << spc.d[i].name << " " << spc.d[j].name << " "
                << pmfd[i][j].res << " "
                << pmfd[i][j].xmax << " "
                << pmfd[i][j].comment << endl; 
        o << endl; 
        return o.str();
      }
  };
}
#endif
