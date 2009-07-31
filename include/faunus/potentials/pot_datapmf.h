#ifndef FAU_POT_DATAPMF_H
#define FAU_POT_DATAPMF_H

#include "faunus/species.h"
#include "faunus/xydata.h"
#include "faunus/potentials/base.h"
#include "faunus/potentials/pot_hscoulomb.h"

namespace Faunus {
  /*!
   * \brief Load pair potentials from disk
   * \author Mikael Lund
   * \date Prague, 2007
   *
   * Class to load potential of mean force (PMF) data from file(s). If
   * the distance between the particles is outside the data range a 
   * simple Coulomb/Hardsphere potential will be applied.
   */
  class pot_datapmf : public pot_hscoulomb {
    private:
      //xydata<double> pmfd[particle::LAST][particle::LAST];
      xydata<double> pmfd[20][20];
    public:
      string pmfdir; //!< Directory containing PMF data
      pot_datapmf(inputfile &in) : pot_hscoulomb(in) {
        name+="/Empirical data potential";
        pmfdir=in.getstr("pmfdir","./");
      }
      double pairpot (const particle &p1, const particle &p2) {
        unsigned short i=p1.id,j=p2.id;
        if (i>j)
          std::swap(i,j);
        if (pmfd[i][j].xmax<0.01)
          return pot_hscoulomb::pairpot(p1,p2);
        double r2=p1.sqdist(p2);
        return ( r2 > pmfd[i][j].xmax * pmfd[i][j].xmax ) ?
          pot_hscoulomb::pairpot(p1,p2) :  // use Coulomb pot. outside data
          pmfd[i][j].x2y(sqrt(r2));        // ...else use table.
      }

      /*!
       * \param id particle type to search for.
       */
      void loadpmf(unsigned char id) {
        string n1,n2;
        char i,j=id;
        for (i=0; i<atom.list.size(); i++) {
          n1=atom[i].name+"-"+atom[j].name;
          n2=atom[j].name+"-"+atom[i].name;
          if (loadpmf(n1 + ".dat")==false)
            loadpmf(n2 + ".dat");
        }
      }

      /* Search through particle vector and attempt to load pmf's for all possible
       * particle combinations.
       */
      void loadpmf(container &c) {
        vector<unsigned char> id;
        for (unsigned short i=0; i<c.p.size(); i++) {
          vector<unsigned char>::iterator iter = std::find(id.begin(), id.end(), c.p[i].id);
          if (iter==id.end())
            id.push_back(c.p[i].id);
        }
        for (unsigned short i=0; i<id.size(); i++)
          for (unsigned short j=i; j<id.size(); j++)
            loadpmf(atom[id[i]].name+"-"+atom[id[j]].name+".dat"); 
      }

      /*! Load PMF(s) from a file. File format: Each set starts
       * with "#$ type1 type2 length". Several sets can be present
       * in the same file.
       */
      bool loadpmf(string filename) {
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
              a = atom[a_str].id;
              b = atom[b_str].id;
              if (a>b)
                std::swap(a,b);
              for (int i=0; i<len; i++) {
                fh >> x[i] >> y[i];
                y[i]=y[i]/f; // unit: kT/f
              }
              if (pmfd[a][b].xmax==0) {
                pmfd[a][b].add( x, y );
                pmfd[a][b].comment=filename;
              }
            }
          }
          fh.close();
          return true;
        }
        return false;
      }

      /*! Show info + list of loaded pmf data */
      string info() {
        double n=atom.list.size();
        std::ostringstream o;
        o << pot_hscoulomb::info()
          << "#   PMF directory     = " << pmfdir << endl
          << "#   PMF Info: (a,b,resolution,r_max,file)" << std::endl; 
        for (unsigned short i=0; i<n; i++)
          for (unsigned short j=0; n; j++)
            if (pmfd[i][j].xmax>0) 
              o << "#     " << atom[i].name << " " << atom[j].name << " "
                << pmfd[i][j].res << " "
                << pmfd[i][j].xmax << " "
                << pmfd[i][j].comment << endl; 
        o << endl; 
        return o.str();
      }
  };
}
#endif
