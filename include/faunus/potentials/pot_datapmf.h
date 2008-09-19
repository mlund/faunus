#ifndef FAU_POT_DATAPMF_H
#define FAU_POT_DATAPMF_H
#include "faunus/potentials/base.h"
#include "faunus/xydata.h"
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
  class pot_datapmf : public pot_lj {
    private:
      xydata<double> pmfd[particle::LAST][particle::LAST];
    public:
      //! \param pot.lB Bjerrum length
      //! \param pot.eps L-J parameter
      pot_datapmf(pot_setup &pot) : pot_lj(pot.eps/pot.lB) {
        f=pot.lB;
        name+="/Empirical data potential";
      }
      bool loadpmf(species &, string);          // load pmf's from disk
      void loadpmf(species &, string,string);   // -//-
      void showpmf(species &);                  //!< Lists loaded pmf's.
      string info();
      double pairpot (const particle &p1, const particle &p2) {
        unsigned short i=p1.id,j=p2.id;
        if (i>j) std::swap(i,j);
        double r2=p1.sqdist(p2);
        if (pmfd[i][j].xmax==0) {               // if no data exists:
          double u,x=p1.charge*p2.charge;       // use Coulomb + lj pot.
          u=lj(p1,p2,r2);
          return (x!=0) ? u+x/sqrt(r2) : u; 
        }; 
        r2=sqrt(r2);
        return (r2>pmfd[i][j].xmax) ? 
          p1.charge*p2.charge/r2 : // use Coulomb pot. outside data 
          pmfd[i][j].x2y(r2);      // ...else use table. 
      }
  };
  /*!\param spc Species class.
   * \param pmfir Directory in which to search for PMF's
   * \param type Particle name to search for. I.e. "NA" or "CL"
   */
  void pot_datapmf::loadpmf(species &spc, string pmfdir, string type) {
    string n1,n2;                                                      
    unsigned short i,j=spc.id(type); 
    for (i=particle::FIRST; i<particle::LAST; i++) { 
      n1=spc.d[i].name+"-"+spc.d[j].name;                                      
      n2=spc.d[j].name+"-"+spc.d[i].name;                                      
      if (loadpmf( spc, pmfdir + n1 + ".dat")==false)                       
        loadpmf( spc, pmfdir + n2 + ".dat");                                
    };                                   
  }
  // Show table of loaded PMF's
  void pot_datapmf::showpmf(species &spc) {
    std::cout << "# --- LOADED PMF's ----------------------------------\n"; 
    std::cout << "# (a,b,resolution,r_max,file)\n"; 
    for (unsigned int i=particle::FIRST; i<particle::LAST; i++)
      for (unsigned int j=particle::FIRST; j<particle::LAST; j++)
        if (pmfd[i][j].xmax>0) 
          std::cout << "# " << spc.d[i].name << " " << spc.d[j].name << " "
            << pmfd[i][j].res << " "
            << pmfd[i][j].xmax << " "
            << pmfd[i][j].comment << endl; 
    std::cout << endl; 
  }
  /*! Load PMF(s) from a file. File format: Each set starts
   * with "#$ type1 type2 length". Several sets can be present
   * in the same file.
   */
  bool pot_datapmf::loadpmf(species &spc, string filename) {
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
}
#endif
