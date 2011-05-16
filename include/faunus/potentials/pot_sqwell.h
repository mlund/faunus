#ifndef FAU_POT_SQWELL_H
#define FAU_POT_SQWELL_H
#include "faunus/potentials/base.h"
#include "faunus/species.h"

namespace Faunus {
  /*!
   * \brief Gereral class for square well liquid with abritrary number of components
   * \author Bjoern Persson
   * \date 2010
   */
  class pot_sqwell {
    private:
    /* This class will use the faunusatoms properties. Note that parameters wd_ij and iw2_ij with
       i or j larger than 2 will always be set to zero. If there is need for more than three components
       with wells/shoulders please adjust this info accordingly before submitting.  
    */
      double hce, box, halfbox;
      vector< vector<double> >
        iw2,    // ij square well interaction dist.
        wd,     // ij well/shoulder depth/height
        hw2;    // ij square mean diameter
    public:
      string name;
      double f;
      pot_sqwell( inputfile &in ) {
        name+="Square well/shoulder particles";
        f=1.0;
        hce=100000.;
        box=in.getflt("boxlen",-1);
        halfbox=box/2;
        init(atom);
        std::cout<<atom.list.size()<<std::endl;
        iw2[1][2]=iw2[2][1]=in.getflt("iw2_21",0);
        iw2[1][3]=iw2[3][1]=in.getflt("iw2_31",0);
        iw2[3][2]=iw2[2][3]=in.getflt("iw2_23",0);
        wd[1][2]=wd[2][1]=in.getflt("wd_21",0);
        wd[1][3]=wd[3][1]=in.getflt("wd_31",0);
        wd[3][2]=wd[2][3]=in.getflt("wd_23",0);
      }
      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist(p2);
      }
      inline double pairpot(const particle &p1, const particle &p2) {
        register double a;
        a=p1.sqdist_mi_xyz(p2,box,halfbox);
        if(a<iw2[p1.id][p2.id])
          return a<hw2[p1.id][p2.id] ? hce : wd[p1.id][p2.id];
        else 
          return 0;
      }
      string info() {
        std::ostringstream o;
        o.precision(3);
        o << "#   Type              = " << name << std::endl
          << "#     Interaction matrix: depth "<< std::endl
          << "#       "<<std::setw(6)<<wd[1][1]<<"  "<<std::setw(6)<<wd[1][2]<<"  "<<std::setw(6)<<wd[1][3]<<std::endl
          << "#       "<<std::setw(6)<<wd[2][1]<<"  "<<std::setw(6)<<wd[2][2]<<"  "<<std::setw(6)<<wd[2][3]<<std::endl
          << "#       "<<std::setw(6)<<wd[3][1]<<"  "<<std::setw(6)<<wd[3][2]<<"  "<<std::setw(6)<<wd[3][3]<<std::endl
          << "#     Interaction matrix: range squared "<< std::endl
          << "#       "<<std::setw(6)<<iw2[1][1]<<"  "<<std::setw(6)<<iw2[1][2]<<"  "<<std::setw(6)<<iw2[1][3]<<std::endl
          << "#       "<<std::setw(6)<<iw2[2][1]<<"  "<<std::setw(6)<<iw2[2][2]<<"  "<<std::setw(6)<<iw2[2][3]<<std::endl
          << "#       "<<std::setw(6)<<iw2[3][1]<<"  "<<std::setw(6)<<iw2[3][2]<<"  "<<std::setw(6)<<iw2[3][3]<<std::endl
          << "#     Hard wall matrix  : range squared "<< std::endl
          << "#       "<<std::setw(6)<<hw2[1][1]<<"  "<<std::setw(6)<<hw2[1][2]<<"  "<<std::setw(6)<<hw2[1][3]<<std::endl
          << "#       "<<std::setw(6)<<hw2[2][1]<<"  "<<std::setw(6)<<hw2[2][2]<<"  "<<std::setw(6)<<hw2[2][3]<<std::endl
          << "#       "<<std::setw(6)<<hw2[3][1]<<"  "<<std::setw(6)<<hw2[3][2]<<"  "<<std::setw(6)<<hw2[3][3]<<std::endl;
        return o.str();
      }
      void init(atoms &a) {
        short i,j,n=a.list.size();
        hw2.resize(n);
        iw2.resize(n);
        wd.resize(n);
        for (i=0; i<n; i++) {
          hw2[i].resize(n);
          iw2[i].resize(n);
          wd[i].resize(n);
          for (j=0; j<n; j++)
          {
          wd[i][j]=0.0;
          //iw2[i][j]=0.0;
          iw2[i][j]=pow(a.list[i].radius+a.list[j].radius,2.0);
          hw2[i][j]=pow(a.list[i].radius+a.list[j].radius,2.0);
          }
        }
      }
      inline double bond(particle &a, particle&b)
       { return 0.0;}

      void setvolume(double v) {}
  };

}//namespace
#endif
