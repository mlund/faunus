#ifndef FAU_ENERGY_RF_H
#define FAU_ENERGY_RF_H

#include "faunus/energy/base.h"
#include "faunus/potentials/pot_test_rf.h"

namespace Faunus {
  /*!
   * \brief .... 
   * \todo Not finished
   * \warning Untested!
   * \author Mikael Lund and Bjorn Persson
   * \date Lund 2009
   *
   * some description
   */
  class reactionfield : public interaction<pot_test_rf> {
  private:
  public:
    reactionfield(inputfile &in) : interaction<pot_test_rf>(in) {
      reactionfield::name+="SOMETHING!";
    }
    // TOTAL ENERGY
    double potential(const vector<particle> &p, int i) {
      double phi=0;
      for (int j=0; j<i; j++)
        phi += pair.rf.phi(p[j], p[i]);
      for (int j=0; j<i; j++)
        phi += pair.rf.phi(p[j], p[i]);
      return pair.rf.f*phi;  //*lB ??
    }
    double energy(vector<particle> &p ) {
      double u=interaction<pot_test_rf>::energy(p);
      for (int i=0; i<p.size(); i++)
        u+=pair.rf.f*pair.rf.selfenergy(p[i]);
      return u;
    }
    double energy(vector<particle> &p, int i) {
      return pair.rf.f*pair.rf.selfenergy(p[i]) + interaction<pot_test_rf>::energy(p,i);
    }
    double energy(vector<particle> &p, const group &g ) {
      double u=interaction<pot_test_rf>::energy(p,g);
      for (int i=g.beg; i<=g.end; i++)
        u += pair.rf.selfenergy(p[i]);
      return pair.rf.f*u;
    }
    string info() {
      std::ostringstream o;
      o << interaction<pot_test_rf>::info()
      << pair.info();
      return o.str();
    }
  };
  
  /*!
   * \brief Interaction class that includes image charges outside a spherical cell
   * \todo Not optimized
   * \warning Untested! (behaves properly but not independently confirmed correct)
   *          This class is only INTENDED to be used with the potetnential test_pot!
   * \author Mikael Lund, Bjorn Persson
   * \date Lund 2008-2009
   *
   * Calculates inter-particle interactions in a dielectric sphere (epsi) and
   * adds image charge interactions with the dielectric surroundings (epso). Useful
   * for simulating explicit water in a spherical container. Note that after each
   * MC move -- both accepted and rejected -- the img vector *must* be updated
   * with the updateimg() function. The cavity origin is assumed to be 0,0,0. Note that
   * image charge and distance are divergent at (0,0,0). Also keep in mind that the energy
   * is NOT the product of the potential at a charge times the charge it self but rather 
   * half of the product since it is an reaction potential (as what is conserned for
   * the reaction part). There are there for additional 
   * functions sphericalimage::elenergy to obtain the electrostatic energy for particles 
   * and groups. The functions sphericalimage::potential returns the total potential.
   *
   * For reference: Harold L Friedman, Molecular Physics, 1975 vol. 29 ppg. 1533-1543.
   * Extended method: Wei Cai, XXXX.
   */
  template<class T> class sphericalimage : public interaction<T> {
  private:
    average<double> ratio;
    double scale, radius, radius2, epso, epsi, ui, ur;
    vector<point> img;   //!< Contain image charge particles
    vector<double> ich;  //!< Complementary vector for image charges
  public:
    sphericalimage(inputfile &in) : interaction<T>(in) {
      interaction<T>::name+=" w. spherical image charges";
      epso=in.getflt("epso",80);   //permitivity of outside medium
      epsi=in.getflt("epsi",1);  //permitivity of inner medium
      radius=in.getflt("cellradius",0); 
      radius2=radius*radius;
      scale=-(epso-epsi)/(epso+epsi)*radius/2*in.getflt("bjerrum", 560.2);  // Friedman 1975, Mol. Phys 29 pp. 1533
    }                                                                 // warning!!! temperature units are given through scale!
    
    // RE-CALC IMAGE POSITIONS
    inline void updateimg(const particle &a, int i) {
      img[i] = a* (radius2 / a.dot(a)); 
      ich[i] = a.charge/a.len();
    }
    void updateimg(vector<particle> &p) {
      img.resize( p.size() );
      ich.resize( p.size() );
      for (int i=0; i<img.size(); ++i)
        updateimg(p[i],i);
    }
    void updateimg(vector<particle> &p, const group &g) {
      for (int i=g.beg; i<=g.end; ++i)
        updateimg(p[i],i);
    }
    void updateimg(vector<particle> &p, molecules &m, vector<int> &n) {
      for (int i=0; i<n.size(); i++) {
        updateimg(p,m[i]);
        updateimg(p,m[n[i]]);
      }
    }
    // IMAGE POTENTIAL  (half the potential)
    double impot(double &ich, const point &pr, point &pi) {
      return ich / sqrt(pr.sqdist(pi));
    }
    
    // IMAGE ENERGY
    double image(vector<particle> &p) {
      int n=p.size();
      double u=0;
#pragma omp parallel for reduction (+:u) schedule (dynamic)
      for (int i=0; i<n; ++i)
        u+=image(p,i);
      return u;
    }
    double image(vector<particle> &p, int i) {
      double u=0;
      int t=img.size();
#pragma omp parallel for reduction (+:u) schedule (dynamic)
      for (int j=0; j<t; ++j)
        u += impot(ich[j], p[i], img[j] );
      return p[i].charge*scale*u*2. - p[i].charge*scale*impot(ich[i],p[i],img[i]) ;        //make sure too and not too double count
    }
    double image(vector<particle> &p, const group &g) {
      double u=0;
      for (int i=g.beg; i<=g.end; ++i)
        u+=image(p, i, g.beg, g.end);
      return u;
    }
    double image(vector<particle> &p, int i, int j, int k) {
      int n=p.size();
      double uin=0,uex=0;
#pragma omp parallel for reduction (+:uex) schedule (dynamic)
      for (int s=0; s<j; s++)
        uex += impot(ich[s], p[i], img[s] );  //make sure to double count
#pragma omp parallel for reduction (+:uex) schedule (dynamic)
      for (int t=k+1; t<n; t++)
        uex += impot(ich[t], p[i], img[t] );  //make sure to double count
      for (int u=j; u<=k; u++)
        uin += impot(ich[u], p[i], img[u] );  // internal interactions will be double counted implicitly
      // not paralleized on purpose
      // the self term will not be double counted
      return p[i].charge*scale*(uex*2+uin); 
    }
    double image(vector<particle> &p, const group &g1, const group &g2) {
      double u=0;
      for (int i=g1.beg; i<=g1.end; i++)               //Dielectric and g1
        for (int j=g1.beg; j<=g1.end; j++)
          u += p[i].charge*impot(ich[j], p[i], img[j]);
      for (int i=g2.beg; i<=g2.end; i++)               //Delectric and g2
        for (int j=g2.beg; j<=g2.end; j++)
          u += p[i].charge*impot(ich[j], p[i], img[j]);
      for (int i=g1.beg; i<=g1.end; i++)               //g1 and g2 through dielectric
        for (int j=g2.beg; j<=g2.end; j++)
          u += 2*p[i].charge*impot(ich[j], p[i], img[j]);
      return u*scale;
    } 
    double imageint(vector<particle> &p, group) {
      double u=0;
      return u;
    }
    // TOTAL ENERGY
    double potential(vector<particle> &p, int i) {
      int n=p.size();
      double ur=0,ui=0;
      updateimg(p[i],i);
      ur=interaction<T>::potential(p,i);
#pragma omp parallel for reduction (+:ui) schedule (dynamic)
      for (int s=0; s<n; s++)
        ui += impot(ich[s], p[i], img[s] );  
      ui*=2;
      return ur+ui*scale;  // This is the POTENTIAL, this should not be used to calculate the 
    }                      // interaction since that would 'double' count the self term.
    double potential(vector<particle> &p, point i) {
      ur=ui=0;
      ur=interaction<T>::potential(p,i);
      //#pragma omp parallel for reduction (+:ui) schedule (dynamic)
      for (int s=0; s<p.size(); s++)
        ui += impot(ich[s], i, img[s] );  
      ui*=2;               // Due to the definition of scale
      return ur+ui*scale;  // This is the POTENTIAL in any given point inside the cavity
    }                     
    double elenergy(vector<particle> &p, int i) {
      ur=ui=0;
      updateimg(p[i],i);
      ur=p[i].charge*interaction<T>::potential(p,i);
      ui=image(p,i);
      return ur+ui;        // Returns the electrostatic interaction energy of i with p
    }
    double elenergy(vector<particle> &p, const group &g) {
      ur=ui=0;
      updateimg(p,g);
      for (int i=g.beg; i<=g.end; i++) {
        ur+=p[i].charge*interaction<T>::potential(p,i);
        ui+=image(p,i);
      }
      return ur+ui;        // Returns the electrostatic energy of g with p
    }
    double energy(vector<particle> &p ) {
      updateimg(p);
      ur=interaction<T>::energy(p);
      ui=image(p);
      ratio+=std::abs(ui/(ur+ui));
      return ur+ui;
    }
    double energy(vector<particle> &p, int i) {
      updateimg(p[i],i);
      ur=interaction<T>::energy(p,i);
      ui=image(p,i);
      ratio+=std::abs(ui/(ur+ui));
      return ur+ui;
    }
    double energy(vector<particle> &p, const group &g ) {
      updateimg(p,g); // update all images in group
      ur=interaction<T>::energy(p,g);
      ui=image(p,g);
      ratio+=std::abs(ui/(ur+ui));
      return ur+ui;
    }
    double energy(vector<particle> &p, const group&g1, const group &g2) {
      updateimg(p,g1), updateimg(p,g2);
      ur=interaction<T>::energy(p, g1, g2);
      ui=image(p,g1,g2);
      return ur+ui;
    }
    double energy(vector<particle> &p, molecules &m, vector<int> &n) {
      ur=ui=0;
      group g;
      g.beg=m[0].beg;
      g.end=m[n.size()-1].end;
      updateimg(p,m,n);
      ur=interaction<T>::energy(p,g);
      for (int i=0; i<n.size()-1; i++)
        for (int j=i+1; j<n.size(); j++)
          ur+=interaction<T>::energy(p,m[i],m[j]);
      //for (int i=0; i<n.size(); i++)
      //  ur+=interaction<T>::energy(p, m[i]);
      ui=image(p,g);
      ratio+=std::abs(ui/(ur+ui));
      return ur+ui;
    }
    // INFO
    string info() {
      std::ostringstream o;
      o << interaction<T>::info()
      << "#   Dielectric const. (in, out) = " << epsi << " " << epso << endl
      << "#   Cavity radius               = " << radius << endl
      << "#   Avg. image energy ratio     = " << ratio.avg() <<" , stdev "<< ratio.stdev()<< endl
      << "#   Number of images            = " << img.size() << endl
      << "#   Scaling const.              = " << scale << endl
      << "#      -(epso-1)/(epso+1)*radius/2*lb "<<endl;
      return o.str();
    }
    string printimg() {
      std::ostringstream o;
      for (int i =0; i<img.size(); i++) 
        o << img[i] <<"  "<<ich[i]<<endl;
      o << endl;
      return o.str();
    }
  };
  
}//namespace
#endif
