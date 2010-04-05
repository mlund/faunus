#ifndef FAUNUS_HYDROPHOBIC_H
#define FAUNUS_HYDROPHOBIC_H

namespace Faunus {
  
  /*!
   * \brief Hydrophobic interaction between ions and molecular surfaces
   * \author Mikael Lund
   * \date Canberra 2008
   * \todo Not optimized - inelegant "end_of_protein_one" hack. Use vector instead...
   *
   * This class will use the specified pair potential as usual but in addition add
   * a hydrophobic interaction between ions (or any specified species) and the
   * nearest hydrophobic particle. This requires an expanded pair-potential that
   * contains a function hypairpot(). If you need the ions to interact with the
   * hydrophobic groups on TWO proteins, you must set the end_of_protein_one variable.
   * In this way the minimum distance search is repeated on the remaining particles.
   *
   * \warning Calling an energy() function means that hyenergy() is added more
   *          than once. Carefully check your MC moves to see if this is the
   *          case.
   */
  template<class T> class int_hydrophobic : public interaction<T> {
  private:
    vector<unsigned short> hy; //!< List of hydrophobic particles
    vector<unsigned short> pa; //!< List of particles with empirical pmf.
    double hyenergy(vector<particle> &);
    double hyenergy(vector<particle> &, int);
    
  public:
    unsigned int end_of_protein_one;     //!< Last particle in protein one (set if appropriate)
    
    int_hydrophobic(inputfile &in) : interaction<T>(in) {
      end_of_protein_one=int(1e7);
    }
    
    void search(vector<particle> &);     //!< Locate hydrophobic groups and ions
    
    double energy(vector<particle> &p) {
      return interaction<T>::energy(p) + hyenergy(p);
    }
    
    double energy(vector<particle> &p, int i) {
      return interaction<T>::energy(p,i) + hyenergy(p,i);
    }
    
    double energy(vector<particle> &p, const group &g ) {
      return interaction<T>::energy(p,g) + hyenergy(p);
    }
    
    double energy(vector<particle> &p, const group &g1, const group &g2 ) {
      return interaction<T>::energy(p,g1,g2) + hyenergy(p); //
    }
  };
  
  template<class T> void int_hydrophobic<T>::search(vector<particle> &p) {
    pa.clear();
    hy.clear();
    for (int i=0; i<p.size(); i++)
      if (p[i].hydrophobic==true)
        hy.push_back(i);
      else if (p[i].id==atom["NA"].id ||
               p[i].id==atom["CL"].id ||
               p[i].id==atom["I"].id) pa.push_back(i);
    
    cout << "# Hydrophobic sites = " << hy.size() << endl
    << "# Ions              = " << pa.size() << endl;
  }
  
  template<class T> double int_hydrophobic<T>::hyenergy(vector<particle> &p) {
    double u=0;
    int n=pa.size();
#pragma omp parallel for reduction (+:u)
    for (int i=0; i<n; i++)     // loop over ions
      u+=hyenergy(p, pa[i]);    // energy with hydrophobic groups
    return u; // in kT
  }
  
  template<class T> double int_hydrophobic<T>::hyenergy(vector<particle> &p, int i) {
    if (p[i].hydrophobic==true)
      return 0;
    int hymin=0;
    double d,dmin=1e7,u=0;
    for (int j=0; j<hy.size(); j++) { // loop over hydrophobic groups
      if (hy[j]>end_of_protein_one) { // test if we move into second protein
        u=interaction<T>::pair.hypairpot( p[i], p[hymin], sqrt(dmin) );
        dmin=1e7;                     // reset min. distane
      }
      d=p[i].sqdist( p[hy[j]]);       // find min. distance
      if (d<dmin) {
        dmin=d;      // save min dist.
        hymin=hy[j]; // ...and particle number
      }
    }
    return interaction<T>::pair.f *
    (u + interaction<T>::pair.hypairpot( p[i], p[hymin], sqrt(dmin) ) );
  }
} //namespace
#endif