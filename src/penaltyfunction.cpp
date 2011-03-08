#include "faunus/energy/penaltyfunction.h"

namespace Faunus {

penaltyfunction::penaltyfunction(double xmin, double xmax, double xres, double upenalty, double scalingfactor)
  : v(xres,xmin,xmax), cnt(xres,xmin,xmax)
  {
    number_of_updates=0;
    du=upenalty;
    sf=scalingfactor;
    for (int i=0; i<v.y.size(); i++) {
      v.y.at(i)=0;
      cnt.y.at(i)=0;
    }
  }
  
  double penaltyfunction::update(double x) {
    number_of_updates++;
    cnt(x)++;
    v(x)+=du;
    return du;
  }

  double penaltyfunction::energy(double x) {
    return v(x);
  }
  
  void penaltyfunction::write(string file) {
    std::ofstream f(file.c_str());
    if (f) {
      f.precision(6);
      for (double x=v.xmin; x<v.xmax(); x+=v.xres) {
        f << x << " " << -v(x) << " " << -std::log( cnt(x)+1e-4 ) << endl;
      }
      f.close();
    }
  }
  
  void penaltyfunction::dump(string filename) {
    v.dumptodisk(filename);
  }
  
  void penaltyfunction::dump(string file, int num, string ext) {
    v.dumptodisk(file, num, ext);
  }
  
  void penaltyfunction::gofrdump(string filename) {
    cnt.dumptodisk(filename);
  }

  void penaltyfunction::gofrdump(string file, int num, string ext) {
    cnt.dumptodisk(file, num, ext);
  }

  bool penaltyfunction::load(string filename) {
    if (v.loadfromdisk(filename)==true) {
      cout << "# Loaded penaltyfunction from " << filename << endl;
      return true;
    } else {
      cout << "# No penaltyfunction loaded from disk" << endl;
      return false;
    }
  }
  
  bool penaltyfunction::gofrload(string filename) {
    if (cnt.loadfromdisk(filename)==true) {
      cout << "# Loaded gofr from " << filename << endl;
      return true;
    } else {
      cout << "# No gofr loaded from disk" << endl;
      return false;
    }
  }  
  
  double penaltyfunction::scaledu() {
    du=du*sf;
    return du;
  }
   
  string penaltyfunction::info() {
    std::ostringstream o;
    o << "#  Penalty function:\n"
      << "#     Penalty energy (kT) = " << du << std::endl
      << "#     Energy scaling factor = " << sf << std::endl
      << "#     Number of updates   = " << number_of_updates << std::endl;
    return o.str();
  }
}//namespace

