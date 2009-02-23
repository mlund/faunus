#include "faunus/species.h"
namespace Faunus {
  atoms::atoms() {
    filename="faunatoms.dat";
    data a = {0,0,0,0,1.,0,0,false,"UNK"};
    list.push_back(a);
    init();
  }
// dumtest
  void atoms::init() {
    int n=list.size();
    eps.resize(n);
    sigma.resize(n);
    for (int i=0; i<n; i++) {
      eps[i].resize(n);
      sigma[i].resize(n);
      for (int j=0; j<n; j++) {
        eps[i][j]=sqrt(list[i].eps*list[j].eps);
        sigma[i][j]=(list[i].sigma+list[j].sigma)/2;
      }
    }
  }

  particle atoms::get(char i) {
    particle a;
    return set(a,i);
  }

  particle & atoms::set(particle &p, char i) {
    p.charge=list[i].charge;
    p.mw=list[i].mw;
    p.id=static_cast<particle::type>(list[i].id);
    p.radius=list[i].radius;
    return p;
  }

  atoms::data & atoms::operator[] (char i) { return list[i]; }
  atoms::data & atoms::operator[] (string s) { return list[ find(s) ]; }
  particle atoms::operator() (char i) { return get(i); }
  particle atoms::operator() (string s) { return get(find(s)); }

  char atoms::find(string s) {
    for (char i=0; i<list.size(); i++)
      if (s==list[i].name)
        return list[i].id;
    return 0;
  }

  bool atoms::load(inputfile &in) {
    return load(in.getstr("atomfile",filename));
  }

  bool atoms::load(string file) {
    data a;
    string t;
    filename=file;
    std::ifstream f(filename.c_str());
    if (f) {
      while (!f.eof()) {
        f >> t;
        if (t=="Atom") {
          f >> a.name >> a.charge >> a.radius >> a.eps >> a.mw >> a.pka >> t;
          a.sigma=2*a.radius;
          a.hydrophobic = (t=="yes") ? true : false;
          a.id=list.size();
          list.push_back(a);
        }
      }
      init();
      f.close();
      return true;
    }
    std::cerr << "# WARNING! Parameter file " << filename << " NOT found.\n";
    filename+=" (not found)";
    return false;
  }

  string atoms::info() {
    std::ostringstream o;
    o << endl
      << "# ATOM PARAMETERS:\n"
      << "#   Number of species     = " << list.size() << endl
      << "#   Parameter file        = " << filename << endl
      << "#   Species found:";
    for (int i=0; i<list.size(); i++) {
      if (i%10==0) o << endl << "#     ";
      o << setw(6) << std::left << list[i].name;
    }
    o << endl;
    return o.str();
  }

  void atoms::reset_properties(vector<particle> &p) {
    for (int i=0; i<p.size(); i++)
      set(p[i],p[i].id);
  }
}//namespace
