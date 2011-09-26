#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/point.h>
#include <faunus/textio.h>

namespace Faunus {
  atoms atom; // Instantiate global copy

  atoms::atoms() {
    filename="atomlist.inp";
    specdata a;
    a.id=list.size();
    a.sigma=0;      
    a.eps=0;
    a.radius=0;
    a.mw=0.1;
    a.charge=0;
    a.chempot=0;
    a.dp=0;
    a.mean=0;
    a.variance=0;
    a.hydrophobic=0;
    a.name="UNK";
    list.push_back(a);
    init();
  }

  void atoms::init() {
    int n=list.size();
    qq.resize(n);
    eps.resize(n);
    sigma.resize(n);
    for (int i=0; i<n; i++) {
      qq[i].resize(n);
      eps[i].resize(n);
      sigma[i].resize(n);
      for (int j=0; j<n; j++) {
        qq[i][j]    = list[i].charge * list[j].charge;
        eps[i][j]   = sqrt(list[i].eps * list[j].eps);
        sigma[i][j] = ( list[i].sigma + list[j].sigma ) / 2;
      }
    }
  }

  specdata & atoms::operator[] (short i) { return list[i]; }

  specdata & atoms::operator[] (string s) {
    for (auto &l_i : list)
      if (s==l_i.name)
        return l_i;
    return list[0];
  }

  bool atoms::includefile(inputfile &in) {
    return includefile(in.getstr("atomfile",filename));
  }

  bool atoms::includefile(string file) {
    specdata a;
    string t;
    filename=file;
    std::ifstream f(filename.c_str());
    if (f) {
      while (!f.eof()) {
        f >> t;
        if (t=="Atom") {
          f >> a.name >> a.charge >> a.radius >> a.eps >> a.mw >> t;
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
    string w="Error! Parameter file " + filename + " not NOT found.",
           id="species";
    std::cerr << "# " << w << endl;
    filename+=" (not found)";
    return false;
  }

  string atoms::info() {
    using namespace textio;
    char w=25;
    std::ostringstream o;
    o << header("Atomic Species")
      << pad(SUB,w,"Number of species") << list.size() << endl
      << pad(SUB,w,"Parameter file") << filename << endl
      << indent(SUB) << "Species:";
    for (int i=0; i<list.size(); i++) {
      if (i%10==0)
        o << endl << indent(SUBSUB);
      o << setw(SUBSUB) << std::left << list[i].name;
    }
    o << endl;
    return o.str();
  }

  void atoms::reset_properties(vector<particle> &p) {
    for (auto &p_i : p )
      p_i = list[p_i.id];
  }

}//namespace
