#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/point.h>
#include <faunus/textio.h>

namespace Faunus {

  AtomData::AtomData() {
    sigma=0;      
    eps=0;
    radius=0;
    mw=0.1;
    charge=0;
    chempot=0;
    dp=0;
    mean=0;
    variance=0;
    hydrophobic=0;
    name="UNK";
  }

  AtomTypes atom; // Instantiate global copy

  AtomTypes::AtomTypes() {
    filename="atomlist.inp";
    AtomData a;
    a.id=list.size();
    a.name="UNK";
    list.push_back(a);
    init();
  }

  void AtomTypes::init() {
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

  AtomData & AtomTypes::operator[] (short i) { return list[i]; }

  AtomData & AtomTypes::operator[] (string s) {
    for (auto &l_i : list)
      if (s==l_i.name)
        return l_i;
    return list[0];
  }

  bool AtomTypes::includefile(InputMap &in) {
    return includefile(in.get<string>("atomfile",filename));
  }

  bool AtomTypes::includefile(string file) {
    AtomData a;
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

  string AtomTypes::info() {
    using namespace textio;
    char w=25;
    std::ostringstream o;
    o << header("Atomic Species")
      << pad(SUB,w,"Number of species") << list.size() << endl
      << pad(SUB,w,"Parameter file") << filename << endl
      << indent(SUB) << "Species:";
    for (size_t i=0; i<list.size(); i++) {
      if (i%10==0)
        o << endl << indent(SUBSUB);
      o << setw(SUBSUB) << std::left << list[i].name;
    }
    o << endl;
    return o.str();
  }

  void AtomTypes::reset_properties(vector<particle> &p) {
    for (auto &p_i : p )
      p_i = list[p_i.id];
  }

}//namespace
