#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/point.h>
#include <faunus/textio.h>
#include <faunus/physconst.h>

namespace Faunus {

  AtomData::AtomData() {
    sigma=0;      
    eps=0;
    radius=0;
    mw=0.1;
    charge=0;
    activity=0;
    dp=0;
    dprot=0;
    mean=0;
    variance=0;
    hydrophobic=0;
    patchtype=0;  
    name="UNK";
  }

  AtomTypes atom; // Instantiate global copy

  AtomTypes::AtomTypes() {
    filename.clear();
    AtomData a;
    a.id=list.size();
    a.name="UNK";
    list.push_back(a);
  }

  AtomData & AtomTypes::operator[] (particle::Tid i) { return list.at(i); }

  AtomData & AtomTypes::operator[] (string s) {
    for (auto &l_i : list)
      if (s==l_i.name)
        return l_i;
    return list.at(0);
  }

  bool AtomTypes::includefile(InputMap &in) {
    return includefile( in.get<string>("atomlist",filename) );
  }

  bool AtomTypes::includefile(string file) {
    AtomData a;
    string t;
    filename=file;
    std::ifstream f(filename.c_str());
    cout << "Reading atom data from '" << filename << "'. ";
    if (f) {
      cout << "OK!\n";
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
      f.close();
      return true;
    }
    cout << "FAILED!\n";
    filename+=" (n/a)";
    return false;
  }

  string AtomTypes::info() {
    using namespace textio;
    char w=25;
    if (filename.empty())
      filename="(undefined)";
    std::ostringstream o;
    o << header("Atomic Species")
      << pad(SUB,w,"Number of species") << list.size() << endl
      << pad(SUB,w,"Parameter file") << filename << endl
      << indent(SUB) << "Species:";
    for (size_t i=0; i<list.size(); i++) {
      if (i%10==0)
        o << endl << indent(SUBSUB);
      o << setw(SUBSUB+1) << std::left << list[i].name;
    }
    o << endl;
    return o.str();
  }

  void AtomTypes::reset_properties(p_vec &p) {
    for (auto &p_i : p )
      p_i = list.at( p_i.id );
  }

}//namespace
