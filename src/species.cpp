#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/textio.h>
#include <faunus/physconst.h>
#include <faunus/json.h>

namespace Faunus {

  AtomData::AtomData() {
    activity=0;
    alpha=0;
    charge=0;
    dp=0;
    dprot=0;
    eps=0;
    hydrophobic=0;
    mean=0;
    mu=0;
    mw=1.0;
    name="UNK";
    patchtype=0;  
    radius=0;
    sigma=0;      
    variance=0;
  }

  AtomMap::AtomMap() {
    AtomData a;
    a.id=list.size();
    a.name="UNK";
    list.push_back(a);
  }

  AtomData & AtomMap::operator[] (AtomData::Tid i) { return list.at(i); }

  AtomData & AtomMap::operator[] (string s) {
    for (auto &l_i : list)
      if (s==l_i.name)
        return l_i;
    return list.at(0);
  }

  /**
   * This will look for the keyword `atomlist` in the InputMap and
   * use value as the filename.
   */
  bool AtomMap::includefile(InputMap &in) {
    return includefile( in.get<string>("atomlist",filename) );
  }

  bool AtomMap::includeJSON(const string& file) {
    int n=0;
    filename=file;
    auto j=json::open(file);
    for (auto &atom : json::object("atomlist", j)) {
      n++;
      AtomData a;
      a.name = atom.first;
      a.activity = json::value<double>(atom.second, "activity", 0);
      a.alpha = json::value<double>(atom.second, "alpha", 0)*4*pc::pi*pc::e0*(1e-10)*pc::kT()/(pc::e*pc::e); // Get in units of (4\pi \epsilon_0)Å^3
      a.alphamatrix = Eigen::Matrix3d::Identity()*a.alpha;
      a.dp = json::value<double>(atom.second, "dp", 0);
      a.dprot = json::value<double>(atom.second, "dprot", 0) * pc::pi / 180.; // deg->rads
      a.eps = json::value<double>(atom.second, "eps", 0);
      a.hydrophobic = json::value<bool>(atom.second, "hydrophobic", false);
      a.mu = json::value<double>(atom.second, "mu", 0)*pc::D2eA(); // Debye to eÅ
      a.mw = json::value<double>(atom.second, "Mw", 1.);
      a.charge = json::value<double>(atom.second, "q", 0);
      a.radius = json::value<double>(atom.second, "r", 0);
      a.sigma = 2*a.radius;
      a.sigma = json::value<double>(atom.second, "sigma", a.sigma);
      a.radius = a.sigma/2;
      a.id=AtomData::Tid( list.size() );
      a.patchtype = json::value<double>(atom.second, "patchtype", 0);

      list.push_back(a); // add to main particle list
    }
    return (n>0) ? true : false;
  }

  /**
   * This will automatically detect if the file is a JSON
   * file and call includeJSON(). If not, the old restricted
   * file format is used.
   */
  bool AtomMap::includefile(string file) {
    // is it a JSON file?
    if (file.substr(file.find_last_of(".") + 1) == "json")
      return includeJSON(file);
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

  string AtomMap::info() {
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

  AtomMap atom; // Instantiate global copy

}//namespace
