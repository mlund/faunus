#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/textio.h>
#include <faunus/physconst.h>
#include <faunus/json.h>
#include <stdlib.h>


namespace Faunus {

  AtomData::AtomData() {
    activity=0;
    charge=0;
    dp=0;
    dprot=0;
    eps=0;
    hydrophobic=0;
    mean=0;
    muscalar=0;
    mw=1.0;
    name="UNK";
    patchtype=0;  
    radius=0;
    sigma=0;      
    variance=0;
    len=0;
    half_len=0;
    pswitch=0;
    pdis=0;
    pangl=0;
    panglsw=0;
    //rcutwca=0;
    //rcut=0;
    //pcangl=0;
    //pcanglsw=0;
    //pcoshalfi=0;
    //psinhalfi=0;
    chiral_angle=0;
    mu.clear();
    theta.clear();
    alpha.clear();
  }

  AtomMap::AtomMap() {
    AtomData a;
    a.id=list.size();
    a.name="UNK";
    list.push_back(a);
  }

  AtomData & AtomMap::operator[] (AtomData::Tid i) {
    assert(i>=0 && i<list.size() && "Particle id not found!");
    return list.at(i);
  }

  AtomData & AtomMap::operator[] (string s) {
    for (auto &l_i : list)
      if (s==l_i.name)
        return l_i;
    return list.at(0);
  }

  /**
   * This will look for the keyword `atomlist` in the InputMap and
   * use value as the filename.fs
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
      //Tensor<double> fallback;
      //fallback.clear();
      a.name = atom.first;
      a.activity = json::value<double>(atom.second, "activity", 0);
      a.alpha << json::value<std::string>(atom.second, "alpha", "");
      a.alpha *= 4*pc::pi*pc::e0*(1e-10)*pc::kT()/(pc::e*pc::e);
      a.theta << json::value<std::string>(atom.second, "theta", "");
      a.theta *= 0.20819434; // Debye Å -> e Å^2
      a.dp = json::value<double>(atom.second, "dp", 0);
      a.dprot = json::value<double>(atom.second, "dprot", 0) * pc::pi / 180.; // deg->rads
      a.eps = json::value<double>(atom.second, "eps", 0);
      a.hydrophobic = json::value<bool>(atom.second, "hydrophobic", false);
      a.mu << json::value<std::string>(atom.second, "mu", "0 0 0");
      a.muscalar = a.mu.len()*pc::D2eA();
      if (a.mu.len()>1e-6)
        a.mu = a.mu/a.mu.len();
      a.mw = json::value<double>(atom.second, "Mw", 1.);
      a.charge = json::value<double>(atom.second, "q", 0);
      a.radius = json::value<double>(atom.second, "r", 0);
      a.sigma = 2*a.radius;
      a.sigma = json::value<double>(atom.second, "sigma", a.sigma);
      a.radius = a.sigma/2;
      a.id=AtomData::Tid( list.size() );

      a.half_len = 0.5 * json::value<double>(atom.second, "len", 0);
      a.patchtype = json::value<double>(atom.second, "patchtype", 0);
      a.pswitch = json::value<double>(atom.second, "patchswitch", 0);
      a.pdis = json::value<double>(atom.second, "patchdistance", 0);
      a.pangl = json::value<double>(atom.second, "patchangle", 0)/180.0*pc::pi;
      a.panglsw = json::value<double>(atom.second, "patchangleswitch", 0)/180.0*pc::pi;
      a.chiral_angle = json::value<double>(atom.second, "patchchiralangle", 0)/180.0*pc::pi;

      // add to particle list 
      bool insert=true;
      for (auto &i : list)
        if (i.name==a.name) {
          a.id=i.id; // keep old id and
          i=a; // override existing
          insert=false;
          break;
        }
      if (insert)
        list.push_back(a);
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
