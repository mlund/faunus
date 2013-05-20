#include <faunus/io.h>
#include <faunus/species.h>
#include <faunus/geometry.h>
#include <faunus/point.h>
#include <faunus/group.h>
#include <faunus/inputfile.h>
#include <faunus/space.h>
#include <faunus/energy.h>

namespace Faunus {
  
  particle FormatGRO::s2p(string &s) {
    std::stringstream o;
    string name;
    double x,y,z;
    o << s.substr(10,5) << s.substr(20,8) << s.substr(28,8) << s.substr(36,8);
    o >> name >> x >> y >> z;
    particle p;
    p=atom[name]; 
    p.x()=x*10; // nm->angstrom
    p.y()=y*10;
    p.z()=z*10;
    return p;
  }

  bool FormatGRO::load(string file) {
    p.clear();
    v.resize(0);
    if (IO::readFile(file,v)==true) {
      int last=atoi(v[1].c_str())+1;
      for (int i=2; i<=last; i++)
        p.push_back( s2p(v[i]) );
      return true;
    }
    return false;
  }

  //----------------- IOQTRAJ ----------------------
  FormatTopology::FormatTopology() : rescnt(0) {}
  //bool FormatTopology::open(string);                         //!< Open file for writing
  //void FormatTopology::writeDefaults();
  //
  //name     mass        charge    ptype         sigma           epsilon
}  //namespace
