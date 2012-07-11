#include <faunus/inputfile.h>
#include <faunus/textio.h>
#include <faunus/species.h>

namespace Faunus {

  InputMap::InputMap() {}

  InputMap::InputMap(string filename) {
    include(filename);
  }

  /*!
   * \todo Remove ugly atom loader.
   */
  bool InputMap::include(string filename) {
    string line,key,val;
    std::ifstream f( filename.c_str() );
    if (f) {
      incfiles.push_back(filename);
      while (std::getline(f,line)) {
        std::istringstream i(line);
        i >> key >> val;
        if (!key.empty() && !val.empty() ) {
          if (key.find("#")==string::npos) {
            if (val=="yes") val="1";
            if (val=="no") val="0";
            map[key]=val;
          }
        }
      }
      string atomfile = get<string>("atomlist", "");
      if (!atomfile.empty())
        atom.includefile(atomfile);
      return true;
    }
    return false;
  }

  string InputMap::info() {
    short w=25;
    using namespace Faunus::textio;
    std::ostringstream o;
    o << header("User Input Parameters (InputMap)");
    o << pad(SUB,w,"Number of parameters") << map.size() << endl;
    for (auto &m : map)
      o << std::left << setw(25) << m.first << m.second << endl;
    return o.str();
  }

  bool InputMap::save(string file) {
    std::ofstream f(file);
    if (f) {
      f << "# Generated on " << __DATE__ << " " << __TIME__
#ifdef __VERSION__
        << " using " << __VERSION__
#endif
#ifdef __SVN_REV__
        << " (SVN revision: " << __SVN_REV__ << ")" << endl
#endif
      ;
      for (auto &m : map) {
        f << std::left << setw(35) << m.first << setw(20) << m.second;
        if ( !keyinfo[m.first].empty() )
          f << "# " << keyinfo[m.first];
        f << endl;
      }
      return true;
    }
    return false;
  }

  UnitTest::UnitTest(string testfile, bool state) : InputMap(testfile) {
    cnt=0;
    stable=state;
    file=testfile;
  }

  UnitTest::UnitTest(InputMap &in) : InputMap(in.get<string>("test_file","")) {
    cnt=0;
    file=in.get<string>("test_file","");
    stable=in.get<bool>("test_stable", true);
  }

  bool UnitTest::operator()(const string &key, double value, double precision) {
    cnt++;
    if (stable==true) {
      add(key,value);
      save(file);
    } else {
      double ref = get<double>(key, 0.);
      double diff = std::abs( (ref-value)/ref );
      if (diff>precision) {
        failed[key] = std::pair<double,double>(ref, value);
        return false;
      }
    }
    return true;
  }

  string UnitTest::info() {
    short w=37;
    using namespace Faunus::textio;
    std::ostringstream o;
    o << header("Unittests")
      << pad(SUB,w,"Unittest state") << ( (stable==true) ? "stable" : "unstable") << endl
      << pad(SUB,w,"Test file") << file << ( (stable==true) ? " (created)" : " (read)") << endl
      << pad(SUB,w,"Number of tests") << cnt << endl;
    if (stable==false) {
      o << pad(SUB,w,"Number of failed tests") << failed.size() << endl;
      if (!failed.empty()) {
        o << endl << std::left << setw(w+2) << "" << setw(12) << "Stable" << setw(12) << "Current" << "Difference" << endl;
        for (auto &m : failed) {
          double current=m.second.first;
          double ref=m.second.second;
          o << indent(SUBSUB) << std::left << setw(w-2) << m.first
            << setw(12) << m.second.first << setw(12) << m.second.second 
            << std::abs( (ref-current)/ref*100 ) << percent << endl;
        }
      }
    }
    return o.str();
  }

}//namespace
