#ifndef FAU_INPUTFILE_H
#define FAU_INPUTFILE_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/textio.h>
#include <faunus/species.h>
#include <faunus/json.h>
#endif

namespace Faunus {

  /**
   * @brief Retrieve parameters from a formatted input file
   *
   * The input file must have this format:
   *
   * -  `KEYWORD` `VALUE` (space separated)
   *
   * Blank lines or lines containing `[` or `#` are ignored.
   *
   * Example:
   *
   *     InputMap in("myinput.in");
   *     double T = in.get<double>("temperature", 298.15);
   */
  class InputMap {
    private:
      std::map<string,string> map, keyinfo;
      vector<string> usedkeys;
      vector<string> incfiles;
      bool isjson;
      json::Tval js;
      json::Tval jsdir;

    public:
      InputMap();
      InputMap(string);      //!< Construct and include input file.

      bool include(string);  //!< Include keyword/value file.
      bool save(string);     //!< Save map to disk
      string info();         //!< Information string about read files and keywords

      /**
       * @brief CD (change dir) like functionality
       * 
       * Example:
       *
       * ~~~~
       *   InputMap in( "input.json" );
       *   in.cd( "general/system" );
       *   double T = in.get( "temperature", 298.15 );
       *   in.cd(); // back to upper level (root)
       * ~~~~
       *
       * where the json file may look like this,
       *
       * ~~~~
       *   {
       *     "general" : {
       *        "system" : { "temperature":341. },
       *        "other subsection" : { ... }
       *     }
       *   }
       * ~~~~
       */
      inline void cd( std::string dir="" ) {
        if ( isjson ) { 
          if ( dir.empty() ) {
            jsdir = js;
          } else {
            std::replace( dir.begin(), dir.end(), '/', ' ');
            jsdir = json::cd( js, textio::words2vec<std::string>( dir ) ); 
          }
        }
      } 

      inline json::Tobj getMap() {
        if (isjson)
          if (jsdir.is<json::Tobj>())
            return jsdir.get<json::Tobj>();
        return json::Tobj();
      }

      //!< Add a keyword and an associated value
      template<typename T>
        void add(const string &key, T value, string infostring=string()) {
          std::ostringstream o;
          o << value;
          map[key]=o.str();
          if (!infostring.empty())
            keyinfo[key]=infostring;
        }

      //!< Get value associated with keyword
      template<typename T>
        T get(const string &key, T fallback=T(), string infostring=string()) {
          if (isjson) {
            typedef typename std::conditional< std::is_same<int, T>::value,
                    double,T>::type Tjson;
            return json::value<Tjson>( jsdir, key, fallback );
          }

          if ( !infostring.empty() )
            keyinfo[key] = infostring;               // save information string (if any)
          if ( map.find(key) != map.end() ) {
            std::istringstream i( map[key] );
            i >> fallback;
          }
          return fallback;
        }

      //!< Set value associated with keyword
      template<typename T>
        void set(const string &key, T value, string infostring=string()) {
          if ( map.find(key) != map.end() ) {
            std::ostringstream o;
            o << value;
            map[key] = o.str();
          } else {
            add(key,value,infostring);
          }
        }

      //!< Get value associated with keyword
      template<typename T>
        T operator()(const string &key, T fallback=T()) {
          return get<T>(key,fallback);
        }
  };

  inline InputMap::InputMap() : isjson(false) {}

  inline InputMap::InputMap(string filename) : isjson(false) { include(filename); }

  /** @todo Remove ugly atom loader */
  inline bool InputMap::include(string filename) {
    string line,key,val;

    isjson = std::regex_match( filename,
        std::regex( "(.+?)(\\.json)", std::regex_constants::icase) );

    if ( isjson ) {
      jsdir = js = json::open( filename );
      assert( js != json::Tval() && "Error loading json file" );
    }

    if ( !isjson ) {
      std::ifstream f( filename.c_str() );
      if (f) {
        incfiles.push_back(filename);
        while (std::getline(f,line)) {
          std::istringstream i(line);
          i >> key >> val;
          if (!key.empty() && !val.empty() ) {
            if (key.find("#")==string::npos) {
              if (val=="yes" || val=="true") val="1";
              if (val=="no" || val=="false") val="0";
              map[key]=val;
            }
          }
        }
      } else return false;
    }

    cd("general");
    string atomfile = get<string>("atomlist", "");
    if ( !atomfile.empty() )
      atom.includefile(atomfile);

    if ( atom.empty() )
     std::cerr << "Warning: AtomMap is empty." << endl; 

    return true;
  }

  inline string InputMap::info() {
    short w=25;
    using namespace Faunus::textio;
    std::ostringstream o;
    o << header("User Input Parameters (InputMap)");
    o << pad(SUB,w,"Number of parameters") << map.size() << endl;
    for (auto &m : map)
      o << std::left << setw(25) << m.first << m.second << endl;
    return o.str();
  }

  inline bool InputMap::save(string file) {
    std::ofstream f(file);
    if (f) {
      f << "# Generated on " << __DATE__ << " " << __TIME__
#ifdef __VERSION__
        << " using " << __VERSION__
#endif
        << endl;
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

  /**
   * @brief Class for checking generated output against stored data
   *
   * This class is used to test generated output by comparing against
   * a database file generated for a "stable" run, known to work -
   * a reference system.
   * The input file passed to the constructor is expected to have the following
   * keywords:
   *
   * - `test_stable` Boolean for state. Default: NO = `unstable`
   * - `test_file`   Name of test file to load (if `unstable`) or
   *                 generate (if `stable`).
   *
   * To check a value, simply call unittest(name,value,precision%).
   * If `stable==true` this will generate a new reference testfile.
   * If `stable==false` this will check the given value against the
   * loaded testfile.
   *
   */
  class UnitTest : private InputMap {
    private:
      std::map<string, std::pair<double, double> > failed;
      int cnt;
      string file; // test file
      bool stable; //!< True if passed value is OK to be saved to disk
    public:
      UnitTest(string, bool=true);  //!< Constructor
      UnitTest(InputMap&);          //!< Constructor using InputMap
      bool operator()(const string&, double, double=0.1); //!< Check or set value
      string info(); //!< Info on passed and failed tests
      int numFailed() const; //!< Number of failed tests - i.e. zero if flawless
  };

  inline UnitTest::UnitTest(string testfile, bool state) : InputMap(testfile) {
    cnt=0;
    stable=state;
    file=testfile;
  }

  inline UnitTest::UnitTest(InputMap &in) {
    cnt=0;
    in.cd ("general/unittest");
    file = in.get<string>( "testfile","" );
    stable = in.get<bool>( "stable", true );
    include( file );
  }

  inline bool UnitTest::operator()(const string &key, double value, double precision) {
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

  inline int UnitTest::numFailed() const { return (int)failed.size(); }

  inline string UnitTest::info() {
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
          double current=m.second.second;
          double ref=m.second.first;
          o << indent(SUBSUB) << std::left << setw(w-2) << m.first
            << setw(12) << ref << setw(12) << current
            << std::abs( (ref-current)/ref*100 ) << percent << endl;
        }
      }
    }
    return o.str();
  }
}//namespace
#endif

