#ifndef FAU_INPUTFILE_H
#define FAU_INPUTFILE_H

#ifndef SWIG
#include <faunus/common.h>
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
    public:
      InputMap();
      InputMap(string);      //!< Construct and include input file.
      bool include(string);  //!< Include keyword/value file.
      bool save(string);     //!< Save map to disk
      string info();         //!< Information string about read files and keywords
      //!< Add a keyword and an associated value
      template<typename T> void add(const string &key, T value, string infostring=string()) {
        std::ostringstream o;
        o << value;
        map[key]=o.str();
        if (!infostring.empty())
          keyinfo[key]=infostring;
      }

      //!< Get value associated with keyword
      template<typename T> T get(const string &key, T fallback=T(), string infostring=string()) {
        if ( !infostring.empty() )
          keyinfo[key] = infostring;               // save information string (if any)
        if ( map.find(key) != map.end() ) {
          std::istringstream i( map[key] );
          i >> fallback;
        }
        return fallback;
      }
      
      //!< Get value associated with keyword
      template<typename T>
      T operator()(const string &key, T fallback=T()) {
        return get<T>(key,fallback);
      }
  };

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

}//namespace
#endif

