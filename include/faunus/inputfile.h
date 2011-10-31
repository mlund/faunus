#ifndef FAU_INPUTFILE_H
#define FAU_INPUTFILE_H

#include <faunus/common.h>

namespace Faunus {

  /*!
   * \brief Retrieve parameters from a formatted input file
   * \author Mikael Lund
   *
   * The input file is expected to have the following
   * format:
   * \li  KEYWORD VALUE (space separated)
   *
   * Blank lines or lines containing [ or # are ignored.
   *
   * Example:
   * \code
   * InputMap in("myinput.in");
   * double T = in.get<double>("temperature", 298.15);
   * \endcode
   */
  class InputMap {
    private:
      std::map<string,string> map;
      vector<string> usedkeys;
      vector<string> incfiles;
    protected:
      bool save(string);     //!< Save map to disk
    public:
      InputMap();
      InputMap(string);      //!< Construct and include input file.
      bool include(string);  //!< Include keyword/value file.
      string info();         //!< Information string about read files and keywords
      //!< Add a keyword and an associated value
      template<typename T> void add(const string &key, T value) {
        std::ostringstream o;
        o << value;
        map[key]=o.str();
      }

      //!< Get value associated with keyword
      template<typename T> T get(const string &key, T fallback) {
        if ( map.find(key)!=map.end() ) {
          std::istringstream i( map[key] );
          i >> fallback;
        }
        return fallback;
      }
  };

  /*!
   * \brief Class for checking generated output against stored data
   * \author Mikael Lund
   *
   * This class is used to test generated output by comparing against
   * a database file generated for a "stable" run that are known to work -
   * a reference system.
   * The input file passed to the constructor is expected to have the following
   * keywords:
   * format:
   * \li \c test_stable - Specifies if the system is "stable" or need testing.
   * \li \c test_file   - Name of test file to load (if unstable) or generate (if stable).
   *
   * To check a value, simply call unittest(name,value,precision%). If stable==true this will generate
   * a new reference testfile. If stable==false this will check the given value against the
   * loaded testfile.
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
      string info(); //!< Information about passed and failed tests
  };

}//namespace
#endif

