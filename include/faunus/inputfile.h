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
   */
  class inputfile {
    protected:
      struct dataformat {
        string name;
        vector<string> val;
      };
      vector<dataformat> matrix;
      vector<string> calls;
      int findKey(string &);
      void record_call(string);
      string file;
    public:
      inputfile();
      inputfile(string);                     //!< Constructor
      bool load(string);                     //!< Load inputfile from disk
      string getstr(string, string="");      //!< Get string value
      double getflt(string, double=0);       //!< Get double value
      int getint(string, int=0);             //!< Get integer value
      bool getboo(string, bool=false);       //!< Get boolean value
      vector<string> getvec(string,string);  //!< Get vector of strings
      void add(string,string);               //!< Add an entry to the loaded list
      void add(string,double);               //!< Add an entry to the loaded list
      string info();                         //!< Show info
      string print();                        //!< Print a string of the inputfile
      void updateval(string, string);        //!< Update the inputfile for next run
      void updateval(string, double);        //!< Update the inputfile for next run
      bool write();                          //!< Write to input file
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
   * \li  testsuite_stable   - Specifies if the system is "stable" or need testing.
   * \li  testsuite_testfile - Name of test file to load (if unstable) or generate (if stable).
   *
   * To check a value, simply call Check(name,value,threshold). If the stable==true this will generate
   * a new reference testfile. If stable==false this will check the given value against the
   * loaded testfile.
   */
  class checkValue : private inputfile {
    private:
      vector<bool> result;
    public:
      bool stable;                       //!< True if test suite is stable (=reference)
      checkValue(inputfile &);           //!< Read parameters from inputfile
      bool check(string, double, double=1e-2); //!< Check or store value depending on the stable state.
      string report();                   //!< Print report.
      int returnCode();                  //!< Zero if no errors, one otherwise.
  };
}
#endif

