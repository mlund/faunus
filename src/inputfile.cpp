#include "faunus/inputfile.h"

namespace Faunus {
  inputfile::inputfile() { }

  inputfile::inputfile(string filename) {
    load(filename);
  }

  bool inputfile::load(string filename) {
    matrix.clear();
    calls.clear();
    string s;
    char cstr[256];
    file=filename;
    std::ifstream f( filename.c_str() );
    if (f) {
      while (!f.eof()) {
        dataformat tmp;
        f.getline(cstr,256);
        std::istringstream i( cstr );
        i >> tmp.name;
        if (tmp.name.find("#")==string::npos &&
            tmp.name.find("[")==string::npos && tmp.name.empty()==false )
        {
          while (i >> s)
            tmp.val.push_back(s);
          matrix.push_back(tmp);
          if (tmp.val.size()==0)
            std::cerr << "*** FATAL ERROR in Inputfile: '" << tmp.name << "' is defined but has no value ***" << endl;
        }
      }
      f.close();
      std::cout << "# Input parameters read from: " << filename << endl;
      return true;
    }
    else {
      std::cerr << "*** Failed to open inputfile ***" << endl;
      return false;
    }
  }

  int inputfile::findKey(string &key) {
    record_call(key);
    for (int i=0; i<matrix.size(); i++)
      if (matrix[i].name.compare(key)==0) return i;
    return -1;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  string inputfile::getstr(string key, string def) {
    int i = findKey(key);
    return (i!=-1) ? matrix[i].val[0] : def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  double inputfile::getflt(string key, double def) {
    int i = findKey(key);
    return (i!=-1) ? atof(matrix[i].val[0].c_str()) : def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  int inputfile::getint(string key, int def) {
    int i = findKey(key);
    return (i!=-1) ? atoi(matrix[i].val[0].c_str()) : def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  bool inputfile::getboo(string key, bool def) {
    int i = findKey(key);
    if (i!=-1) {
      if (matrix[i].val[0].compare("yes")==0)
        return true;
      else
        return false;
    }
    return def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  vector<string> inputfile::getvec(string key, string def) {
    int i = findKey(key);
    return (i!=-1) ? matrix[i].val : vector<string>(1,def);
  }

  //! \param key Name of the new keyword
  //! \param val String value
  void inputfile::add(string key, string val) {
    dataformat tmp;
    tmp.name=key;
    tmp.val.push_back(val);
    matrix.push_back(tmp);
  }

  //! \param key Name of the new keyword
  //! \param val Floating point number
  void inputfile::add(string key, double val) {
    std::ostringstream o;
    o << val;
    add(key, o.str());
  }

  void inputfile::record_call(string key) {
    bool newcall=true;
    for (int i=0; i<calls.size(); i++)
      if (calls[i].compare(key)==0) {
        newcall=false;
        break;
      }
    if (newcall==true)
      calls.push_back(key);
  }

  string inputfile::info() {
    std::ostringstream o;
    o << "\n"
      << "# INPUT FILE INFORMATION:\n"
      << "#   Config file = " << file << "\n"
      << "#   Accessed parameters:\n";
    for (int i=1; i<=calls.size(); i++) {
      if ( i%2!=0) o << "#   ";
      o << std::setw(19) << std::left << calls[i-1]
        << std::setw(19) << getstr(calls[i-1], "n/a");
      if ( i%2==0 ) o << endl;
    }
    o << endl;
    return o.str();
  }

  void inputfile::updateval(string key, string v) {
    int i=findKey( key);
    if (i<0) {
      std::cerr <<"# Could not find "<<key<<" during update!!!"<<endl;
    } else {
      matrix[i].val[0]=v;
    }
  }

  void inputfile::updateval(string key, double x) {
    std::ostringstream o;
    o << x;
    updateval(key,o.str());
  }

  string inputfile::print() {
    std::ostringstream o;
    //o <<"\n";
    for (int i=0; i<matrix.size(); i++) {
      o <<matrix[i].name;
      for (int j=0; j<matrix[i].val.size(); j++)
        o <<"   "<<matrix[i].val[j];
      o <<"\n";
    }
    return o.str();
  }

  //
  // CHECKVALUE CLASS
  //
  //
  unittest::unittest(inputfile &in) {
    stable = in.getboo("testsuite_stable", true);
    file = in.getstr("testsuite_testfile", "test.stable");
    if (stable==false)
      load(file);
    else {
      std::ostringstream o;
      o << "# Stable output generated on " << __DATE__ << " " << __TIME__
#ifdef __VERSION__
        << " using " << __VERSION__
#endif
#ifdef __SVN_REV__
        << " (SVN revision: " << __SVN_REV__ << ")"
#endif
        ;
      add( o.str(), "..." );
    }
  }

  /*!
   * \param name Name of test entry - spaces not allowed.
   * \param val Value to check (unstable) or save (stable)
   * \param threshold Maximum allowed absolute relative difference between stable and unstable
   */
  bool unittest::check(string name, double val, double threshold) {
    bool rc=true;
    // Stable: Save value.
    if (stable==true) {
      add(name, val);
      std::ofstream f( file.c_str() );
      if (f) {
        for (int i=0; i<matrix.size(); i++)
          f << matrix[i].name << " " << matrix[i].val[0] << std::endl;
        f.close();
      }
    }
    else {
      // Test value
      double ref=getflt(name,1e9),
             reldiff = std::abs( (ref-val)/ref );
      if (reldiff>threshold) {
        std::cerr.unsetf( std::ios_base::floatfield );
        std::cerr << "!!! Test " << std::setw(10) << name << " failed (ref,new,err): "
          << std::setprecision(4) 
             << std::setw(8)
             << ref << " "
             << std::setw(8) << val << " "
             << std::setw(8) << reldiff << " !!!" << std::endl;
        rc=false;
      }
    }
    result.push_back(rc); // Save outcome
    return rc;
  }

  bool unittest::smallerThan(string name, double x, double ref) {
    bool rc=false;
    if (x<ref)
      rc=true;
    else
      std::cerr << "!!! " << name << " smaller than test failed (" << x << " !< " << ref << ")" << endl;
    result.push_back(rc);
    return rc;
  }

  int unittest::returnCode() {
    for (int i=0; i<result.size(); i++)
      if (result[i]==false) return 1;
    return 0;
  }

  string unittest::report() {
    std::ostringstream o;
    unsigned int numerr = 0;
#ifdef __SUNPRO_CC
    for (int i=0; i<result.size(); i++)
      if (result[i]==false) numerr++;
#else
    for (int i=0; i<result.size(); i++)
      if (result[i]==false) numerr++;
    //numerr = std::count( result.begin(), result.end(), false);  
#endif
    o << std::endl << "# TEST SUITE:" << std::endl;
    if (stable==true)
      o << "#   Generated reference file: " << file << " with " << matrix.size()-1 << " item(s)." << std::endl;
    else
      o << "#   Test performed on " << result.size() << " item(s) with "
        << numerr << " errors." << std::endl;
    return o.str();
  }
}//namespace
