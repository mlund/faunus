#include "faunus/inputfile.h"

namespace Faunus {
  inputfile::inputfile(string filename) {
    file=filename;
    dataformat tmp;
    std::ifstream f( filename.c_str() );
    if (f) {
      while (!f.eof()) {
        f >> tmp.name;
        if (tmp.name.find("#")!=string::npos ||
            tmp.name.find("[")!=string::npos)
          f.ignore(256, '\n');
        else {
          f >> tmp.val;
          matrix.push_back(tmp);
        }
      };
      f.close();
      std::cout << "# Input parameters read from: " << filename << endl;
    } else {
      std::cout << "*** Failed to open inputfile ***" << endl;
      throw;
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
    return (i!=-1) ? matrix[i].val : def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  double inputfile::getflt(string key, double def) {
    int i = findKey(key);
    return (i!=-1) ? atof(matrix[i].val.c_str()) : def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  int inputfile::getint(string key, int def) {
    int i = findKey(key);
    return (i!=-1) ? atoi(matrix[i].val.c_str()) : def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  bool inputfile::getboo(string key, bool def) {
    int i = findKey(key);
    if (i!=-1) {
      if (matrix[i].val.compare("yes")==0)
        return true;
      else
        return false;
    }
    return def;
  }

  //! \param key Name of the new keyword
  //! \param val String value
  void inputfile::add(string key, string val) {
    dataformat tmp;
    tmp.name=key;
    tmp.val=val;
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
      << "#  INPUT FILE INFORMATION:\n"
      << "#    Config file = " << file << "\n"
      << "#    Accessed parameters:\n";
    for (int i=1; i<=calls.size(); i++) {
      if ( i%2!=0) o << "#      ";
      o << std::setw(15) << std::left << calls[i-1]
        << std::setw(20) << getstr(calls[i-1], "n/a");
      if ( i%2==0 ) o << endl;
    }
    o << endl;
    return o.str();
  }
}
