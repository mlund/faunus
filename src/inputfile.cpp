#include "faunus/inputfile.h"

namespace Faunus {
  inputfile::inputfile(string filename) {
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

  int inputfile::findKey(string &key) const {
    for (int i=0; i<matrix.size(); i++)
      if (matrix[i].name.compare(key)==0) return i;
    return -1;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  string inputfile::getstr(string key, string def) const {
    int i = findKey(key);
    return (i!=-1) ? matrix[i].val : def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  double inputfile::getflt(string key, double def) const {
    int i = findKey(key);
    return (i!=-1) ? atof(matrix[i].val.c_str()) : def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  int inputfile::getint(string key, int def) const {
    int i = findKey(key);
    return (i!=-1) ? atoi(matrix[i].val.c_str()) : def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  bool inputfile::getboo(string key, bool def) const {
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
}
