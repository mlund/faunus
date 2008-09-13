#include "faunus/inputfile.h"

namespace Faunus {
  //! \param filename Input file to scan
  inputfile::inputfile(string filename) {
    matrix.resize(50);
    int i=0;
    std::ifstream f( filename.c_str() );
    if (f) {
      while (!f.eof()) {
        f >> matrix[i].name >> matrix[i].val;
        i++;
      };
      matrix.resize(i-1);
      f.close();
      std::cout << "# Configuration read from: " << filename << endl;
    } else std::cout << "*** Failed to open inputfile ***" << endl;
  }

  int inputfile::findKey(string &key) {
    for (int i=0; i<matrix.size(); i++)
      if (matrix[i].name.compare(key)==0) return i;
    return -1;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  string inputfile::getstr(string key, string def) {
    int i = findKey(key);
    if (i!=-1)
      return matrix[i].val;
    return def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  double inputfile::getflt(string key, double def) {
    int i = findKey(key);
    if (i!=-1) return atof(matrix[i].val.c_str());
    return def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  int inputfile::getint(string key, int def) {
    int i = findKey(key);
    if (i!=-1) return atoi(matrix[i].val.c_str());
    return def;
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
}
