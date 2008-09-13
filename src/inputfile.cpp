#include "faunus/inputfile.h"

namespace Faunus {
  //! \param filename Input file to scan
  inputfile::inputfile(std::string filename) {
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
      std::cout << "# Configuration read from: " << filename << std::endl;
    } else std::cout << "*** Failed to open inputfile ***" << std::endl;
  }

  int inputfile::findKey(std::string &key) {
    for (int i=0; i<matrix.size(); i++)
      if (matrix[i].name.compare(key)==0) return i;
    return -1;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  std::string inputfile::getstr(std::string key, std::string def) {
    int i = findKey(key);
    if (i!=-1)
      return matrix[i].val;
    return def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  double inputfile::getflt(std::string key, double def) {
    int i = findKey(key);
    if (i!=-1) return atof(matrix[i].val.c_str());
    return def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  int inputfile::getint(std::string key, int def) {
    int i = findKey(key);
    if (i!=-1) return atoi(matrix[i].val.c_str());
    return def;
  }

  //! \param key Keyword to look for
  //! \param def Default value if keyword is not found
  bool inputfile::getboo(std::string key, bool def) {
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
