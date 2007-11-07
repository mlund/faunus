#include "inputfile.h"

inputfile::inputfile(string filename) {
  matrix.resize(50);
  int i=0;
  ifstream f( filename.c_str() );
  if (f) {
    while (!f.eof()) {
      f >> matrix[i].name >> matrix[i].val;
      i++;
    };
    matrix.resize(i-1);
    f.close();
  } else cout << "*** Failed to open inputfile ***" << endl;
};

//find matrix index of key
int inputfile::findKey(string &key) {
   for (int i=0; i<matrix.size(); i++)
      if (matrix[i].name.compare(key)==0) return i;
   return -1;
};

string inputfile::getStr(string key) {
  int i = findKey(key);
  if (i!=-1)
    return matrix[i].val;
  else
    cout << "Keyword '" << key << "' NOT found!\n";
  return "";
};

double inputfile::getDbl(string key, double def) {
  int i = findKey(key);
  if (i!=-1) return atof(matrix[i].val.c_str());
  else {
    cout << "Keyword '" << key << "' not found! - using "<<def<<"\n";
    return def;
  };
};

int inputfile::getInt(string key, int def) {
  int i = findKey(key);
  if (i!=-1) return atoi(matrix[i].val.c_str());
  else {
    cout << "Keyword '" << key << "' not found! -using "<<def<<"\n";
    return def;
  };
};

bool inputfile::getBool(string key, bool def) {
  int i = findKey(key);
  if (i!=-1) {
    if (matrix[i].val.compare("yes")==0)
      return true;
    else
      return false;
  } else {
    cout << "Keyword '" << key << "' not found! - using "<<def<<"\n";
    return def;
  };
};

