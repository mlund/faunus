#include "inputfile.h"

inputfile::inputfile(string filename) {
  matrix.resize(80);
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

void inputfile::updateInp(double nPos, double nNeg, string inp) {
  int s=matrix.size();
  int j=0;
  ofstream f( inp.c_str() );
  if (f) {
    f << setfill(' ') << setiosflags(ios::left) << endl;
    for (int i=0; i<s; i++) {
      if (matrix[i].name == "nion2" || matrix[i].name == "nion1") {
        if (matrix[i].name == "nion1")
          f << setw(13) <<matrix[i].name  << nPos <<endl;
        if (matrix[i].name == "nion2")
          f << setw(13) <<matrix[i].name  << nNeg <<endl;
        }else{
        f << setw(13) <<matrix[i].name  << matrix[i].val <<endl;
      };
    };
    f.close();
  cout << "# --- Input updated ---"<<endl;
  } else cout << "*** Failed to update input ***"<<endl;
};

void inputfile::setcatPot( string& fileName ) {  //Helpfunction to update the input file doing GC
  ifstream f( fileName.c_str());
  if (f) {
    vector<string> catPotval(10);
    int i=0;
    while (!f.eof()) {
      f >> catPotval[i];
      ++i;
    };
    f.close();
    int s=matrix.size();
    int j=0;
    while (j<s){
      if (matrix[j].name == "catPot")
        matrix[j].val = catPotval[0];
      if (matrix[j].name == "nion1")
        matrix[j].val =catPotval[1];
      if (matrix[j].name == "nion2")
        matrix[j].val =catPotval[2];
      j++;
    };
    cout << "# Chemical potential updated from previus run!" <<endl
         << "# You are now entering the Semi Canonical domain!!!"<<endl;
  } else {
    cout << "*** Failed to update chemical potential***"<<endl;
  };
};
/*  ofstream f( fileName.c_str() );
  if (f) {
    f << setfill(' ') << setiosflags(ios::left) << endl;
    for (int i=0; i<s; i++) {
      f << setw(13) <<matrix[i].name  << matrix[i].val <<endl;
    
    };
    f.close();
  } else cout << "*** Failed to update chemical potential ***"<<endl;*/


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

