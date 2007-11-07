// Class for specifying a data potential
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "point.h"

class datapot {
  private:
    double res;
    vector<double> a[2][2];
  public:
    datapot(double);
    bool loaddata(string);
    void printdata(int, int);
    double energy(particle &p1, particle &p2) {
    };
};

datapot::datapot(double resolution=0.2) {
  res=resolution; //data resolution in angstrom
};

bool datapot::loaddata(string filename) {
  //a[1][2].push_back(2.);
  string s;
  vector<string> data;
  ifstream f(filename.c_str());
  if (f) {
    //read entire file
    while (!f.eof()) {
      getline(f,s);
      data.push_back(s); 
    };
    //find data sets, starting with "#"
    for (int i=0; i<data.size(); i++) {
      if (data[i].find("#")==0) //look for new set
        cout << data[i] << endl;
    };
    return true;
  }
  else
    cout << "*** Failed to open empirical potential file ***\n";
  return false;
};
