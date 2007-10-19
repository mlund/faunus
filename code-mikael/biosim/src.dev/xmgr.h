#include "point.h"

class xmgr {
  struct head {
    string title;
    string xlabel;
    string ylabel;
    string comment;
    float xmin,xmax,ymin,ymax;
  };
  xmgr();
  bool save(string, vector<point> &, vector<point> &);  
};

bool xmgr::save(string filename,
    vector<point> x, vector<point> y) {

  cout << "@TYPE xy\n"
    << "legend on\n"
    << "
};
