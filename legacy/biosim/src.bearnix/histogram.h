#ifndef _histogram_h
#define _histogram_h

// HISTOGRAM CLASS - Manages histogram data
// M. Lund, 2002
#include <iostream>
#include <valarray>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>

using namespace std;

class histogram {
 public:
  int bin, binmax;          //bin identifier and its max. value
  enum typeCode {AVERAGE, PROBABILITY, PENALTY};

  struct data {
    typeCode type;            //type of histogram
    string name;              //user-defined name
    long int cnt;             //total number of counts
    vector<long long int> hist;  //integer histogram
    vector<double> dhist;   //double histogram
  };

  double width; //width of each bin

  histogram(double,double);                //init: specify width and maxvalue
  void init(data &);                       //initializes data structure
  void init(vector<data> &);               // -//-
  inline double get(data &, double);       //get a value from histogram
  inline double getavg(data &, double);    //get average value
  inline void add(data &, double);         //adds a value to the histogram
  inline void add(data &, int, int);
  inline void add(data &, double,double);  //adds a value to the histogram
  void show(vector<data> &, double,double);//print histogram to screen
  bool save(data &, string);               //save histogram to file
  bool load(data &, string);               //load from file
};

//add to probability
inline void histogram::add(histogram::data &d, double val) {
  d.hist[int(val/width+0.5)]++;
  d.cnt++;
};
inline void histogram::add(histogram::data &d, int val, int check) { //check is a dummy
  d.hist[val]++;
  d.cnt++;
};

//add to average
inline void histogram::add(histogram::data &d, double val, double x) {
  int i=int(val/width+0.5);
  d.hist[i]++;
  d.dhist[i] += x;
  d.cnt++;
};

inline double histogram::get(histogram::data &d, double val) {
  return d.dhist[int(val/width+0.5)];
};

inline double histogram::getavg(histogram::data &d, double val) {
  return d.dhist[int(val/width+0.5)] / double( d.hist[int(val/width+0.5)] );
};
#endif
