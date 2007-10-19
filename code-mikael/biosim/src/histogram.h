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
  unsigned int bin, binmax;          //bin identifier and its max. value
  enum typeCode {AVERAGE, PROBABILITY, PENALTY};

  struct data {
    typeCode type;            //type of histogram
    string name;              //user-defined name
    unsigned long int cnt;             //total number of counts
    vector<unsigned long long int> hist;  //integer histogram
    vector<double> dhist;   //double histogram
  };

  double width; //width of each bin

  histogram(double,double);                //init: specify width and maxvalue
  void init(data &);                       //initializes data structure
  void init(vector<data> &);               // -//-
  double get(data &, double);       //get a value from histogram
  double getavg(data &, double);    //get average value
  void add(data &, double);         //adds a value to the histogram
  void add(data &, double,double);  //adds a value to the histogram
  void show(vector<data> &, double,double);//print histogram to screen
  bool save(data &, string);               //save histogram to file
  bool load(data &, string);               //load from file
};
#endif
