#ifndef _montecarlo_h
#define _montecarlo_h


/*  MONTE CARLO CLASS
    -------------------------------------------*/
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

class trial {
 private:
  void add(bool);
 public:
  long long int accepted,rejected,cnt;
  string name;
  trial();
  float percent(long long int=-1);
  void operator+=(bool);
};


class montecarlo {
private:
  double progress;
  long long int accepted, rejected, energy_rejected, hc_rejected, total;
  int elapsed_time, remaining;
  time_t rawtime;
  struct tm *timeinfo;

public:
  enum trialtype {ION=0, ROTATE, TRANSLATE, TITRATE, MONOMER, CLUSTER, SITE, INFLOW, OUTFLOW, ENDENUM};
  enum rejectcause {ENERGY=0, HC};
  vector<trial> t;

  int macroSteps, microSteps;
  int time_i;
  montecarlo(int, int);
  void showStatus(int);
  void showStatistics();
  void accept(trialtype);
  void reject(trialtype, rejectcause);
  void adjust_dp(trialtype, double &, double=30, double=40); //specify acceptance range
};
#endif
