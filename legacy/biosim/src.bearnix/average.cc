#include "average.h"

average::average(int max) {
  maxcnt=max;
  sum=0; cnt=0;
  v.clear();
};

double average::operator+(average &a) { return avg() + a.avg(); };
double average::operator-(average &a) { return avg() - a.avg(); };
void average::operator+=(double x) { add(x); };

void average::add(double x) {
  sum+=x;
  cnt++;
  if (cnt==maxcnt) {
    v.push_back( sum/cnt );
    sum=0; cnt=0;
  };
};

double average::avg() {
  if (v.size()==0) return sum/cnt;
  return accumulate(v.begin(),v.end(),0.)/v.size();
};

double average::stddev() {
  double vav=avg(), sum=0;
  for (int i=0; i<v.size(); i++)
    sum+=pow(v[i]-vav, 2);
  return sqrt( sum / ( v.size() * (v.size()-1) ) );  
};

//! Print average of all sets and more than one set, also the std. deviation.
ostream &operator<<( ostream &output, average &a) {
  output << a.avg();
  if (a.v.size()>1) output << " " << a.stddev();
  return output;
};

void average::info() {
  cout << "# " << v.size() << " sets each containing " << maxcnt << " points." << endl;
};
