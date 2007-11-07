#include "histogram.h"

histogram::histogram(double binwidth, double maxvalue) {
  width = binwidth;
  binmax = int(maxvalue/width+.5);
};

void histogram::init(histogram::data &d) {
  d.hist.reserve(binmax);
  d.hist.resize(binmax);
  d.dhist.reserve(binmax);
  d.dhist.resize(binmax);
  d.cnt=0;
  for (unsigned int i=0; i<d.hist.size(); i++) {
    d.hist[i] = 0;
    d.dhist[i] = 0;
  };
};
void histogram::init(vector<histogram::data> &dvec) {
  for (int i=0; i<dvec.size(); i++)
    histogram::init(dvec[i]);
};

//Shows a vector of histogram data (same size!).
void histogram::show(vector<histogram::data> &d,
		     double minvalue=0, double maxvalue=0) {
  int min=int(minvalue/width+.5); //first hist to show
  int max=int(maxvalue/width+.5); //last hist to show
  if (maxvalue==0 || max>=d[0].hist.size() )
    max=binmax;

  cout << "#DISTRIBUTION FUNCTIONS:\n"
       << "#1 distance" << endl;

  for (int j=0; j<d.size();j++)
    cout << "#" << j+2 << " " << d[j].name << endl;

  for (int i=min; i<max; i++) {
    cout << i*width << " ";
    for (int j=0; j<d.size(); j++) {
      if (d[j].type==AVERAGE)
        cout << d[j].dhist[i]/double(d[j].hist[i]) << " ";
      if (d[j].type==PROBABILITY)
        cout << d[j].hist[i]/double(d[j].cnt) << " ";
      if (d[j].type==PENALTY)
        cout << d[j].dhist[i] << " ";
    };
    cout << endl;
  };
};

bool histogram::save(histogram::data &d, string file) {
  ofstream f( file.c_str() );
  if (f) {
    f.precision(30);
    for (int i=0; i<d.hist.size(); i++)
      f << d.hist[i] << " " << d.dhist[i] << endl;
    f.close();
    cout << "# Histogram: " << d.name << " saved to disk.\n";
    return true;
  };
  return false;
};

bool histogram::load(histogram::data &d, string file) {
  ifstream f(file.c_str());
  if (f) {
    f.precision(30);
    for (int i=0; i<d.hist.size(); i++)
      f >> d.hist[i] >> d.dhist[i];
    f.close();
    cout << "# Histogram: " << d.name << " loaded from disk.\n";
    return true;
  };
  return false;
};

//add to probability
void histogram::add(histogram::data &d, double val) {
  d.hist[int(val/width+.5)]++;
  d.cnt++;
};

//add to average
void histogram::add(histogram::data &d, double val, double x) {
  int i=int(val/width+.5);
  d.hist[i]++;
  d.dhist[i] += x;
  d.cnt++;
};

double histogram::get(histogram::data &d, double val) {
  return d.dhist[int(val/width+.5)];
};

double histogram::getavg(histogram::data &d, double val) {
  return d.dhist[int(val/width+.5)] / double( d.hist[int(val/width+.5)] );
};
