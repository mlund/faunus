#include "montecarlo.h"

trial::trial() {
  accepted=0;
  rejected=0;
  cnt=0;
};

float trial::percent(long long int total) {
  long long int n=cnt;
  if (total!=-1)
    n=total;
  return static_cast<float>(accepted)/n*100.;
};

void trial::operator+=(bool b) { add(b); }

void trial::add(bool b) {
  if (b==true)
    accepted++;
  if (b==false)
    rejected++;
  cnt++;
};

montecarlo::montecarlo(int macro, int micro) {
  //cout << "# Number of MC steps = " << micro*macro << endl;
  macroSteps = macro; //outer loop
  microSteps = micro; //inner loop
  accepted=0;
  rejected=0;
  energy_rejected=0;
  hc_rejected=0;
  total=0;
  time_i=time(0);

  t.resize(ENDENUM);

};

void montecarlo::accept(trialtype type) {
  t[type]+=true;
  accepted++;
};
void montecarlo::reject(trialtype type, rejectcause cause) {
  t[type]+=false;
  if (cause==ENERGY)
    energy_rejected++;
  if (cause==HC)
    hc_rejected++;
  rejected++;
};

void montecarlo::showStatus(int macroCnt) {
  progress =
  static_cast<double>(macroCnt)/macroSteps;

  elapsed_time = time(0)-time_i;

  remaining = int(static_cast<double>(elapsed_time)/
                  macroCnt*(macroSteps-macroCnt));

  rawtime = time(NULL) + remaining;
  timeinfo = localtime(&rawtime);

  string finishStr = asctime(timeinfo);

  cout << "# " << progress*100 << "%. Complete at: "
    << finishStr;
};

//Attempts to tune displacement parameter so that the acceptance
//ration is in the range [min;max]. Default is 30-40% (see header file)
void montecarlo::adjust_dp(trialtype type, double &dp, double min, double max) {
  double d=0,a;
  a=t[type].percent();
  if (type==ION) d=2.;
  if (type==TRANSLATE) d=0.5;
  if (type==ROTATE) d=0.1;
  if (type==MONOMER) d=0.4;
  if (type==CLUSTER) d=1.;
  
  if (a>max) dp += d; // too much accepted, increase dp.
  if (a<min) dp -= d; // too little accepted, decrease dp.
  if (dp<0) dp=-dp;   // we do not allow negative dp values.
  //if (dp>200.) dp=200.;
  if (type==ROTATE && dp>3.14) //no point going any further...
    dp=3.;
};

void montecarlo::showStatistics() {
  total=accepted+rejected;
  //for (int i=0; i<t.size(); i++)
  //  total+=t[i].cnt;

  total = accepted+rejected;
  cout << "MONTE-CARLO STATISTICS:" << endl
       << "-----------------------" << endl
       << "Elapsed time (h)  " << elapsed_time/3600. << endl
       << "Total steps       " << microSteps*macroSteps << endl
       << "  Macro           " << macroSteps << endl
       << "  Micro           " << microSteps
       << endl
       << "Total configs.    " << total << endl
       << "Accepted          " << accepted << " "
       << static_cast<double>(accepted)/total*100 <<"%\n"
       << "Rejected          " << rejected << " "
       << static_cast<double>(rejected)/total*100 <<"%\n"
       << "  Energy          " << energy_rejected << " "
       << static_cast<double>(energy_rejected)/total*100 <<"%\n"
       << "  Hardcore        " << hc_rejected << " "
       << static_cast<double>(hc_rejected)/total*100 << "%\n\n";

  cout << "Ion move       " << t[ION].percent() << endl
       << "Rotatation     " << t[ROTATE].percent() << endl
       << "Translate      " << t[TRANSLATE].percent() << endl
       << "Monomer move   " << t[MONOMER].percent() << endl
       << "Titration      " << t[TITRATE].percent() << endl
       << "Cluster move   " << t[CLUSTER].percent() << endl
       << "Site move      " << t[SITE].percent() << endl
       << endl;
};
