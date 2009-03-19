/*! \file Small test of uniform deviates
 *  
 *
 */
#include "faunus/faunus.h"

using namespace Faunus;
using namespace std;

int main() {
  inputfile in("random.conf");
  int macro=in.getint("macrosteps"), micro=in.getint("microsteps");
  average<double> one, half, onetimeshalf;
  average<double> onesq, halfsq;
  double O=0,H=0;
  histogram oneh(0.01, 0, 1);
  histogram halfh(0.01, -0.5, 0.5);
  vector<double> trace;
  trace.clear();
  io fio;
  slump slump;

  slump.random_seed(in.getint("seed", -7));

  cout<<in.info();

  cout << "# First 'slumpclass'::random_one() = "<<slump.random_one()<<endl<<endl;

  for (int mac=0; mac<macro; mac++) {
    for (int mic=0; mic<micro; mic++) {
      O=slump.random_one();
      H=slump.random_half();
      one+=O;
      half+=H;
      onetimeshalf+=O*H;
      onesq+=O*O;
      halfsq+=H*H;
      oneh.add(O);
      halfh.add(H);
      if (trace.size()<10000)
        trace.push_back(O);
    }
    cout <<"#Macrostep "<<mac+1<<endl;
  }
  cout <<endl
       <<"# Finished: macro*micro = "<<macro*micro<<" steps"<<endl
       <<"# Average 'slumpclass'::random_one()               = "<< one.avg()<<endl
       <<"# Average 'slumpclass'::random_half()              = "<< half.avg()<<endl
       <<"# Average 'slumpclass'::random_one()*random_half() = "<< onetimeshalf.avg()<<endl
       <<"# Average 'slumpclass'::random_one()^2             = "<< onesq.avg()<<endl
       <<"# Average 'slumpclass'::random_half()^2            = "<< halfsq.avg()<<endl;

  oneh.write("onedist.dat");
  halfh.write("halfdist.dat"); 

  string spur;
  int i=trace.size();
  std::ostringstream o;
  o << "# Trace of 'slumpclass'::random_one(), first ten thousand numbers" <<endl;
  for (int j=0; j<10000; j++)
    o <<j+1<<"    "<<trace[j]<<endl;
  fio.writefile("trace.dat", o.str());
}
