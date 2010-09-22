#include "faunus/slump.h"

namespace Faunus {
  // Baseclass
  double random::random_half() {
    return -0.5 + random_one();
  }

  bool random::runtest(float f) {
    return (random_one()<f) ? true : false;
  }

  std::string random::info() {
    std::ostringstream o;
    o << "# RANDOM NUMBER GENERATOR:" << std::endl
      << "#   Scheme: " << name << std::endl;
    return o.str();
  }

  /*
   *  Generic C++ generator
   */
  randomDefault::randomDefault() {
    name = "C++ build in";
    rand_max_inv = 1./(RAND_MAX);
  }

  double randomDefault::random_one() {
    double r=rand_max_inv*rand();
    return (r>=1) ? r-1e-5 : r;  // we don't like *exactly* 1 !!
  }

  void randomDefault::random_seed(int s) {
    (s!=0) ? srand(s) : srand(time(0));
  }

  unsigned int randomDefault::rand() {
    return rand();
  }

#ifdef HAVETR1
  /*
   *  Mersenne Twister generator (STL TR1 build-in)
   */
  randomTwister::randomTwister() : dist(0.0,1.0) {
    name="Mersenne Twister (C++ TR1)";
    maxinv=1./(eng.max()+1.);
  }

  double randomTwister::random_one() {
    double x;
    #pragma omp critical
    x=dist(eng);
    return x*maxinv;
  }

  void randomTwister::random_seed(int s) {
    eng.seed(s);
  }

  unsigned int randomTwister::rand() {
    return dist(eng);
  }

#endif

  // "Ran2" - ran2 from 'Numerical Recipies'
  const double ran2::EPS=3.0e-16;
  ran2::ran2() {
    name="Numerical recipes 'ran2' (not thread safe!)";
    AM=1.0/2147483563.0;
    RNMX=1.0-3.0e-16;
    iy=0;
    idum2=123456789;
  }

  double ran2::random_one()  {
    #pragma omp critical
    {
      random_seed(idum);
      k=idum/IQ1;
      idum=IA1*(idum-k*IQ1)-k*IR1;
      if (idum<0)
        idum+=IM1;
      k=idum2/IQ2;
      idum2=IA2*(idum2-k*IQ2)-k*IR2;
      if (idum2<0)
        idum2+=IM2;
      j=iy/NDIV;
      iy=iv[j]-idum2;
      iv[j]=idum;
      if (iy<1)
        iy+=IMM1;
    }
    double temp;
    if((temp=AM*iy)>RNMX)
      return RNMX;
    else
      return temp;
  }

  void ran2::random_seed(int s) {
    idum=s;    
    if(idum<=0) {
      idum=(idum==0 ? 1: -idum);
      idum2=idum;
      for (j=NTAB+7; j>=0; j--) {
        k=idum/IQ1;
        idum=IA1*(idum-k*IQ1)-k*IR1;
        if (idum<0) 
          idum+=IM1;
        if (j<NTAB)
          iv[j]=idum;
      }
      iy=iv[0];
    }
  }

  unsigned int ran2::rand() {
    return random_one()*RAND_MAX;
  }

  slump slp;
}//namespace

