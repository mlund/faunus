#include "faunus/slump.h"

namespace Faunus {
  // Baseclass
  double RandomBase::randHalf() {
    return -0.5 + randOne();
  }

  bool RandomBase::runtest(float f) {
    return ( randOne()<f) ? true : false;
  }

  std::string RandomBase::info() {
    std::ostringstream o;
    o << "# RANDOM NUMBER GENERATOR:" << std::endl
      << "#   Scheme: " << name << std::endl;
    return o.str();
  }

  /*
   *  Generic C++ generator
   */
  RandomDefault::RandomDefault() {
    name = "C++ build in";
    rand_max_inv = 1./(RAND_MAX);
  }

  double RandomDefault::randOne() {
    double r;
    #pragma omp critical
    r=rand_max_inv*rand();
    return (r>=1) ? r-1e-5 : r;  // we don't like *exactly* 1 !!
  }

  void RandomDefault::seed(int s) {
    (s!=0) ? srand(s) : srand(time(0));
  }

  unsigned int RandomDefault::rand() {
    double x;
    #pragma omp critical
    x=rand();
    return x;
  }

  /*
   *  Mersenne Twister generator (STL TR1 build-in)
   */
  RandomTwister::RandomTwister() : dist(0.0,1.0) {
    name="Mersenne Twister (C++ TR1)";
    maxinv=1./(eng.max()+1.);
  }

  double RandomTwister::randOne() {
    double x;
    #pragma omp critical
    x=dist(eng)*maxinv;
    return x;
  }

  void RandomTwister::seed(int s) {
    #pragma omp critical
    eng.seed(s);
  }

  unsigned int RandomTwister::rand() {
    return dist(eng);
  }

  // "Ran2" - ran2 from 'Numerical Recipies'
  const double RandomRan2::EPS=3.0e-16;
  RandomRan2::RandomRan2() {
    name="Numerical recipes 'ran2' (not thread safe!)";
    AM=1.0/2147483563.0;
    RNMX=1.0-3.0e-16;
    iy=0;
    idum2=123456789;
  }

  double RandomRan2::randOne()  {
    #pragma omp critical
    {
      seed(idum);
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

  void RandomRan2::seed(int s) {
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

  unsigned int RandomRan2::rand() {
    return randOne()*RAND_MAX;
  }

  slump slp_global;
}//namespace

