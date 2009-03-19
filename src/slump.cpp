#include "faunus/slump.h"

namespace Faunus {
  // Baseclass
  double random::random_half() { return -0.5 + random_one(); }
  bool random::runtest(float f) { return (random_one()<f) ? true : false; }

  // "slump" - the default generator
  randomDefault::randomDefault() { rand_max_inv = 1./RAND_MAX; }
  double randomDefault::random_one() { return rand_max_inv*rand(); }
  void randomDefault::random_seed(unsigned int s) { (s!=0) ? srand(s) : srand(time(0)); }

  // "Twister" - Mersenne Twister generator
  double randomTwister::random_one() { return mt.rand(); }
  void randomTwister::random_seed(unsigned int s) { mt.seed(s); }

  // "Ran2" - ran2 from 'Numerical Recipies'
  ran2::ran2() {
    AM=1.0/2147483563.0;
    RNMX=1.0-3.0e-16;
    iy=0;
    idum2=123456789;
  }

  double ran2::random_one()  {
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
    if((temp=AM*iy)>RNMX)
      return RNMX;
    else
      return temp;
  }
  void ran2::random_seed( unsigned int s ) {
    idum=s;    
    if(idum<=0) {
      std::cout << "insied seeder!"<<std::endl;
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
}

