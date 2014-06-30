#include <faunus/slump.h>

namespace Faunus {
  const double RandomRan2::EPS=3.0e-16;

  RandomRan2::RandomRan2() {
    name="Numerical Recipes 'ran2'";
    AM=1.0/2147483563.0;
    RNMX=1.0-3.0e-16;
    iy=0;
    idum2=123456789;
    seed(-13);
  }

  double RandomRan2::_randone()  {
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

  slump slp_global;
}//namespace

