
#include "pot_ewald.h"

Ewald::Ewald(int size, double bjerrum, int ink) {
  lB = bjerrum; // bjerrum length in aangstroms
  alpha=0.001;
  alphasqrt=sqrt(alpha);
  double pi=acos(-1.);
  twopi = 2*pi;
  pisqrt = sqrt(pi);
  kmax=ink;               // default should be kmax=5 and ksqmax=27
  ksqmax=ink*ink+1;
  eix.resize(size);
  eiy.resize(size);
  eiz.resize(size);
  eixold.resize(size);
  eiyold.resize(size);
  eizold.resize(size);
  for(int i=0; i<size; i++) {
    eix[i].resize(kmax+1);
    eiy[i].resize(2*kmax+1);
    eiz[i].resize(2*kmax+1);
    eixold[i].resize(kmax+1);
    eiyold[i].resize(2*kmax+1);
    eizold[i].resize(2*kmax+1);
  };
};

Ewald::Ewald(int size, double bjerrum, double a, int ink) {
  lB = bjerrum; // bjerrum length in aangstroms
  alpha=a;
  alphasqrt=sqrt(a);
  double pi=acos(-1.);
  twopi = 2*pi;
  pisqrt = sqrt(pi);
  kmax=ink;
  ksqmax=ink*ink+1;
  eix.resize(size);
  eiy.resize(size);
  eiz.resize(size);
  eixold.resize(size);
  eiyold.resize(size);
  eizold.resize(size);
  for(int i=0; i<size; i++) {
    eix[i].resize(kmax+1);
    eiy[i].resize(2*kmax+1);
    eiz[i].resize(2*kmax+1);
    eixold[i].resize(kmax+1);
    eiyold[i].resize(2*kmax+1);
    eizold[i].resize(2*kmax+1);
  };
};

/*********************
  SYSTEM ENERGY
 *********************/

double Ewald::realSpaceEwald(vector<Particle> &p, Simbox &s) {
  double u=0.;
  unsigned int n = p.size();
  for (unsigned int i=0; i<n-1; i++)
    for (unsigned int j=i+1; j<n; j++)
      if (p[i].charge * p[j].charge != 0)
        u += realSpaceEwald(p[i], p[j],s);
  return lB*u;
};

double Ewald::realSpaceEwald(vector<Particle> &p, int j, Simbox &s) {
  double u=0.;
  unsigned int n=p.size();
  for(unsigned int i=0; i<j ; i++)
    u += realSpaceEwald(p[i],p[j],s);
  for(unsigned int i=j+1; i<n; i++)
    u += realSpaceEwald(p[i],p[j],s);
  return lB*u;
};

void Ewald::calcAlphaEwald(int size, Simbox &s, double inewaldpression) {
  ewaldpres=inewaldpression;
  int i = 0;
  double incr= 1e-3;
  double A=10.;
  double rcut2;

  rcut2=s.len_half*s.len_half;

  while( A > (1.+1e-10) || A < (1-1e-10)) {
    A=sqrt(alpha) * ewaldpres / exp(-alpha * rcut2);

    if(A>1){
      while(A > 1.){
        alpha-= incr;
        A = sqrt(alpha) * ewaldpres / exp(-alpha*rcut2);
        i++;
        if( i > 10000){
          cout << "ERROR ewaldpres = " << 1/sqrt(alpha)*exp(-alpha*rcut2) << endl;
          exit(1);
          return;
        };
      };
      if( i > 10000){
        cout << "ERROR ewaldpres = " << 1/sqrt(alpha)*exp(-alpha*rcut2) << endl;
        exit(1);
        return;
      };

      if( A < (1.-1e-10)) {
        alpha += 2*incr;
        A = sqrt(alpha)*ewaldpres / exp(-alpha*rcut2);
        incr = incr*0.1;
      };  
    } else {
      while( A < 1.){
        alpha += incr;
        A = sqrt(alpha)*ewaldpres / exp(-alpha*rcut2);
        if( i > 10000){
          cout << "ERROR ewaldpres = " << 1/sqrt(alpha)*exp(-alpha*rcut2) << endl;
          exit(1);
          return;
        };
      };
      if( i > 10000){
        cout << "ERROR ewaldpres = " << 1/sqrt(alpha)*exp(-alpha*rcut2) << endl;
        exit(1);
        return;
      };

      if( A < (1.-1e-10)) {
        alpha += 2*incr;
        A = sqrt(alpha)*ewaldpres / exp(-alpha*rcut2);
        incr = incr*0.1;
      };

    }; 
  };

  if(alpha < 0. ) {
    cout << "#### Alpha negative !!!! ###### " << endl;
    exit(1); 
  };
  cout << "# New alpha fit = " << alpha << endl;
  alphasqrt=sqrt(alpha);

  kmax = 1. + alpha*s.len_half*s.len/twopi*2.;
  cout << "# New kmax fit = " << kmax << endl;
  ksqmax = (1.+alpha*s.len_half*s.len/twopi*2.)*(1.+alpha*s.len_half*s.len/twopi*2.);

  for(int i=0; i<size; i++) {
    eix[i].resize(kmax+1);
    eiy[i].resize(2*kmax+1);
    eiz[i].resize(2*kmax+1);
    eixold[i].resize(kmax+1);
    eiyold[i].resize(2*kmax+1);
    eizold[i].resize(2*kmax+1);
  };

};




double Ewald::selfEwald(vector<Particle> &p){
  double u=0.;
  int size=p.size();
  for (int i=0; i<size; i++)
    u-= p[i].charge*p[i].charge*alphasqrt/pisqrt;
  return u*lB;
};

void Ewald::initKSpaceEwald(Simbox &s){
  double b;
  int maxk=2000;
  int ksq;
  double rkx,rky,rkz;
  double rksq;
  b= 1./4./alpha;
  double twopii=twopi/s.len;
  totk=0;

  for(int kx=0;kx<(kmax+1);kx++) {      // uses symmetry
    rkx = twopii*double(kx);
    for(int ky=-kmax;ky<(kmax+1);ky++) {
      rky = twopii*double(ky);
      for(int kz=-kmax;kz<(kmax+1);kz++) {
        rkz = twopii*double(kz);
        ksq = kx*kx+ky*ky+kz*kz;
        if( ksq < ksqmax && ksq!=0){
          totk+=1;
          if(totk > maxk) {
            cout << "K vector to small" << endl;
            exit(1);
          };

          rksq = rkx*rkx + rky*rky + rkz*rkz;
          kvec.push_back(twopi * exp( -b*rksq)/rksq/s.vol_A);
        };
      };
    };
  }; 
  kvec.resize(totk);
  eikr.resize(totk);
  eikrold.resize(totk);
  cout << "# Number of wavefunctions :" << totk << endl;

};

void Ewald::kSpaceEwald(vector<Particle> &p, Simbox &s) {
  int size=p.size();
  double twopii=twopi/s.len;

  for(int i=0; i<size; i++){
    eix[i][0]=complex<double>(1.0,0.0);
    eiy[i][kmax]=complex<double>(1.0,0.0);      //kmax = defines the '0'
    eiz[i][kmax]=complex<double>(1.0,0.0);

    eix[i][1]= complex<double>(cos(twopii*p[i].x),sin(twopii*p[i].x));
    eiy[i][kmax+1]= complex<double>(cos(twopii*p[i].y),sin(twopii*p[i].y));
    eiz[i][kmax+1]= complex<double>(cos(twopii*p[i].z),sin(twopii*p[i].z));

    eiy[i][kmax-1]= conj(eiy[i][kmax+1]);
    eiz[i][kmax-1]= conj(eiz[i][kmax+1]);
  };

  for(int kk=2;kk<(kmax+1);kk++) {
    for(int i=0; i<size; i++) {
      eix[i][kk]=eix[i][kk-1]*eix[i][1];
      eiy[i][kmax+kk]=eiy[i][kmax+kk-1]*eiy[i][kmax+1];
      eiy[i][kmax-kk]=conj(eiy[i][kmax+kk]);
      eiz[i][kmax+kk]=eiz[i][kmax+kk-1]*eiz[i][kmax+1];
      eiz[i][kmax-kk]=conj(eiz[i][kmax+kk]);
    };
  };

};

void Ewald::kSpaceEwald(vector<Particle> &p, Simbox &s, int j) {
  double twopii=twopi/s.len;

  eix[j][0]=complex<double>(1.0,0.0);
  eiy[j][kmax]=complex<double>(1.0,0.0);
  eiz[j][kmax]=complex<double>(1.0,0.0);

  eix[j][1]= complex<double>(cos(twopii*p[j].x),sin(twopii*p[j].x));
  eiy[j][kmax+1]= complex<double>(cos(twopii*p[j].y),sin(twopii*p[j].y));
  eiz[j][kmax+1]= complex<double>(cos(twopii*p[j].z),sin(twopii*p[j].z));

  eiy[j][kmax-1]= conj(eiy[j][kmax+1]); //1+kmax = -1   2+kmax = -2 osv.
  eiz[j][kmax-1]= conj(eiz[j][kmax+1]);

  for(int kk=2;kk<(kmax+1);kk++) {
    eix[j][kk]=eix[j][kk-1]*eix[j][1];
    eiy[j][kmax+kk]=eiy[j][kmax+kk-1]*eiy[j][kmax+1];
    eiy[j][kmax-kk]=conj(eiy[j][kmax+kk]);
    eiz[j][kmax+kk]=eiz[j][kmax+kk-1]*eiz[j][kmax+1];
    eiz[j][kmax-kk]=conj(eiz[j][kmax+kk]);
  };

};

double Ewald::sumkSpaceEwald(vector<Particle> &p) {
  int size=p.size();
  complex<double> sum;
  int ksq;
  double u=0.;
  double fact;
  totk=0;

  for(int kx=0; kx<(kmax+1);kx++){
    if(kx==0) {
      fact=1.0;
    } else {
      fact=2.0;
    };
    for(int ky=0; ky<(2*kmax+1); ky++) {
      for(int kz=0; kz<(2*kmax+1); kz++) {
        ksq=kx*kx+(ky-kmax)*(ky-kmax)+(kz-kmax)*(kz-kmax);
        if(ksq < ksqmax && ksq!=0){
          sum = complex<double>(0.0,0.0);
          for(int i=0;i<size; i++)
            sum += p[i].charge*eix[i][kx]*eiy[i][ky]*eiz[i][kz];
          eikr[totk] = sum;
          u += fact*kvec[totk]*real(sum*conj(sum));
          totk+=1;
        };
      };
    };
  };
  return u*lB;
};

double Ewald::sumkSpaceEwald(vector<Particle> &p, int j) {
  complex<double> sum;
  int ksq;
  double u=0.;
  double fact;
  totk=0;

  for(int kx=0; kx<(kmax+1);kx++){
    if(kx==0) {
      fact=1.0;
    } else {
      fact=2.0;
    };
    for(int ky=0; ky<(2*kmax+1); ky++) {
      for(int kz=0; kz<(2*kmax+1); kz++) {
        ksq=kx*kx+(ky-kmax)*(ky-kmax)+(kz-kmax)*(kz-kmax);
        if(ksq < ksqmax && ksq!=0){
          eikr[totk] = eikrold[totk]+p[j].charge*(eix[j][kx]*eiy[j][ky]*eiz[j][kz]
              -eixold[j][kx]*eiyold[j][ky]*eizold[j][kz]);
          sum = eikr[totk];

          u += fact*kvec[totk]*real(sum*conj(sum));
          totk+=1;
        };
      };
    };
  };
  return u*lB;
};
