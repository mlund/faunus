#include "widommod2.h"

Widommod::Widommod(int in, double radi, double charg, string nm){
  name = nm;
  charge = charg;
  radius = radi;
  cnt=0;
  ghostin = in;
  chel=0.;
  chhc=0.;
  chex=0.;
  chtot=0.;
  ewden.resize(11);
  ewnom.resize(11);
  chint.resize(11);
  chexw=0.;
  ihc=0;
  irej=0;

  for(int j=0; j<11; j++) {
      ewden[j]=0.;
      ewnom[j]=0.;
      chint[j]=0.;
  };
  
};

void Widommod::ghostInsert(vector<particle> &p,  interact &ip , collision &coll, cell &c) {
  particle ghost;
  double u,cu;      //do they need to be initiated???
  
  cnt++;
  for(int i=0; i < ghostin; i++) {
    c.randomPos( ghost); //ghost.randomPos(s);
    int goverlap=0;
    ghost.radius = radius;
    irej=0;
    int j=0;
    while(!coll.overlap(p, ghost) && j<p.size()){
      j++;
    };
    if(j != p.size()) {     
      ihc++;
      irej=1;
      goverlap++;
    };
   
    
    if(goverlap!=1) {  //gspec, 0 or 1 ???
      cu=0.;
      u=0.;
      double invdi;
      for(int l=0; l < p.size(); l++) {
        invdi=ghost.invdist(p[l]);
        cu+=invdi;
        u+=invdi*p[l].charge;
      };

      
      cu=cu*ip.lB;
      u=u*ip.lB;
      double ew,ewla,ewd;

      if(irej==0) {
        expuw+=exp(-u*charge);
        for(int cint=0; cint<11; cint++) {
          ew=charge*(u-double(cint)*0.1*charge*cu/double(p.size()));
          ewla = ew*double(cint)*0.1;
          ewd=exp(-ewla);
          ewden[cint]+=ewd;
          ewnom[cint]+=ew*ewd;
        };
      };        
    };
  };
};

void Widommod::getWidomResults() {
  double aint4,aint2,aint1;
    for(int cint=0; cint<11; cint++) {
      if(ewden[cint]==0) {
        cout << "Widom denominator equals zero" << endl;
      } else {
        chint[cint]=ewnom[cint]/ewden[cint];
      };
    
    aint4=chint[1]+chint[3]+chint[5]+chint[7]+chint[9];
    aint2=chint[2]+chint[4]+chint[6]+chint[8];
    aint1=chint[0]+chint[10];  
    chel=1./30.*(aint1+2*aint2+4*aint4);
  };
  
  double cnttot;
  cnttot=cnt*ghostin;
  cout << "#" << endl;
  cout << "# SimpleScaling WIDOM: "<<name << endl;
  cout << "# Radius (A) " << radius <<endl
       << "# Charge (e) " << charge <<endl;
  cout << "# cnttot: " << cnttot << endl;
  chhc=-log(double(cnttot-ihc)/cnttot);
  chexw=-log(expuw);
  chex=chhc+chel;
  cout << "# Widom excess:  " << chex  
       << endl;
  cout << "# Widom electro: " << chel << endl;
  cout << "# Widom hc     : " << chhc << endl;
  cout << "#" << endl;
  
  
};

double Widommod::exessChempot() {
  double aint4,aint2,aint1;
    for(int cint=0; cint<11; cint++) {
      if(ewden[cint]==0) {
        cout << "Widom denominator equals zero" << endl;
      } else {
        chint[cint]=ewnom[cint]/ewden[cint];
      };

    aint4=chint[1]+chint[3]+chint[5]+chint[7]+chint[9];
    aint2=chint[2]+chint[4]+chint[6]+chint[8];
    aint1=chint[0]+chint[10];
    chel=1./30.*(aint1+2*aint2+4*aint4);
  };

  double cnttot;
  cnttot=cnt*ghostin;
  chhc=-log(double(cnttot-ihc)/cnttot);
  chexw=-log(expuw);
  chex=chhc+chel;
  return chex;
}

