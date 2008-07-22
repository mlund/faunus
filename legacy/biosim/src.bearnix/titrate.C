#include "titrate.h"

/*
 * CONSTRUCTOR (sets pH, finds titrateable sites)
 *
 * p      = full particle vector
 * mobile = salt group to look for protons (+1 ions)
 * peeage = pH
 * 
 */
titrate::titrate(double peeage) {
  ph=peeage;
};
titrate::titrate(vector<particle> &p, group &mobile, double peeage) {
  ph=peeage;
  init(p, mobile);
};

void titrate::init(vector<particle> &p, group &mobile) {
  for (int i=0; i<p.size(); i++)    //search for titrateable sites
    if (p[i].id!=-1 && d[p[i].id].pka!=0) {
        sites.push_back(i);
        p[i].charge=d[p[i].id].charge;
        eqpos.push_back(p[i]);
      };
  for (int i=mobile.beg; i<mobile.end+1; i++) {
    if (p[i].charge==+1) protons.push_back(i);
    if (p[i].charge==0) neutrons.push_back(i);
  };
  q.resize( sites.size() ); //adjust site average charge vector
};

void titrate::update(vector<particle> &p, group &salt) {
  protons.erase(protons.begin(), protons.end());
  neutrons.erase(neutrons.begin(), neutrons.end());
  for (int i=salt.beg; i<salt.end+1; i++) {
    if (p[i].charge==+1) protons.push_back(i);
    if (p[i].charge==0 ) neutrons.push_back(i);
  };
};

//returns number a random, titrateable
//site in particle vector
short int titrate::random(vector<short int> &iv) {
  short int i = rand() % iv.size();     //and pick randomly
  return iv[i];
};

//Returns protonation status of particle i in particle vector
titrate::keywords titrate::status(vector<particle> &p, short int i) {
  if (p[i].charge>d[p[i].id].charge)
    return PROTONATED;
  return DEPROTONATED;
};

titrate::action titrate::takeFromBulk(vector<particle> &p, short int i, short int j) {
  if (j==-1)                               //if not specified..
    j=random(protons);                     //..pick a random proton
  
  protons.erase(remove(protons.begin(), protons.end(), j));

  p[i].charge+=1;                            //protonate site.
  p[j].charge=0;                             //neutralize proton...
  neutrons.push_back(j);                     //..and add to neutron vector

  action res;
  res.action=PROTONATED;
  res.site=i;
  res.proton=j;
  return res;
};

titrate::action titrate::moveToBulk(vector<particle> &p, short int i, short int j) {
  if (j==-1)             // if not specified..
    j=random(neutrons);  // ..pick random neutron
  neutrons.erase(remove(neutrons.begin(), neutrons.end(), j)); //del a neutron

  p[i].charge -= 1;      // deprotonate site.
  p[j].charge  = 1;      // charge up neutron...
  protons.push_back(j);  // ...add to proton vector

  action res;
  res.action=DEPROTONATED;
  res.site=i;
  res.proton=j;
  return res;
};

titrate::action titrate::exchange(vector<particle> &p) {
  int i=random( sites );        // pick a random, titrateable particle
  if (status(p,i)==PROTONATED)
    return moveToBulk(p,i);
  else
    return takeFromBulk(p,i);
};
titrate::action titrate::exchange(vector<particle> &p, action &a) {
  if (status(p, a.site)==PROTONATED)
    return moveToBulk(p, a.site, a.proton);
  else
    return takeFromBulk(p, a.site, a.proton);
};

double titrate::energy(vector<particle> &p, double du, action &a) {
  short unsigned int i=p[a.site].id;
  if (a.action==PROTONATED)
    return du+( log(10.)*(ph-d[i].pka) );
  else
    return du-( log(10.)*(ph-d[i].pka) );
};

double titrate::energy(vector<particle> &p, double du, double chempotCat, action &a) {
  short unsigned int i=p[a.site].id; //!< Energy function for GC tit.
  if (a.action==PROTONATED)          //!< OBS!!! Include the ideal part
                                     //!< as an argument to metropolis!
                                     //!< 
    return du+( log(10.)*(ph-d[i].pka) )+chempotCat;
  else 
    return du-( log(10.)*(ph-d[i].pka) )-chempotCat;
};

double titrate::idPref(action &a, double V) { //!< Returns the 
  int Na;                           //!< the GC prefactor, 
  Na=protons.size();                //!< N/V V/(N+1)
  if (a.action!=PROTONATED)
    return (V/double(Na));          //!< The proton vector is already 
  else                              //!< updated!!!
    return (double(Na+1)/V);
}

double titrate::check(action &a, double V, double chempotCat) { //!< Returns the
  int Na;                           //!< the GC prefactor, 
  Na=protons.size();                //!< N/V V/(N+1)
  if (a.action!=PROTONATED)
    return double(V/(Na)*exp(chempotCat));          //!< The proton vector is already 
  else                              //!< updated!!!
    return double((Na+1)/V*exp(-chempotCat));
}

void titrate::info() {
  cout << "# Titrateable sites   = " << sites.size() << endl
       << "# Protons             = " << protons.size() << endl
       << "# pH                  = " << ph << endl ;
};

double titrate::sumsites() {
  double Q=0;
  for (int i=0; i<sites.size(); i++)
    Q += q[i].avg();
  return Q;
};

void titrate::samplesites(vector<particle> &p) {
  for (int i=0; i<sites.size(); i++)
    q[i] += p[sites[i]].charge;
  nprot+=protons.size();
};

void titrate::showsites(vector<particle> &p) {
  cout << "# --- AVERAGE SITE CHARGES ---------------------\n";
  for (int i=0; i<d.size(); i++)
    if (d[i].pka!=0) {
      cout << d[i].name << ": " << setiosflags(ios::fixed);
      cout.precision(3);
      for (int j=0; j<sites.size(); j++)
        if (p[sites[j]].id==i)
          cout << q[j].avg() << " ";
      cout << endl;
    };
};

/*!
 * This function will take the average charges from titrate::q
 * and apply these to the specified particle vector. Also it will
 * smear out charges on the protons.
 * WARNING! After this function you can no longer perform any
 * titration steps. It is meant to be called before saving coordinates
 * to disk, so as to include the partial charges of the system.
 */
double titrate::applycharges(vector<particle> &p) {
  for (int i=0; i<sites.size(); i++)
    p[sites[i]].charge = q[i].avg();

  double q = nprot.avg() / protons.size();
  for (int i=0; i<protons.size(); i++)
    p[protons[i]].charge = q;
  
  cout << "# Charge, Titrateable sites = " << sumsites() << endl
       << "# Charge, Protons           = " << nprot.avg() << endl
       << "# Average proton charge     = " << q << endl
       << "# Total                     = " << sumsites()+nprot.avg() << endl;
  
  return sumsites();
};
