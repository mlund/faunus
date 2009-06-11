#include "faunus/titrate.h"
namespace Faunus {
  /*!
   * \param atom Species class
   * \param peeage pH
   */
  titrate::titrate(double peeage) {
    ph=peeage;
  }
  titrate::titrate(vector<particle> &p, group &mobile, double peeage) {
    ph=peeage;
    init(p, mobile);
  }

  /*!
   * Search for titrateable sites and locate "protons" and "neutrons"
   * in the particle vector.
   *
   * \note This function will de-protonate ALL titratable sites if no "neutrons" are found.
   *       Usefull if class is called after loading a prviously saved configuration.
   */
  void titrate::init(vector<particle> &p, group &mobile) {
    unsigned short i;
    for (i=mobile.beg; i<mobile.end+1; i++) {
      if (p[i].charge==+1) protons.push_back(i);
      if (abs(p[i].charge)<1e-6) neutrons.push_back(i);
    }
    for (i=0; i<p.size(); i++)    //search for titrateable sites
      if ( atom[p[i].id].pka !=0 ) {
        sites.push_back(i);
        if ( neutrons.size()==0 )
          p[i].charge = atom[p[i].id].charge; // deprotonate everything
      }
    q.resize( sites.size() ); //adjust site average charge vector
  }

  //returns number a random, titrateable
  //site in particle vector
  short int titrate::random(vector<short int> &iv) {
    short int i = rand() % iv.size();     //and pick randomly
    return iv.at(i);
  }

  //Returns protonation status of particle i in particle vector
  titrate::keywords titrate::status(vector<particle> &p, short int i) {
    if ( p[i].charge > atom[p[i].id].charge )
      return PROTONATED;
    return DEPROTONATED;
  }

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
  }

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
  }

  titrate::action titrate::exchange(vector<particle> &p) {
    int i=random( sites );        // pick a random, titrateable particle
    if (status(p,i)==PROTONATED)
      return moveToBulk(p,i);
    else
      return takeFromBulk(p,i);
  }
  titrate::action titrate::exchange(vector<particle> &p, action &a) {
    if (status(p, a.site)==PROTONATED)
      return moveToBulk(p, a.site, a.proton);
    else
      return takeFromBulk(p, a.site, a.proton);
  }

  double titrate::energy(vector<particle> &p,
      double du, action &a) {
    int i=p[a.site].id;
    if (a.action==PROTONATED)
      return du+( log(10.)*( ph - atom[i].pka ) );
    else
      return du-( log(10.)*( ph - atom[i].pka ) );
  }

  void titrate::infos() {
    std::cout << "# Titrateable sites   = " << sites.size() << endl
      << "# Protons             = " << protons.size() << endl
      << "# pH                  = " << ph << endl ;
  }

  double titrate::sumsites() {
    double Q=0;
    for (unsigned int i=0; i<sites.size(); i++)
      Q+= q[i].avg();
    return Q;
  }

  void titrate::samplesites(vector<particle> &p) {
    for (unsigned int i=0; i<sites.size(); i++)
      q[i] += p[sites[i]].charge;
    nprot+=protons.size();
  }

  /*
   * Returns the average charge of a titratable site. If
   * the chosen particle, i, is not titratable its charge
   * from the particle vector will be returned instead.
   */
  double titrate::avgcharge(vector<particle> &p, unsigned int i) {
    for (unsigned int j=0; j<sites.size(); j++)
      if (sites[j]==i) return q[j].avg();
    return p[i].charge;
  }

  void titrate::showsites(vector<particle> &p) {
    std::cout << "# --- AVERAGE SITE CHARGES ---------------------\n";
    for (unsigned int i=0; i<atom.list.size(); i++)
      if ( atom[i].pka!=0) {
        std::cout << atom[i].name << ": " << setiosflags(std::ios::fixed);
        std::cout.precision(3);
        for (unsigned int j=0; j<sites.size(); j++)
          if (p[sites[j]].id==i)
            std::cout << q[j].avg() << " ";
        std::cout << endl;
      };
  }

  /*!
   * This function will take the average charges from titrate::q
   * and apply these to the specified particle vector. Also it will
   * smear out charges on the protons.
   * WARNING! After this function you can no longer perform any
   * titration steps. It is meant to be called before saving coordinates
   * to disk, so as to include the partial charges of the system.
   */
  double titrate::applycharges(vector<particle> &p) {
    for (unsigned int i=0; i<sites.size(); i++)
      p[sites[i]].charge = q[i].avg();

    double q = nprot.avg() / protons.size();
    for (unsigned int i=0; i<protons.size(); i++)
      p[protons[i]].charge = q;

    std::cout << "# Charge, Titrateable sites = " << sumsites() << endl
      << "# Charge, Protons           = " << nprot.avg() << endl
      << "# Average proton charge     = " << q << endl
      << "# Total                     = " << sumsites()+nprot.avg() << endl;
    return sumsites();
  }

  //-------------------------------------------

  titrate_implicit::titrate_implicit(vector<particle> &p, double peeage, double muH) {
    ph=peeage;
    mu_proton=muH;
    for (int i=0; i<p.size(); i++)    //search for titrateable sites
      if ( atom[p[i].id].pka !=0 ) {
        sites.push_back(i);
        p[i].charge = atom[p[i].id].charge; // deprotonate everything
      }
  }

  int titrate_implicit::exchange(vector<particle> &p, int i) {
    if (i<0) i=random();
    double Zdp =atom[p[i].id].charge; // deprotonated charge of residue
    if ( p[i].charge > Zdp ) {
      p[i].charge=Zdp;
      recent = titrate_implicit::DEPROTONATED;
    }
    else {
      p[i].charge=Zdp+1.;
      recent = titrate_implicit::PROTONATED;
    }
    return i;
  }

  unsigned int titrate_implicit::random() {
    unsigned int i = rand() % sites.size();
    return sites.at(i);
  }
  double titrate_implicit::energy(vector<particle> &p, double du, int j) {
    int i=p[j].id;
    if (recent==PROTONATED)
      return du+( log(10.)*( ph - atom[i].pka ) ) - mu_proton;
    else
      return du-( log(10.)*( ph - atom[i].pka ) ) + mu_proton;
  }
}//namespace
