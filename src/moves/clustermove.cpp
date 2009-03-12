#include <faunus/moves/clustermove.h>

namespace Faunus {
  clustertranslate::clustertranslate( ensemble &e,
      container &c, energybase &i, double S ) : markovmove(e,c,i) {
    name = "MACROMOLECULAR CLUSTER TRANSLATION";
    runfraction=0.5;
    deltadp=1.;
    dp=10.;
    sep=S;
    nacc.assign(30,0);
    ntrial.assign(30,0);
  }

  double clustertranslate::move(vector<macromolecule> &g) {
    du=0;
    if (slp.runtest(runfraction)==false)
      return du;
    unew=0;
    uold=0;
    FLOW=0;
    cnt++;
    point p;
    p.x=dp*slp.random_half();
    p.y=dp*slp.random_half();
    p.z=dp*slp.random_half();
    sort(g);

    short int iend=cluster.size();
    short int jend=free.size();

    ntrial[iend]++;

    for (short i=0; i<iend; i++)
      g[cluster[i]].move(*con, p);
    flowcheck(g);
    for (short i=0; i<iend; i++) {
      if (con->collision(g[i].cm_trial)==true || FLOW!=0) {
        rc=HC;
        dpsqr+=0;
        for (short i=0; i<iend; i++)
          g[i].undo(*con);
        return du;
      }
    }
    double Uold;
    double Unew;
    for (short i=0; i<iend; i++)
      for (short j=0; j<jend; j++) {
        Uold+= pot->energy( con->p, g[cluster[i]], g[free[j]] );
        Unew+= pot->energy( con->trial, g[cluster[i]], g[free[j]] );
      }
    du = Unew-Uold;
    if (ens->metropolis(du)==true ) {
      rc=OK;
      utot+=du;
      dpsqr+=con->sqdist( g[0].cm, g[0].cm_trial );
      naccept++;
      nacc[iend]++;
      for (short i=0; i<iend; i++)
        g[cluster[i]].accept(*con);
      return du;
    } else rc=ENERGY;
    du=0;
    dpsqr+=0;
    for (short i=0; i<iend; i++)
      g[cluster[i]].undo(*con);
    return du;
  }

  bool clustertranslate::decreasing(int i, int j) {
    return (i>j);
  } 

  void clustertranslate::flowcheck(vector<macromolecule> &g) {
    short int iend=cluster.size();
    short int jend=free.size();
    for (short i=0; i<iend; i++)
      for (short j=0; j<jend; j++) {
        if(coll.overlap( con->p, g[cluster[i]], g[free[j]], sep))
          FLOW++;
      }
  }

  void clustertranslate::sort(vector<macromolecule> &g) {
    int nprot=g.size();
    free.clear();
    cluster.clear();
    for (short i=0; i<nprot; i++) {
      free.push_back(i);
    }

    int seed=slp.random_one()*nprot;              //pick seed
    cluster.push_back(seed);                    //place it in cluster
    free.erase(free.begin() + seed);            //and remove it from the free
    for (short s=0; s<cluster.size(); s++)
      for (int i=0; i<free.size(); i++) {
        if (coll.overlap(con->p, g[cluster[s]] ,g[free[i]] ,sep)) {
          cluster.push_back(free[i]);
          free.erase(free.begin() + i);
          i--;  
        }
      }
    if (nprot != cluster.size()+free.size())
      cout << "CLUSTER SORT HAS FAILED"<<endl;
  }

  string clustertranslate::info() {
    ostringstream o;
    o << markovmove::info();
    o << "#   CLUSTER SIZES "<<endl
      << "# --------------------------------------------------------------" <<endl;
    for (short int i=0; i<30; i++) {
      if (ntrial[i]!=0)
        o << "#  "<<i<< "-particle cluster  : "<< 100*nacc[i]/ntrial[i] << " % acc., "<<ntrial[i]<< " attempts." <<endl
          << "# -------------------------------------------------------------- " <<endl;
    }
    return o.str();
  }

  //---------- ROTATE CLUSTER OF MACROMOLECULE ----------------
  clusterrotate::clusterrotate( ensemble &e,
      container &c, interaction<T_pairpot> &i ) : markovmove(e,c,i) {
    name = "MACROMOLECULAR CLUSTER ROTATION";
    runfraction=0.8;
    deltadp=1.;
    dp=0.5;
    nacc.assign(30,0);
    ntrial.assign(30,0);
  };

  double clusterrotate::move(vector<macromolecule> &g) {
    du=0;
    uold=0;
    unew=0;
    FLOW=0;
    cnt++;
    point p,cr;
    p.ranunit(slp);
    sort(g);
    angle=slp.random_half()*dp;
    cr=g[cluster[0]].cm;

    short int iend=cluster.size();
    short int jend=free.size();

    ntrial[iend]++;

    for (short i=0; i<iend; i++)
      g[cluster[i]].rotate(*con, cr, p, angle );
    flowcheck(g);
    for (short i=0; i<iend; i++) {
      if (con->collision(g[i].cm_trial)==true || FLOW!=0 || cluster.size()>3) {
        rc=HC;
        dpsqr+=0;
        for (short i=0; i<iend; i++)
          g[i].undo(*con);
        return du;
      }
    }
    double Unew, Uold;
    //  #pragma omp parallel for reduction (+:Uold) reduction (+:Unew) num_threads(4)
    for (short i=0; i<iend; i++)
      for (short j=0; j<jend; j++) {
        Uold+= pot->energy( con->p, g[cluster[i]], g[free[j]] );
        Unew+= pot->energy( con->trial, g[cluster[i]], g[free[j]] );
      }
    du = Unew-Uold;
    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      dpsqr+=con->sqdist( g[0].cm, g[0].cm_trial );
      naccept++;
      nacc[iend]++;
      for (short i=0; i<iend; i++)
        g[cluster[i]].accept(*con);
      return du;
    } else rc=ENERGY;
    du=0;
    dpsqr+=0;
    for (short i=0; i<iend; i++)
      g[cluster[i]].undo(*con);
    return du;
  }

  bool clusterrotate::decreasing(int i, int j) {
    return (i>j);
  } 

  void clusterrotate::flowcheck(vector<macromolecule> &g) {
    short int iend=cluster.size();
    short int jend=free.size();
    for (short i=0; i<iend; i++)
      for (short j=0; j<jend; j++) {
        if(coll.overlap( con->p, g[cluster[i]], g[free[j]], sep))
          FLOW++;
      }
  }

  void clusterrotate::sort(vector<macromolecule> &g) {
    int nprot=g.size();
    free.clear();
    cluster.clear();
    for (short i=0; i<nprot; i++) {
      free.push_back(i);
    }

    int seed=slp.random_one()*nprot;              //pick seed
    cluster.push_back(seed);                    //place it in cluster
    free.erase(free.begin() + seed);            //and remove it from the free
    for (short s=0; s<cluster.size(); s++)
      for (int i=0; i<free.size(); i++) {
        if (coll.overlap(con->p, g[cluster[s]] ,g[free[i]] ,sep)) {
          cluster.push_back(free[i]);
          free.erase(free.begin() + i);
          i--;  
        }
      }
    if (nprot != cluster.size()+free.size())
      cout << "CLUSTER SORT HAS FAILED"<<endl;
  }

  string clusterrotate::info() {
    ostringstream o;
    o << markovmove::info();
    o << "#   CLUSTER SIZES "<<endl
      << "# --------------------------------------------------------------" <<endl;
    for (short int i=0; i<30; i++) {
      if (ntrial[i]!=0)
        o << "#  "<<i<< "-particle cluster  : "<< 100*nacc[i]/ntrial[i] << " % acc., "<<ntrial[i]<< " attempts." <<endl
          << "# -------------------------------------------------------------- " <<endl;
    }
    return o.str();
  }
}//namespace
