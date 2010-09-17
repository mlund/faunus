#include "faunus/moves/clustermove.h"
#include "faunus/ensemble.h"

namespace Faunus {
//  INVERSION  
  clusterinvw::clusterinvw( ensemble &e,
    container &c, sphericalimage<pot_test> &i ) : markovmove(e,c,i) {
      name.assign("WATER CLUSTER INVERSION, cite xxx");
      runfraction=1.;
      ipot=&i;
     }

  string clusterinvw::info() {
    ostringstream o;
    o << markovmove::info()
      << "#   A fraction of "<<movefrac.avg()<<" are moved ( "<<movefrac.stdev()<<")"<<endl;
    return o.str();
  }

  double clusterinvw::move(molecules &m) {
    cnt++;
    du=-pot->energy(con->p,m)-pot->internal(con->p, m, m.numatom);
    double udiff=0, un=0, uo=0;
      moved.clear(), remaining.clear();
    for (int i=0; i<m.size()/m.numatom; i++)
      remaining.push_back(i);
    int f;
    point origin;
    origin.x=0, origin.y=0, origin.z=0;
    f=slp.random_one()*32; //remaining.size();
    moved.push_back(remaining[f]);
    remaining.erase(remaining.begin()+f);    // Pick first index in m to move
    for(int i=0; i<moved.size(); i++) {
      m[moved[i]].invert(con->trial, origin);
      for(int j=0; j<remaining.size();j++) {
        uo=pot->energy(con->p,     m[moved[i]], m[remaining[j]]);
        un=pot->energy(con->trial, m[moved[i]], m[remaining[j]]);
        udiff=un-uo;
        if(slp.random_one() < (1.-exp(-udiff)) ) {
          moved.push_back(remaining[j]);
          remaining.erase(remaining.begin()+j);
          j=j-1;
        }
      }
      m[moved[i]].accept(*con);
      ipot->updateimg(con->p, m[moved[i]]);
    } 
    du+=pot->energy(con->p,m) + pot->internal(con->p,m, m.numatom);
    movefrac+=moved.size()/(moved.size()+remaining.size());
    rc=OK;
    naccept++;
    utot+=du;
    return du;   
  }

  clusterinvwsalt::clusterinvwsalt( ensemble &e,
    container &c, interaction<pot_hsminimage> &i ) : markovmove(e,c,i) {
      name.assign("SALT CLUSTER INVERSION, cite xxx");
      runfraction=1.;
      ipot=&i;
     }

  string clusterinvwsalt::info() {
    ostringstream o;
    o << markovmove::info()
      << "#   A fraction of "<<movefrac.avg()<<" are moved ( "<<movefrac.stdev()<<")"<<endl;
    return o.str();
  }

  double clusterinvwsalt::move(salt &m) {
    cnt++;
    du=-pot->internal(con->p, m);
    double udiff=0, un=0, uo=0;
      moved.clear(), remaining.clear();
    for (int i=m.beg; i<m.end+1; i++)
      remaining.push_back(i);
    int f;
    point origin, imgpos;
    con->randompos(origin);
    f=slp.random_one()*32; //remaining.size();
    moved.push_back(remaining[f]);
    remaining.erase(remaining.begin()+f);    // Pick first index in m to move
    imgpos.x=dp*slp.random_half();
    imgpos.y=dp*slp.random_half();
    imgpos.z=dp*slp.random_half();
    for(int i=0; i<moved.size(); i++) {
//#pragma omp parallel for reduction (+:udiff) schedule (dynamic)
      //imgpos=con->vdist(con->trial[moved[i]], origin);
      //con->trial[moved[i]]=-imgpos+origin;
      con->trial[moved[i]]=con->p[moved[i]]+imgpos;
      con->boundary(con->trial[moved[i]]);
      for(int j=0; j<remaining.size();j++) {
        udiff= -ipot->pair.pairpot(    con->p[moved[i]], con->p[remaining[j]])
               +ipot->pair.pairpot(    con->trial[moved[i]], con->trial[remaining[j]]);
        //udiff=un-uo;
        if(slp.random_one() < (1.-exp(-udiff)) ) {
          moved.push_back(remaining[j]);
          remaining.erase(remaining.begin()+j);
          j=j-1;
          dpsqr+=imgpos.len()*imgpos.len();
        }
      }
      con->p[moved[i]]=con->trial[moved[i]];
    } 
    du+=pot->internal(con->p,m);
    movefrac+=double(moved.size())/double((moved.size()+remaining.size()));
    rc=OK;
    naccept++;
    utot+=du;
    dpsqr+=0.0;
    return du;   
  }

/*  RESTRICTED INVERSION  
  clusterinvwsaltr::clusterinvwsaltr( ensemble &e,
    container &c, interaction<pot_hsminimage> &i, inputfile &in 
    ) : markovmove(e,c,i) {
      name.assign("RESTRICTED WATER CLUSTER INVERSION, cite xxx");
      runfraction=1.;
      ipot=&i;
      r=rad-3;
      cr=cavrad;
      cr2=cr*cr;
     }

  string clusterrinvw::info() {
    ostringstream o;
    o << markovmove::info()
      << "#   A fraction of "<<movefrac.avg()<<" are moved ( "<<movefrac.stdev()<<")"<<endl;
    return o.str();
  }

  double clusterrinvw::move(molecules &m) {
    cnt++;
    du=-pot->energy(con->p,m)-pot->internal(con->p, m, m.numatom);
    point ip;
    ip.ranunit(slp);
    ip=ip*r*slp.random_one();
    double udiff=0, un=0, uo=0;
      moved.clear(), remaining.clear();
    for (int i=0; i<m.size()/m.numatom; i++)
      if (con->sqdist(con->trial[m[i].beg], ip)<cr2)
        remaining.push_back(i);
    group g;
    for (int i=0; i<remaining.size(); i++) {
      g.beg=m[i].beg, g.end=m[i].end;
      m[remaining[i]].swap(*con, g);
    }
    int f;
    f=slp.random_one()*remaining.size();
    moved.push_back(remaining[f]);
    remaining.erase(remaining.begin()+f);    // Pick first index in m to move
    for(int i=0; i<moved.size(); i++) {
      m[moved[i]].invert(con->trial, ip);
      for(int j=0; j<remaining.size();j++) {
        uo=pot->energy(con->p,     m[moved[i]], m[remaining[j]]);
        un=pot->energy(con->trial, m[moved[i]], m[remaining[j]]);
        udiff=un-uo;
        if(slp.random_one() < (1.-exp(-udiff)) ) {
          moved.push_back(remaining[j]);
          remaining.erase(remaining.begin()+j);
          j=j-1;
        }
      }
//      m[moved[i]].accept(*con);
      ipot->updateimg(con->p, m[moved[i]]);
    } 
    du+=pot->energy(con->p,m) + pot->internal(con->p,m, m.numatom);
    movefrac+=moved.size()/(moved.size()+remaining.size());
    rc=OK;
    naccept++;
    utot+=du;
    return du;   
  }
*/
//  CLUSTER TRANSLATION  
  clustertrans::clustertrans( ensemble &e,
    container &c, energybase &i , vector<macromolecule> &g ) : markovmove(e,c,i) {
      name.assign("REJECTION FREE CLUSTER TRANSLATION");
      cite.assign("doi:10.1103/PhysRevLett.92.035504");
      runfraction=1.;
      skipEnergyUpdate=false;
      //g=&G;
      //distributions d(1., 1., g.size());
      //dist=d;
  }

  string clustertrans::info() {
    ostringstream o;
    o << markovmove::info()
      << "#   dp                        = " << dp <<endl
      << "#   Particle move fraction    = " << movefrac.avg() << " " << movefrac.stdev() << endl
      << "#   Skip energy updates:        "
      << ((skipEnergyUpdate==false) ? "no" : "yes (energy drift!)") << endl;
    return o.str();
  }

  void clustertrans::check(checkValue &t) {
    markovmove::check(t);
    t.check("clustertrans_movefraction", movefrac.avg() );
  }

  double clustertrans::move(vector<macromolecule> &g) {
    double du=0;
    if (slp.runtest(runfraction)==false)
      return du;
    cnt++;
    if (skipEnergyUpdate==false)
#pragma omp parallel for reduction (+:du) schedule (dynamic)
      for (int i=0; i<g.size()-1; i++)
        for (int j=i+1; j<g.size(); j++)
          du-=pot->energy(con->p, g[i], g[j]);
    point ip;
    ip.x=dp*slp.random_half();
    ip.y=dp*slp.random_half();
    ip.z=dp*slp.random_half();
    double udiff=0, un=0, uo=0;
    moved.clear(), remaining.clear();
    for (int i=0; i<g.size(); i++)
      remaining.push_back(i);
    int f;
    f=slp.random_one()*remaining.size();
    moved.push_back(remaining[f]);
    remaining.erase(remaining.begin()+f);    // Pick first index in m to move
    for(int i=0; i<moved.size(); i++) {
      g[moved[i]].move(*con, ip);
      for(int j=0; j<remaining.size();j++) {
        uo=pot->energy(con->p,     g[moved[i]], g[remaining[j]]);
        un=pot->energy(con->trial, g[moved[i]], g[remaining[j]]);
        udiff=un-uo;
        if(slp.random_one() < (1.-exp(-udiff)) ) {
          moved.push_back(remaining[j]);
          remaining.erase(remaining.begin()+j);
          j=j-1;
        }
      }
      g[moved[i]].accept(*con);
    }
    if (skipEnergyUpdate==false)
#pragma omp parallel for reduction (+:du) schedule (dynamic)
      for (int i=0; i<g.size()-1; i++)
        for (int j=i+1; j<g.size(); j++)
          du+=pot->energy(con->p, g[i], g[j]);
    movefrac+=double(moved.size())/double((moved.size()+remaining.size()));
    rc=OK;
    naccept++;
    utot+=du;
    return du;   
  }

  /*/ TRANSLATION
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
*/
  //---------- ROTATE CLUSTER OF MACROMOLECULE ----------------
  clusterrotate::clusterrotate( ensemble &e,
      container &c, energybase &i ) : markovmove(e,c,i) {
    name = "MACROMOLECULAR CLUSTER ROTATION";
    runfraction=0.8;
    deltadp=1.;
    dp=0.5;
    nacc.assign(30,0);
    ntrial.assign(30,0);
  }

  void clusterrotate::print() {
    std::cout <<" Index's in 'cluster': ";
    for (int i=0; i<cluster.size(); i++)
      std::cout <<cluster[i]<<", ";
    std::cout <<endl;
    std::cout <<" Index's in 'free': ";
    for (int i=0; i<free.size(); i++)
      std::cout <<free[i]<<", ";
    std::cout <<endl<<endl;
  }

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

/*    if (cluster.size()>1) {
      std::cout <<" Distance between centers of mass before rotation"<<endl;
      for (short i=0; i<iend-1; i++) {
        point cmp;
        for (short j=i+1; j<iend; j++) {
          std::cout<<con->dist(g[cluster[i]].cm,g[cluster[j]].cm)<<"   ";
          //std::cout <<cmp<<"   ";
        }
      }
      std::cout<<endl;
    }*/


    for (short i=0; i<iend; i++)
      g[cluster[i]].rotate(*con, cr, p, angle );
    flowcheck(g);
    for (short i=0; i<iend; i++) {
      if (con->collision(g[cluster[i]].cm_trial)==true || FLOW!=0) {  //   || cluster.size()>3) {
        rc=HC;
        dpsqr+=0;
        for (short i=0; i<iend; i++)
          g[cluster[i]].undo(*con);
        return du;
      }
    }
    double Unew, Uold, dUcheck;
    Unew=Uold=dUcheck=0.0;
    for (short i=0; i<iend; i++)
      for (short j=0; j<jend; j++) {
        Uold+= pot->energy( con->p, g[cluster[i]], g[free[j]] );
        Unew+= pot->energy( con->trial, g[cluster[i]], g[free[j]] );
      }
    du = Unew-Uold;
/*    if (cluster.size()>1) {
      std::cout <<" Distance between centers of mass after rotation"<<endl;
      for (short i=0; i<iend-1; i++) {
        point cmp;
        for (short j=i+1; j<iend; j++) {
          dUcheck+=pot->energy( con->trial, g[cluster[i]], g[cluster[j]] ) - pot->energy( con->p, g[cluster[i]], g[cluster[j]] ) ;
          std::cout<<con->dist(g[cluster[i]].cm_trial,g[cluster[j]].cm_trial)<<"   ";
          //std::cout <<cmp<<"   ";
        }
      }
      std::cout<<endl;
      std::cout <<" Energy difference in cluster = "<<dUcheck<<endl<<endl;
    }*/
    



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

/*  bool clusterrotate::reversible(point &p, group &g) {
    for (int i=g.beg; i<=g.end; i++)

  }*/

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

    int seed=slp.random_one()*nprot;            //pick seed
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
    o << "#   CLUSTER SIZES || Cluster separation is set to "<<sep<<" AA"<<endl
      << "# --------------------------------------------------------------" <<endl;
    for (short int i=0; i<30; i++) {
      if (ntrial[i]!=0)
        o << "#  "<<i<< "-particle cluster  : "<< 100*nacc[i]/ntrial[i] << " % acc., "<<ntrial[i]<< " attempts." <<endl
          << "# -------------------------------------------------------------- " <<endl;
    }
    return o.str();
  }
 /*  clusterinv::clusterinv( ensemble &e, container &c, interaction<T_pairpot> &i) : markovmove(e,c,i) {
    name = "REJECTION FREE INVERTATIONS, Explicitly writen for spce water!!!";
    runfraction=0;
  }
  double clusterinv::move(molecules &m) {
    udiff=0;
    remaining.clear();
    moved.clear();
    for (int i=0; i<m.size/m.numatom; i++)
      remaining.push_back(i);
    int f;
    f=slp.random_one()*remaining.size();
    moved.push_back(remaining[f]);
    remaining.erase(remaining.begin()+int(f-1));    // Pick first index in m to move

    for(int i=0; i<moved.size(); i++) {
      con->p[m[i].beg].x=-con->p[m[i].beg].x
        con->p[m[i].beg].y=-con->p[m[i].beg].y
        con->p[m[i].beg].z=-con->p[m[i].beg].z
        con->p[m[i].beg+1].x=-con->p[m[i].beg+1].x
        con->p[m[i].beg+1].y=-con->p[m[i].beg+1].y
        con->p[m[i].beg+1].z=-con->p[m[i].beg+1].z
        con->p[m[i].beg+2].x=-con->p[m[i].beg+2].x
        con->p[m[i].beg+2].y=-con->p[m[i].beg+2].y
        con->p[m[i].beg+2].z=-con->p[m[i].beg+2].z
        for(int j=0; i<moved.size(); j++) {
          udiff=pot->energy(con->p, m[i], m[j]) - pot->energy(con->trial, m[i], m[j]);
          if(slp.random_one() < (1.-exp(-udiff))) {
            con->p[m[j].beg].x=-con->p[m[j].beg].x
              con->p[m[j].beg].y=-con->p[m[j].beg].y
              con->p[m[j].beg].z=-con->p[m[j].beg].z
              con->p[m[j].beg+1].x=-con->p[m[j].beg+1].x
              con->p[m[j].beg+1].y=-con->p[m[j].beg+1].y
              con->p[m[j].beg+1].z=-con->p[m[j].beg+1].z
              con->p[m[j].beg+2].x=-con->p[m[j].beg+2].x
              con->p[m[j].beg+2].y=-con->p[m[j].beg+2].y
              con->p[m[j].beg+2].z=-con->p[m[j].beg+2].z

              flytta over mellan vektorer...

          } 
        }
    }
    m[i].accept(*con)
  }*/
}//namespace
