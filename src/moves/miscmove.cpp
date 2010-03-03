#include <faunus/moves/miscmove.h>
namespace Faunus {
  // Combined translation and rotation
  transrot::transrot( ensemble &e,
      container &c, energybase &i ) : markovmove(e,c,i)
  {
    name.assign("COMBINED ROTATION and TRANSLATION");
    runfraction=1.0;
    deltadp=0.1;
    dp=0;
    dpt=dpr=0.5;  
  };
  string transrot::info() {
    std::ostringstream o;
    o <<  markovmove::info()
      << "#   Dp (translation)      = " <<  dpt << endl
      << "#   Dp (rotation)         = " <<  dpr << endl;
    return o.str();
  } 

  /*!
   * \todo Cell overlap test missing, dp is used both for dr and the angle
   */
  double transrot::move(macromolecule &g) {
    du=0;
    cnt++;
    g.rotate(*con, dpr, dpt); //dpt in ang, dpr in rad./2
    //insert cell overlap test
    for (int i=g.beg; i<=g.end; i++) { 
      if (con->collision( con->trial[i] )==true) { 
        rc=HC;
        dpsqr+=0.;
        g.undo(*con);
        return du; 
      }
    }
    //#pragma omp parallel
    {
      //#pragma omp sections
      {
        //#pragma omp section
        { uold = pot->energy(con->p, g);   }
        //#pragma omp section
        { unew = pot->energy(con->trial, g);   }
      }
    }
    du = unew-uold; 
    
    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      naccept++;
      dpsqr+=con->sqdist(g.cm, g.cm_trial);
      g.accept(*con);
      return du;
    } else rc=ENERGY;
    du=0;
    dpsqr+=0.;
    g.undo(*con);
    return du;
  }

  /*!
   * \todo Cell overlap test missing, dp is used both for dr and the angle
   */
  double transrot::move(vector<macromolecule> &g, int n) {
    du=0;
    cnt++;
    g[n].rotate(*con, dpr, dpt); //dpt in ang, dpr in rad./2
    //insert cell overlap test
    for (int i=g[n].beg; i<=g[n].end; i++) { 
      if (con->collision( con->trial[i] )==true) { 
        rc=HC;
        dpsqr+=0.;
        g[n].undo(*con);
        return du; 
      }
    }
    double deltau=0;
#pragma omp parallel for reduction (+:deltau) 
    for (int i=0; i<g.size(); i++)
      if (i!=n)
        deltau += pot->energy(con->trial, g[i], g[n]) - pot->energy(con->p, g[i], g[n]);
    du=deltau;
    
    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      naccept++;
      dpsqr+=con->sqdist(g[n].cm, g[n].cm_trial);
      g[n].accept(*con);
      return du;
    } else rc=ENERGY;
    du=0;
    dpsqr+=0.;
    g[n].undo(*con);
    return du;
  }

  // Multiple paticle transrot
  multtr::multtr( ensemble &e,
      container &c, energybase &i, int M ) : markovmove(e,c,i)
  {
    name.assign("MULTIPLE COMBINED ROTATION and TRANSLATION");
    runfraction=1.0;
    deltadp=0.1;
    dp=1.0;
    dpt=dpr=0.5;  
    m=M;
  };

  string multtr::info() {
    std::ostringstream o;
    o <<  markovmove::info()
      << "#   Number of particles   = " <<  m   << endl
      << "#   Dp (translation)      = " <<  dpt << endl
      << "#   Dp (rotation)         = " <<  dpr << endl;
    return o.str();
  } 
  /*!
   * \todo Cell overlap test missing, dp is used both for dr and the angle
   */
  double multtr::move(molecules &g, vector<int> &n) {
    if (m!=n.size())
      std::cerr << "# MULTTR::MOVE Iterator out of sync with 'm'"<<endl;
    du=0;
    cnt++;
    //rotate
    macromolecule p;
    for (int i =0; i<m; i++) {
      p.beg=g[n[i]].beg;
      p.end=g[n[i]].end;
      p.rotate(*con, dpr, dpt); //dpt in ang, dpr in rad./2
    }
    //cast values 
    point q;
    for (int i=0; i<m; i++)
      for (int j=0; j<g.numatom; j++) {
        q.x=con->trial[g[n[i]].beg+j].x;
        q.y=con->trial[g[n[i]].beg+j].y;
        q.z=con->trial[g[n[i]].beg+j].z;
        con->trial[g[n[i]].beg+j].x=con->trial[g[i].beg+j].x;
        con->trial[g[n[i]].beg+j].y=con->trial[g[i].beg+j].y;
        con->trial[g[n[i]].beg+j].z=con->trial[g[i].beg+j].z;
        con->trial[g[i].beg+j].x=q.x;
        con->trial[g[i].beg+j].y=q.y;
        con->trial[g[i].beg+j].z=q.z;
        q.x=con->p[g[n[i]].beg+j].x;
        q.y=con->p[g[n[i]].beg+j].y;
        q.z=con->p[g[n[i]].beg+j].z;
        con->p[g[n[i]].beg+j].x=con->p[g[i].beg+j].x;
        con->p[g[n[i]].beg+j].y=con->p[g[i].beg+j].y;
        con->p[g[n[i]].beg+j].z=con->p[g[i].beg+j].z;
        con->p[g[i].beg+j].x=q.x;
        con->p[g[i].beg+j].y=q.y;
        con->p[g[i].beg+j].z=q.z;
      }
    p.beg=g[0].beg;
    p.end=g[m-1].end;
    //cell overlap test
    for (int i=p.beg; i<=p.end; i++) { 
      if (con->collision( con->trial[i] )==true) { 
        rc=HC;
        dpsqr+=0.;
        g.undo(*con);
        return du; 
      }
    }
    //#pragma omp parallel
    {
      //#pragma omp sections
      {
        //#pragma omp section
        { uold = pot->energy(con->p, g, n);   }
        //#pragma omp section
        { unew = pot->energy(con->trial, g, n);   }
      }
    }
    du = unew-uold; 
    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      naccept++;
      for (int i=0; i<n.size(); i++)
        dpsqr+=con->sqdist(g[i].masscenter(con->p), g[i].masscenter(con->trial));
      g.accept(*con);
      return du;
    } else rc=ENERGY;
    du=0;
    dpsqr+=0.;
    g.undo(*con);
    return du;
  }
}//namespace
