#include "faunus/moves/membrane.h"
#include "faunus/io.h"

namespace Faunus {
  membranemove::membranemove(ensemble &e, container &c, energybase &i, popscmembrane &m) 
  : markovmove(e,c,i) {
    name.assign("POPS/C MEMBRANE-CHAIN");
    //Balance between the moves and set up vectors
    endacc=lattacc=popsacc=popcacc=0;
    endcnt=lattcnt=popscnt=popccnt=0;
    double sum;
    sum=2.*m.pops.size()+1.*m.popc.size();
    rfend =m.pops.size()/sum; // disabled due to error energy drift ->, rflatt=(m.pops.size()+m.popc.size())/sum;
    rfpopc=m.popc.size()/sum, rfpops =rfend;
    steps =m.popc.size()+m.pops.size();
    dpv.x =1.0, dpv.y =1.0, dpv.z =1.0;
    dpvg.x=1.0, dpvg.y=1.0, dpvg.z=0.0;
  }

  string membranemove::info() {
    std::ostringstream o;
    o << markovmove::info()
      << "#   Displacement parameters (A):"<<endl
      << "#      Monomer_dp         = "<<dpmm<<endl
      << "#      Graftpoint_dp      = "<<dpgp<<endl
      << "#      Latticetrans._dp   = "<<dplatt<<endl;
    o.precision(3);
    o << "#   Acceptance of substeps:"<<endl
      << "#      Pops rattling      = "<<popsacc/popscnt*100.<< " \%" << endl
      << "#      Popc rattling      = "<<popcacc/popccnt*100.<< " \%" << endl
      << "#      End jumping        = "<<endacc/endcnt*100.<< " \%" << endl;
    if (lattcnt>0)
      o << "#      Latt. translation  = "<<lattacc/lattcnt*100.<<endl;
    return o.str();
  }

  double membranemove::move(popscmembrane &m) {
    utemp=du=0;
    int n=-1, pm, pm2, i=0;
    double rand, sum=0;
    point dv;
//    double backreq;
//    backreq=pot->req, pot->req=5.0;
    while (i<steps) {
      rand=slp.random_one();
      cnt++;
      i++;
      if(rand<rfpops) {//Pops
        popscnt++;
        pm=m.pops.size()*slp.random_one();  //Pick a random polymer
        n=m.pops[pm].beg+slp.random_one()*4.;
        if (n==m.pops[pm].beg)
          dv=dpvg*dpgp;
        else
          dv=dpv*dpmm;
        con->trial[n].x = con->p[n].x + dv.x*slp.random_half();
        con->trial[n].y = con->p[n].y + dv.y*slp.random_half();
        con->trial[n].z = con->p[n].z + dv.z*slp.random_half();
        if (con->collision( con->trial[n] )==true) {
          rc=HC;
          con->trial[n] = con->p[n];
        } else {
          uold = pot->u_monomer(con->p, m.pops[pm], n);   
          unew = pot->u_monomer(con->trial, m.pops[pm], n);   
          du = unew - uold;
          if (ens->metropolis(du)==true) {
            rc=OK;
            utot+=du;                               // track energy changes
            double d2=con->sqdist(con->p[n],con->trial[n]); // track avg. displacement
            dpsqr+=d2;
            //rsqr+=d2/pow(g.size(),2);               // Track mean square displacement per particle
            naccept++;                              // accept counter
            con->p[n] = con->trial[n];              // Accept move
            m.pops[pm].masscenter(*con);                     // Recalc mass center
            utemp+=du;
            popsacc++;
          } else {
            rc=ENERGY;
            dpsqr+=0;
            con->trial[n] = con->p[n];
          }
        }
        du=0;
      }
      if(rand>rfpops && rand<rfpops+rfpopc) {//Popc
        popccnt++;
        pm=m.popc.size()*slp.random_one();  //Pick a random polymer
        n=m.popc[pm].beg+slp.random_one()*3.;
        if (n==m.popc[pm].beg)
          dv=dpvg*dpgp;
        else
          dv=dpv*dpmm;
        con->trial[n].x = con->p[n].x + dv.x*slp.random_half();
        con->trial[n].y = con->p[n].y + dv.y*slp.random_half();
        con->trial[n].z = con->p[n].z + dv.z*slp.random_half();
        if (con->collision( con->trial[n] )==true) {
          rc=HC;
          con->trial[n] = con->p[n];
        } else {
          uold = pot->u_monomer(con->p, m.popc[pm], n);   
          unew = pot->u_monomer(con->trial, m.popc[pm], n);   
          du = unew - uold;
          if (ens->metropolis(du)==true) {
            rc=OK;
            utot+=du;                               // track energy changes
            double d2=con->sqdist(con->p[n],con->trial[n]); // track avg. displacement
            dpsqr+=d2;
            //rsqr+=d2/pow(g.size(),2);             // Track mean square displacement per particle
            naccept++;                              // accept counter
            con->p[n] = con->trial[n];              // Accept move
            m.popc[pm].masscenter(*con);            // Recalc mass center
            utemp+=du;
            popcacc++;
          } else {
            rc=ENERGY;
            dpsqr+=0;
            con->trial[n] = con->p[n];
          }
        }
        du=0;
      }
      if(rand>rfpops+rfpopc && rand<rfpops+rfpopc+rfend) {//End jumps
        endcnt++;
        point endp;
        pm  =m.pops.size()*slp.random_one();
        pm2 =m.popc.size()*slp.random_one();
        endp=con->vdist(con->p[m.pops[pm].end], con->p[m.pops[pm].beg+2]);
        con->trial[m.pops[pm].beg]   =con->p[m.popc[pm2].beg];
        con->trial[m.pops[pm].beg+1] =con->p[m.popc[pm2].beg+1];
        con->trial[m.pops[pm].beg+2] =con->p[m.popc[pm2].beg+2];
        con->trial[m.pops[pm].end]   =con->p[m.popc[pm2].beg+2]+endp;
        con->trial[m.popc[pm2].beg]  =con->p[m.pops[pm].beg];
        con->trial[m.popc[pm2].beg+1]=con->p[m.pops[pm].beg+1];
        con->trial[m.popc[pm2].beg+2]=con->p[m.pops[pm].beg+2];
        if (con->collision( con->trial[m.pops[pm].end] )==true) {
          rc=HC;
          con->trial[m.pops[pm].beg]   =con->p[m.pops[pm].beg];  
          con->trial[m.pops[pm].beg+1] =con->p[m.pops[pm].beg+1];
          con->trial[m.pops[pm].beg+2] =con->p[m.pops[pm].beg+2];
          con->trial[m.pops[pm].end]   =con->p[m.pops[pm].end];
          con->trial[m.popc[pm2].beg]  =con->p[m.popc[pm2].beg];   
          con->trial[m.popc[pm2].beg+1]=con->p[m.popc[pm2].beg+1];
          con->trial[m.popc[pm2].beg+2]=con->p[m.popc[pm2].beg+2];
        } else {
//          uold = pot->u_monomer(con->p, m.pops[pm], m.pops[pm].end);   
//          unew = pot->u_monomer(con->trial, m.pops[pm], m.pops[pm].end);   
          uold = pot->energy(con->p, m.pops[pm].end);   
          unew = pot->energy(con->trial, m.pops[pm].end);   
          du = unew - uold;
//          std::cout <<unew<<"  "<<uold<<"  "<<du<<std::endl;
//          std::cout <<con->p[m.pops[pm].end]<<"  "<<con->trial[m.pops[pm].end]<<std::endl;
//          std::cout <<pot->energy(con->p, m.pops[pm].end)<<"  "<<pot->energy(con->trial, m.pops[pm].end)<<std::endl<<std::endl;
          if (ens->metropolis(du)==true) {
            rc=OK;
            utot+=du;                               // track energy changes
            double d2=con->sqdist(con->p[m.pops[pm].end],con->trial[m.pops[pm].end]); // track avg. displacement
            dpsqr+=d2;
            //rsqr+=d2/pow(g.size(),2);               // Track mean square displacement per particle
            naccept++;                              // accept counter
            con->p[m.pops[pm].beg]   =con->trial[m.pops[pm].beg];  
            con->p[m.pops[pm].beg+1] =con->trial[m.pops[pm].beg+1];
            con->p[m.pops[pm].beg+2] =con->trial[m.pops[pm].beg+2];
            con->p[m.pops[pm].end]   =con->trial[m.pops[pm].end];
            con->p[m.popc[pm2].beg]  =con->trial[m.popc[pm2].beg];   
            con->p[m.popc[pm2].beg+1]=con->trial[m.popc[pm2].beg+1];
            con->p[m.popc[pm2].beg+2]=con->trial[m.popc[pm2].beg+2];
            utemp+=du;
            m.popc[pm2].masscenter(*con);            // Recalc mass center
            m.pops[pm].masscenter(*con);             // Recalc mass center
            endacc++;
          } else {
            rc=ENERGY;
            dpsqr+=0;
            con->trial[m.pops[pm].beg]   =con->p[m.pops[pm].beg];  
            con->trial[m.pops[pm].beg+1] =con->p[m.pops[pm].beg+1];
            con->trial[m.pops[pm].beg+2] =con->p[m.pops[pm].beg+2];
            con->trial[m.pops[pm].end]   =con->p[m.pops[pm].end];
            con->trial[m.popc[pm2].beg]  =con->p[m.popc[pm2].beg];   
            con->trial[m.popc[pm2].beg+1]=con->p[m.popc[pm2].beg+1];
            con->trial[m.popc[pm2].beg+2]=con->p[m.popc[pm2].beg+2];
          }
        }
        du=0;
      }
      if(rand-1>rfpops+rfpopc+rfend) {// Translate polymers
        lattcnt++;                    // This does NOT work, gives rise to energy drift
        group tempp;
        point slptrans;
        slptrans.x=dplatt*slp.random_half();
        slptrans.y=dplatt*slp.random_half();
        slptrans.z=0;
        if(slp.random_one()>double(m.pops.size())/double(m.pops.size()+m.popc.size())) { //Popc it is
          tempp=m.popc[slp.random_one()*m.popc.size()];
        } else {                                                                       //Pops it is
          tempp=m.pops[slp.random_one()*m.pops.size()];
        }
        tempp.move(*con, slptrans);
        uold = pot->energy(con->p, tempp);   
        unew = pot->energy(con->trial, tempp);   
        du = unew - uold;
        if (ens->metropolis(du)==true) {
          rc=OK;
          utot+=du;                                // Track energy changes
          double d2=con->sqdist(con->p[n],con->trial[n]); // Track avg. displacement
          dpsqr+=d2;
          naccept++;                               // accept counter
          tempp.accept(*con);                      // Accept move
          utemp+=du;
          lattacc++;                               // Local counter
        } else {
          rc=ENERGY;
          dpsqr+=0;
          tempp.undo(*con);
        }
        du=0;
      }
    }
//    pot->req=backreq;
    return utemp;
  }
}//namespace
