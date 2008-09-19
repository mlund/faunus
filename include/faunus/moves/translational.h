#ifndef FAU_TRANSLATIONAL_H
#define FAU_TRANSLATIONAL_H
#include "faunus/moves/base.h"
namespace Faunus {

  //-----------Z-MOVE--------------------------------------
  /*! \brief Class to make moves of a single particle along z-axis
   *  
   *  This class will make a random walk along the z-axis for a
   *  macromolece
   *
   *  \todo Needs some testing and perhaps some optimization
   *  \note Replaced by dualmove(?)
   *  \author Bjoern Persson
   */
  class zmove : public markovmove {
    public:
      zmove( ensemble&, container&, energybase&, macromolecule&, float);
      bool move(macromolecule &);
      float z;
  };

  zmove::zmove( ensemble &e, container &c, energybase &i, macromolecule &g,float rf) : markovmove(e,c,i) {
    runfraction=rf;
    dp=8;
    deltadp=1;
    name="MACROMOLECULE Z-DISPLACEMENTS";
  }     

  bool zmove::move(macromolecule &g) {
    du=0;
    if (slp.runtest(runfraction)==false)
      return false;
    cnt++;
    z=2*dp*slp.random_half();
    g.zmove(*con, z);
    for (int i=g.beg; i<(g.size()+g.beg); i++) { 
      if (con->collision( con->trial[i] )==true) 
        rc=HC;
    }
    if (rc==HC) {
      for (int i=g.beg; i<(g.size()+g.beg); i++) 
        con->trial[i] = con->p[i];
      return false; }
    else {
      //#pragma omp parallel
      {
        //#pragma omp sections
        {
          //#pragma omp section
          { uold = pot->energy(con->p, g);   }
          //#pragma omp section
          { unew = pot->energy(con->trial,g);   }
        }
      }
      du += unew - uold; }
      if (ens->metropolis(du)==true) {
        rc=OK;
        utot+=du;
        naccept++;
        for (int i=g.beg; i<(g.beg+g.size()); i++)
          con->p[i] = con->trial[i];
        return true;
      } else rc=ENERGY;

      du=0;
      for (int i=g.beg; i<(g.size()+g.beg); i++) 
        con->trial[i] = con->p[i];
      return false;
  }

  //----------------- DUAL MOVE --------------------------
  /*!
   * This move will symmetrically translate two macromolecules
   * along the line connecting their mass centers.
   *
   * \brief Symmetrically move two groups along z-axis
   * \author Mikael Lund
   */
  class dualmove : public markovmove {
    private:
      point v;
    public:
      histogram gofr;             //!< g(r) of the two group mass centers
      double r;                   //!< Current distance between group mass centers
      double rmax;                //!< Maximum allowed mass-center distance
      double rmin;                //!< Minimum allowed mass-center distance
      dualmove( ensemble&, container&, energybase&);
      void setup(inputfile &);
      void load(inputfile &, vector<macromolecule> &g, float=0);
      void direction(double, double, double);
      double move(macromolecule &, macromolecule &);
      string info();
  };
  dualmove::dualmove( ensemble &e,
      container &c, energybase &i ) : markovmove(e,c,i), gofr(0.1,0.,100.)
  {
    name = "SYMMETRIC 1D GROUP TRANSLATION";
    cite = "Biophys J. 2003, 85, 2940";
    runfraction=1.0;
    deltadp=1.;
    dp=1.;
    direction(0,0,1); // move only in z direction
    rmin=0;
    rmax=pow(double(c.getvolume()), 1./3)/2.; // rough estimate from volume
  }

  /*! Load the following parameters from an input
   * object: displacement (dm_dp), minimum separation
   * (dm_minsep), maximum separation (dm_maxsep)
   */
  void dualmove::setup( inputfile &in ) {
    dp = in.getflt("dm_dp",dp);
    rmin = in.getflt("dm_minsep",rmin);
    rmax = in.getflt("dm_maxsep",rmax);
  }

  /*! Specify unit vector that determines which coordinates
   *  that will be moved. The default is (0,0,1) meaning that
   *  the macromolecules will be moved in the z direction only.
   */
  void dualmove::direction(double x, double y, double z) {
    v.x=x;
    v.y=y;
    v.z=z;
  }

  /*!
   * Loads two macromolecules from disk and place them symmetrically
   * around the cell origin, separated by a specified distance.
   *
   * \param in Inputfile object that contains the protein filenames
   * \param g Macromolecular vector (will be erased!)
   * \param dist Initial distance between the protein mass-centers
   * \note The macromolecule vector is erased before any proteins are loaded!
   */
  void dualmove::load(inputfile &in, vector<macromolecule> &g, float dist) {
    if (dist==0)
      dist=rmax;
    g.clear();
    point a = v*(dist/2.);
    ioaam aam(*con);
    aam.load(*con, in, g); 
    g[0].move(*con, -( g[0].cm + a ));    // ...around the cell origo
    g[1].move(*con, -( g[1].cm - a ));    // ...along the z-axis
    g[0].accept(*con);
    g[1].accept(*con);
  }

  string dualmove::info() {
    std::ostringstream o;
    o <<  markovmove::info()
      << "#   Min/max separation  = " << rmin << " " << rmax << endl;
    return o.str();
  } 

  double dualmove::move(macromolecule &g1, macromolecule &g2) {
    du=0;
    if (dp==0)
      return du;
    cnt++;
    r=con->dist(g1.cm, g2.cm);
    point p;
    group g12=g1+g2;
    p.x=v.x*dp*slp.random_half();
    p.y=v.y*dp*slp.random_half();
    p.z=v.z*dp*slp.random_half();
    g1.move(*con, p);
    g2.move(*con,-p);
    double rtrial = con->dist(g1.cm_trial,g2.cm_trial);
    if (con->collision(g1.cm_trial)==true ||
        con->collision(g2.cm_trial)==true || 
        rtrial > rmax || rtrial < rmin) {
      rc=HC;
      dpsqr+=0;
      g1.undo(*con);
      g2.undo(*con);
      gofr.add(r);
      return du;
    }
    //#pragma omp parallel
    {
      //#pragma omp sections
      {
        //#pragma omp section
        { uold = pot->energy(con->p,g12) + pot->energy(con->p,g1,g2);   }
        //#pragma omp section
        { unew = pot->energy(con->trial,g12) + pot->energy(con->trial,g1,g2);   }
      }
    }
    du = unew-uold;
    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      r=con->dist(g1.cm, g2.cm);
      gofr.add(r);
      dpsqr+=4.*con->sqdist(g1.cm,g1.cm_trial);
      naccept++;
      g1.accept(*con);
      g2.accept(*con);
      return du;
    } else
      rc=ENERGY;
    du=0;
    dpsqr+=0;
    g1.undo(*con);
    g2.undo(*con);
    gofr.add(r);
    return du;
  }

  //---------- TRANSLATE GROUP ----------------
  class translate : public markovmove {
    public: 
      translate( ensemble&, container&, energybase&);
      double move(macromolecule &); 
  };

  translate::translate( ensemble &e,
      container &c, energybase &i ) : markovmove(e,c,i) {
    name = "MACROMOLECULAR TRANSLATION";
    runfraction=1.0;
    deltadp=1.;
    dp=10.;
  };

  double translate::move(macromolecule &g) {
    du=0;
    cnt++;
    point p;
    p.x=dp*slp.random_half();
    p.y=dp*slp.random_half();
    p.z=dp*slp.random_half();
    g.move(*con, p); 
    if (con->collision(g.cm_trial)==true) {
      rc=HC;
      dpsqr+=0;
      g.undo(*con);
      return du;
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
      dpsqr+=con->sqdist( g.cm, g.cm_trial );
      naccept++;
      g.accept(*con);
      return du;
    } else rc=ENERGY;
    du=0;
    dpsqr+=0;
    g.undo(*con);
    return du;
  }

  //-----------MOVE----------------------------------------
  /*! \brief Random move of a macromolecule in a container 
   *  \todo generalize, it is cubix periodic boundry specific
   *  \author Bjoern Persson
   */
  class move : public markovmove {
    private:
      box *b;
    public:
      move( ensemble&, container&, energybase &, box &, float);
      bool mOve(macromolecule &);
  };
  /*! \brief Class to prefom a random walk of a macromolecule in space
   *  \note Replaced by translate(?)
   */
  move::move(
      ensemble &e, container &c, energybase &i, box &BOX,float rf) : markovmove(e,c,i) {
    runfraction=rf;
    dp=100;
    deltadp=1;
    name="MACROMOLECULE MOVE";
    b=&BOX;
  }

  bool move::mOve(macromolecule &g) {
    du=0;
    if (slp.runtest(runfraction)==false)
      return false;
    cnt++;
    point rand;
    rand.x = dp*slp.random_half();
    rand.y = dp*slp.random_half();
    rand.z = dp*slp.random_half();
    g.move(*con, rand);
    for (int i=g.beg; i<(g.size()+g.beg); i++) { 
      b->boundary(con->trial[i]);
      if (con->collision( con->trial[i] )==true) 
        rc=HC;
    }
    if (rc==HC) {
      for (int i=g.beg; i<(g.size()+g.beg); i++) 
        con->trial[i] = con->p[i];
      std::cout << "rejected"<<endl;
      return false; }
    else {
      //#pragma omp parallel
      {
        //#pragma omp sections
        {
          //#pragma omp section
          { uold = pot->energy(con->p, g);   }
          //#pragma omp section
          { unew = pot->energy(con->trial,g);   }
        }
      }
      du += unew - uold; }
      if (ens->metropolis(du)==true) {
        rc=OK;
        utot+=du;
        naccept++;
        for (int i=g.beg; i<(g.beg+g.size()); i++)
          con->p[i] = con->trial[i];
        return true;
      } else rc=ENERGY;

      du=0;
      for (int i=g.beg; i<(g.size()+g.beg); i++) 
        con->trial[i] = con->p[i];
      return false;
  }


  //-------------- SALT MOVE ---------------------------------
  /*! \brief Move salt particles
   *  \author Mikael Lund
   */
  class saltmove : public markovmove {
    public:
      saltmove( ensemble &, container&, energybase& );
      double move(group &, int);  //!< Move a single particle
      double move(group &);         //!< Loop over group particles (randomly)
  };

  saltmove::saltmove(
      ensemble &e, container &c, energybase &i ) : markovmove(e,c,i) {
    dp=30;
    deltadp=2;
    name="SALT DISPLACEMENTS";
    runfraction=0.5;
  }

  /*! \param group Group containing mobile ions
  */
  double saltmove::move(group &g) {
    du=0;
    if (slp.runtest( runfraction )==false)
      return du;
    double sum=0;
    for (unsigned short i=0; i<g.size(); i++) {
      move(g, 1);
      sum+=du;
    }
    du=sum;
    return du;
  }
  double saltmove::move(group &g, int n) {
    du=0;
    cnt++;
    n=g.displace(*con, dp); 
    if (con->collision( con->trial[n] )==true)
      rc=HC;
    else {
      //#pragma omp parallel
      {
        //#pragma omp sections
        {
          //#pragma omp section
          { uold = pot->energy(con->p, n);   }
          //#pragma omp section
          { unew = pot->energy(con->trial,n);   }
        }
      }
      du = unew - uold;
      if (ens->metropolis(du)==true) {
        rc=OK;
        utot+=du;                               // track energy changes
        dpsqr+=con->sqdist(con->p[n],con->trial[n]); // track avg. displacement
        naccept++;                              // accept counter
        con->p[n] = con->trial[n];              // Accept move
        return du;
      } else rc=ENERGY;
    }
    du=0;
    dpsqr+=0;
    con->trial[n] = con->p[n];
    return du;
  }
}
#endif
