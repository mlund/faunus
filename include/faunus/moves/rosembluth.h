#ifndef FAU_ROSEMBLUTH_H
#define FAU_ROSEMBLUTH_H
#include "faunus/moves/base.h"
namespace Faunus {
  /*! \brief Rotate group around its mass-center.
   *  \author Mikael Lund
   *  \date Prague 2007
   *  \warning ....
   *  \todo a lot...
   */
  class rosembluth : public markovmove { 
    private:
      short index;
      string bondtype;
      unsigned int k;
      unsigned int gpolymer; //!< GC group location of polymer
      unsigned int gcounter; //!< GC group location of counter ions
      double mu;
      vector<string> polymer, counter;
      grandcanonical *gcPtr;
      void insert();
      void remove();
    public:
      rosembluth( grandcanonical&, container&, energybase&, inputfile &, int);
      double move();
      string info();
  };
  
  rosembluth::rosembluth( grandcanonical &gc,
      container &c, energybase &i, inputfile &in, int idx) : markovmove(gc,c,i)
  {
    name = "ROSEMBLUTH INSERTION";
    runfraction=1.0;
    deltadp=0;
    dp=0;
    index=idx;
    gcPtr=&gc;

    std::ostringstream rb_counter, rb_k, rb_polymer, rb_bond, rb_mu;
    rb_counter << "RB" << index << "_counterions";
    rb_polymer << "RB" << index << "_polymer";
    rb_k << "RB" << index << "_ktrials";
    rb_bond << "RB" << index << "_bond";
    rb_mu << "RB" << index << "_mu";

    mu = in.getflt( rb_mu.str(), 1e6 );
    if (mu==1e6) {
      runfraction=0;
      return;
    }

    k        = in.getint( rb_k.str(), 1 );
    bondtype = in.getstr( rb_bond.str(), "none");
    polymer  = in.getvec( rb_polymer.str(), "NA" );
    gpolymer = gcPtr->findgroup( polymer[0] );
    counter  = in.getvec( rb_counter.str(), "CL" );
    gcounter = gcPtr->findgroup( counter[0] );
  }

  void rosembluth::insert() {
    // just an example:
    // insert a monomer at end of monomer group
    gcPtr->insert(con->trial, gcPtr->g[gpolymer].end, con->atom(polymer[0]), polymer.size() );
  }

  void rosembluth::remove() {
    // just an example:
    int m=gcPtr->g[gpolymer].random(); // Pick random monomer
    int c=gcPtr->g[gcounter].random(); // Pick random counterion

    gcPtr->erase( con->trial, m, polymer.size() );
  }

  /*!
   * \todo Cell overlap test missing
   */
  double rosembluth::move() {
    du=0;
    cnt++;
    // 1. Randomly insert or remove.
    // 2. calc. initial energy ("p" vector).
    // 3. calc. final energy ("trial" vector).
    {
      {
        //{ uold = pot->energy(con->p, g);   }
        //{ unew = pot->energy(con->trial, g);   }
      }
    }
    du = unew-uold;
    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      naccept++;
      //accept move!
      return du;
    } else rc=ENERGY;
    du=0;
    //undo move!
    return du;
  }

  string rosembluth::info() {
    std::ostringstream o;
    if (runfraction>0) {
      o << markovmove::info()
        << "#   Index                     = " << index << endl
        << "#   Chemical potential (kT)   = " << mu << endl
        << "#   No. of monomers           = " << polymer.size() << " " << gcPtr->g[gpolymer].size() << endl
        << "#   No. of counter ions       = " << counter.size() << " " << gcPtr->g[gcounter].size() << endl
        << "#   Number of trials          = " << k << endl
        << "#   Bond type                 = " << bondtype << endl;
    }
    return o.str();
  }
} // namespace
#endif
