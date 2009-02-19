#ifndef FAU_ROSEMBLUTH_H
#define FAU_ROSEMBLUTH_H
#include "faunus/moves/base.h"
namespace Faunus {
  /*! \brief Rosembluth polymer and salt insertion
   *  \author Christophe Labbez and Mikael Lund
   *  \date Dijon / Lund, 2009
   *  \warning ....
   *  \todo a lot...
   */
  class rosembluth : public markovmove { 
    private:
      grandcanonical *gcPtr;    // Pointer to GC class
      struct data {
        group* gPtr;            // Pointer to group containing this species
        unsigned short valency; // Stoechiometry
        short charge;           // Particle charge
        vector<string> seq;     // Sequence of particle names
      };
      data polymer, counter;    // Store data for "polymer" and counter ions
      short index;              // Rosembluth index
      string bondtype;          // Insertion scheme or bond type
      unsigned int k;           // Rosembluth k-value
      double mu;                // Chemical potential
      void insert();
      void remove();
    public:
      // Keep as much as possible private! Easier for the user of the class.
      rosembluth( grandcanonical&, container&, energybase&, inputfile &, int);
      double move();
      string info();
  };
  
  rosembluth::rosembluth( grandcanonical &gc,
      container &c, energybase &i, inputfile &in, int idx) : markovmove(gc,c,i)
  {
    name = "ROSEMBLUTH INSERTION";
    cite = "Rosembluth reference...or Frenkel?";
    runfraction=1.0;
    deltadp=0;
    dp=0;
    index=idx;
    gcPtr=&gc;

    // Fetch from input
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

    // Set parameters
    k               = in.getint( rb_k.str(), 1 );
    bondtype        = in.getstr( rb_bond.str(), "none");
    polymer.seq     = in.getvec( rb_polymer.str(), "NA" );
    counter.seq     = in.getvec( rb_counter.str(), "CL" );
    polymer.valency = polymer.seq.size();
    counter.valency = counter.seq.size();
    polymer.gPtr    = & gcPtr->g[  gcPtr->findgroup( polymer.seq[0] ) ];
    counter.gPtr    = & gcPtr->g[  gcPtr->findgroup( counter.seq[0] ) ];
  }

  void rosembluth::insert() {
    // just an example:
    // insert a monomer at end of monomer group
    gcPtr->insert(con->trial, polymer.gPtr->end, con->atom(polymer.seq[0]), polymer.valency );
  }

  void rosembluth::remove() {
    // just an example:
    int m=polymer.gPtr->random(); // Pick random monomer
    int c=counter.gPtr->random(); // Pick random counterion

    gcPtr->erase( con->trial, m, polymer.valency );
    gcPtr->erase( con->trial, c, counter.valency );
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
        << "#   No. of monomers           = " << polymer.valency << " " << polymer.gPtr->size() << endl
        << "#   No. of counter ions       = " << counter.valency << " " << counter.gPtr->size() << endl
        << "#   Number of trials          = " << k << endl
        << "#   Bond type                 = " << bondtype << endl;
    }
    return o.str();
  }
} // namespace
#endif
