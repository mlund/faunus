#include <faunus/moves/rosenbluth.h>
namespace Faunus {
  // Constructor
  rosenbluth::rosenbluth( grandcanonical &gc,
      container &c, energybase &i, inputfile &in, int idx) : markovmove(gc,c,i)
  {
    name = "ROSENBLUTH INSERTION";
    cite = "Rosenbluth reference...or Frenkel?";
    runfraction=1.0;
    deltadp=0;
    dp=0;
    index=idx;
    gcPtr=&gc;
    ins.resize(0);

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

  void rosenbluth::insert() {
    // just an example:
    // insert a monomer at end of monomer group
    int i = polymer.gPtr->end; // last particle in group
    ins.push_back(i);          // keep track of what was moved.
    gcPtr->insert(con->trial, polymer.gPtr->end, con->atom(polymer.seq[0]), polymer.valency );

    // let's assign a random position (within the container) to the inserted particle:
    con->randompos( con->trial[i] );
    // ...or manually:
    con->trial[i].x = 2.3;
    // ...hmm maybe this took us outside the container boundaries. Let's fix that:
    con->boundary( con->trial[i] );
  }

  void rosenbluth::remove() {
    // just an example:
    int m=polymer.gPtr->random(); // Pick random monomer
    int c=counter.gPtr->random(); // Pick random counterion

    gcPtr->erase( con->trial, m, polymer.valency );
    gcPtr->erase( con->trial, c, counter.valency );

    ins.push_back(m);
    ins.push_back(c);
  }

  /*!
   * \todo Cell overlap test missing
   */
  double rosenbluth::move() {
    ins.clear();
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
      con->p = con->trial; // sync. particle vectors
      rc=OK;
      utot+=du;
      naccept++;
      //accept move!
      return du;
    } else rc=ENERGY;
    du=0;
    //undo move!
    //1) insert or remove particle(s)
    con->trial = con->p; // sync. particle vectors
    return du;
  }

  string rosenbluth::info() {
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
}//namespace
