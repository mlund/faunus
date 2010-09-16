#ifndef FAU_NPT_MOLECULAR_H
#define FAU_NPT_MOLECULAR_H

#include "faunus/bottles/base.h"

namespace Faunus {

  /*!
   * \brief Simulation system for many molecules in the Npt ensemble
   * \author Mikael Lund
   * \date Lund 2010
   */
  template<typename Tcon, typename Tpot> class npt_molecular : public simBottle {
    private:
      string dumpfile;  //!< Filename of container dumpfile
    protected:
      canonical nvt;
      io fileio;

      macrorot *mmRot;
      translate *mmTrans;
      transrot *mmTransrot;
      isobaric<Tpot> *mmVol;

    public:
      vector<macromolecule> g;
      Tcon con; //!< Container
      Tpot pot; //!< Interaction scheme (based on energybase)

      bool boolTranslate, boolVolume, boolRotate, boolCluster;

      npt_molecular(string pfx) : simBottle(pfx), con(in), pot(in) {
        cPtr=&con;
        mmVol = new isobaric<Tpot>(nvt,con,pot,in);
        mmRot = new macrorot(nvt,con,pot);
        mmTrans = new translate(nvt,con,pot,in);
        mmTransrot = new transrot(nvt,con,pot);
        dumpfile=prefix+".dump";
      }

      void prepare() {
        // Load macromolecules
        if (in.getboo("lattice")==true)
          aam.loadlattice(con, in, g);
        else
          aam.load(con, in, g);
        if ( con.loadFromDisk(dumpfile) ) {
          for (int i=0; i<g.size(); i++)
            g[i].masscenter(con);
          cout << "# Initial configuration read from " << dumpfile << endl;
        }  

        usys.initialize( systemEnergy() );
      }

      double systemEnergy() {
        int n=g.size();
        double u=0;
        for (int i=0; i<n-1; i++)
          for (int j=i+1; j<n; j++)
            u+=pot.energy(con.trial, g[i], g[j]);
        return u;
      }

      void microloop() {
        // Combined translation and rotation N times
        int N=g.size();
        for (int i=0; i<N; i++) {
          mmTransrot->dpt = 0.008 * pow(con.len,2);
          usys+=mmTransrot->move( g, N*slp.random_one() );
        }
      }

      void macroloop() {
        usys.update( systemEnergy() );
      }

      void save() {
        pqr.save(prefix+"confout.pqr", con.p);
        con.saveToDisk(dumpfile);     
      }

      string preinfo() {
        std::ostringstream o;
        o << in.info() << con.info() << pot.info();
        return o.str();
      }

      string postinfo() {
        std::ostringstream o;
        o << usys.info() << mmTransrot->info();
        return o.str();
      }
  };

} // end of namespace
#endif
