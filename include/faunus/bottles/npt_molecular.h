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
      io fio;
    
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
        pPtr=&pot;
        mmVol = new isobaric<Tpot>(nvt,con,pot,in);
        mmRot = new macrorot(nvt,con,pot);
        mmTrans = new translate(nvt,con,pot,in);
        mmTransrot = new transrot(nvt,con,pot);
        dumpfile=prefix+".dump";
        P=mmVol->P;
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
        for (int i=0; i<n; i++)
          g[i].masscenter(con);
        double u=0;
        for (int i=0; i<n-1; i++)
          for (int j=i+1; j<n; j++)
            u+=pot.energy(con.p, g[i], g[j]);
        return u;
      }

      void microloop(int ntimes) {
        while ((ntimes--)>0) {
          int N=g.size();
          for (int i=0; i<N; i++)
            g[i].masscenter(con);

          switch (rand() % 2) {
            case 0: // combined translation and rotation N times
              for (int i=0; i<N; i++) {
                mmTransrot->dpt = 0.008 * pow(con.len,2);
                usys+=mmTransrot->move( g, N*slp.random_one() );
              }
              break;
            case 1: // volume move
              usys+=mmVol->move(g);
              break;
          }
        }
      }

      void macroloop() {
        usys.update( systemEnergy() );
      }

      void save() {
        pqr.save(prefix+".pqr", con.p);
        con.saveToDisk(dumpfile);
        mmVol->Ldist.write(prefix+".lendist.dat");
      }

      string preinfo() {
        std::ostringstream o;
        o << in.info() << con.info() << pot.info();
        return o.str();
      }

      string postinfo() {
        std::ostringstream o;
        o << usys.info() << mmTransrot->info() << mmVol->info();
        return o.str();
      }
  };

} // end of namespace
#endif
