#ifndef FAU_NPT_MOLECULAR_H
#define FAU_NPT_MOLECULAR_H

#include "faunus/bottles/base.h"

namespace Faunus {

  /*!
   * \brief Simulation system for many molecules in the NpT ensemble
   * \author Mikael Lund
   * \date Lund 2010
   */
  template<typename Tcon, typename Tpot> class npt_molecular : public bottle {
    private:
      string dumpfile;  //!< Filename of container dumpfile
      unsigned int microcnt, macrocnt;
      ioxtc xtc;
    protected:
      canonical nvt;
      xyfile boxtrj, energytrj;
    
      macrorot *mmRot;
      translate *mmTrans;
      transrot *mmTransrot;
      isobaric<Tpot> *mmVol;

    public:
      vector<macromolecule> g;
      Tcon con; //!< Container
      Tpot pot; //!< Interaction scheme (based on energybase)

      bool boolTranslate, boolVolume, boolRotate, boolCluster;
      bool movie; //!< Set to true if an xtc tracjectory file should be saved

      npt_molecular(string pfx) : bottle(pfx), con(in), pot(in),
                    boxtrj(prefix+".boxlen.dat"), energytrj(prefix+".energy.dat"), xtc(con.len) {
        cPtr=&con;
        pPtr=&pot;
        microcnt=0;
        macrocnt=0;
        movie=in.getboo("movie", false);

        // load molecules
        dumpfile=prefix+".dump";
        if (in.getboo("lattice")==true)
          aam.loadlattice(con, in, g);
        else
          aam.load(con, in, g);
        if ( con.loadFromDisk(dumpfile) ) {
          for (int i=0; i<g.size(); ++i)
            g[i].masscenter(con);
          fout << "# Configuration read from: " << dumpfile << endl;
        } 

        // prepare move routines
        mmVol = new isobaric<Tpot>(nvt,con,pot,in);
        mmRot = new macrorot(nvt,con,pot);
        mmTrans = new translate(nvt,con,pot,in);
        mmTransrot = new transrot(nvt,con,pot);
        P=mmVol->P;

        usys.initialize( systemEnergy() );

        fout << in.info() << con.info() << pot.info();
        fout << "# Macromolecular charges and dipole moments:" << endl;
        fout << "#   First: z=" <<  g[0].charge(con.p) << " mu=" << g[0].dipole(con) << endl;
        fout << std::flush;
      }
    
      ~npt_molecular() {
        delete mmRot;
        delete mmTrans;
        delete mmTransrot;
        delete mmVol;
        boxtrj.close();
        energytrj.close();
      }

      double systemEnergy() {
        int n=g.size();
        for (int i=0; i<n; i++)
          g[i].masscenter(con);
        double u=0;
//#pragma omp parallel for reduction (+:u) schedule (dynamic)
        for (int i=0; i<n-1; i++)
          for (int j=i+1; j<n; j++)
            u+=pot.energy(con.p, g[i], g[j]);
        return u;
      }

      void microloop(int ntimes=1) {
        int N=g.size();
        for (int i=0; i<N; ++i)
          g[i].masscenter(con);
        while ((ntimes--)>0) {
          microcnt++;
          switch (slp.rand() % 2) {
            case 0: // combined translation and rotation N times
              mmTransrot->dpt = 0.0008 * pow(con.len,2);
              for (int i=0; i<N; ++i)
                usys+=mmTransrot->move( g, N*slp.random_one() );
              break;
            case 1: // volume move
              usys+=mmVol->move(g);
              break;
          }
          if (slp.random_one()>0.995) {
            boxtrj.add(microcnt, con.len);
            energytrj.add(microcnt, usys.sum);
            if (movie==true)
              xtc.save(prefix+".xtc", con);
          }
        }
      }

      void macroloop() {
        macrocnt++;
        usys.update( systemEnergy() );
        save();
        fout << "# Macroloop " << macrocnt << " completed." << endl
             << "#   Energy drift = " << usys.drift() << " kT" << endl
             << "#   Acceptances:   vol=" << mmVol->accepted()*100 << " transrot=" << mmTransrot->accepted()*100 << endl
             << std::flush;
      }

      void save() {
        pqr.save(prefix+".pqr", con.p);
        con.saveToDisk(dumpfile);
        mmVol->Ldist.write(prefix+".lendist.dat");
      }

      void finish() {
        fout << usys.info() << mmTransrot->info() << mmVol->info();
      }
  };

} // end of namespace
#endif
