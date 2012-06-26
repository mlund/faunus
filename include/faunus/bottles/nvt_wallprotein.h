#ifndef FAU_NVT_WALLPROTEIN_H
#define FAU_NVT_WALLPROTEIN_H

#include "faunus/bottles/base.h"

namespace Faunus {

  template<typename Tcon, typename Tpot> class nvt_wallprotein : public bottle {
    protected:
      canonical nvt;
      xyfile boxtrj, energytrj;
    
      macrorot *mmRot;
      translate *mmTrans;
      eqtitrate *mmTit;

    private:
      string dumpfile;  //!< Filename of container dumpfile
      unsigned int microcnt, macrocnt;
      ioxtc xtc;

    public:
      polymer pol;
      Tcon con; //!< Container
      Tpot pot; //!< Interaction scheme (based on energybase)

      bool boolTranslate, boolVolume, boolRotate, boolCluster;
      bool movie; //!< Set to true if an xtc tracjectory file should be saved

      nvt_wallprotein(string pfx) :
        bottle(pfx),
        boxtrj(prefix+".boxlen.dat"),
        energytrj(prefix+".energy.dat"),
        con(in),
        xtc(1000.),
        pot(in) {
          cPtr=&con;
          pPtr=&pot;
          microcnt=0;
          macrocnt=0;
          movie=in.getboo("movie", false);

          // load molecule
          pol.babeladd( con, in );          //  add from input
          con.trial=con.p;                  //  synchronize particle vector                          
          pol.masscenter(con);              //  update masscenter                                    
          pol.move(con, -pol.cm+(con.slice_min+con.slice_max)*0.5);  // translate to the middle of the slice or the origo (0,0,0) if no slice defined
          pol.accept(con);                  //  accept translation                                   
          cout << pol.info();           

          pot.expot.update(con);

          // prepare move routines
          mmTit = new eqtitrate(nvt,con,pot,in,"eqtit_");
          mmRot = new macrorot(nvt,con,pot);
          mmTrans = new translate(nvt,con,pot,in);

          mmTrans->dpv.x=mmTrans->dpv.y=0; // move molecule only in z direction
          mmTrans->dp = 10.;

          usys.initialize( systemEnergy() );

          fout << in.info() << con.info() << pot.info();
          fout << std::flush;
        }

      ~nvt_wallprotein() {
        delete mmTit;
        delete mmRot;
        delete mmTrans;
        boxtrj.close();
        energytrj.close();
      }

      double systemEnergy() {
        pol.masscenter(con);
        return pot.energy(con.p,pol) + mmTit->intrinsicenergy(con.p) + pot.uself_polymer(con.p, pol);  
      }

      void microloop(int ntimes=1) {
        pol.masscenter(con);
        while ((ntimes--)>0) {
          microcnt++;
          switch (slp.rand() % 2) {
            case 0: // translation
              usys+=mmTrans->move( pol );
              break;
            case 1:
              usys+=mmTit->move();
              pol.charge(con.p);
              break;
          }
          if (slp.random_one()>0.999) {
            boxtrj.add(microcnt, con.len.x);
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
          << "#   Acceptances: trans=" << mmTrans->accepted()*100 << endl
          << std::flush;
      }

      void save() {
        pqr.save(prefix+".pqr", con.p);
        con.saveToDisk(dumpfile);
      }

      void finish() {
        fout << usys.info() << mmTrans->info() << mmTit->info();
      }
  };

} // end of namespace
#endif
