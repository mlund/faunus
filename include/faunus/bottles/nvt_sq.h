#ifndef FAU_NVT_SQ_H
#define FAU_NVT_SQ_H

#include "faunus/bottles/base.h"
#include "faunus/moves/reptation.h"
#include "faunus/moves/switch.h"

namespace Faunus {

  /*!
   * \brief Simulation system for many molecules in the Npt ensemble
   * \author Mikael Lund
   * \date Lund 2010
   */
  template<typename Tcon, typename Tpot> class nvt_sq : public bottle {
    private:
      unsigned int microcnt, macrocnt;
      string dumpfile;  //!< Filename of container dumpfile
    protected:
      canonical nvt;
      io fio;
//      iopqr pqr;
//      ioaam aam;
      //ioxtc xtc; 

      macrorot     *mrPol;
      translate    *mtPol;
      crankShaft   *csPol;
      reptation    *repPol;
      switch_type  *swPol;
      clustertrans *cltPol;

      aggregation  *agg;

      double randy, picky, rot_f, crank_f, rep_f, trans_f, tit_f, clt_f, sum;

    public:
      vector<macromolecule> mpolymers;
      vector<polymer>        polymers;
      Tcon con; //!< Container
      Tpot pot; //!< Interaction scheme (based on energybase)
      

      bool boolTranslate, boolVolume, boolRotate, boolCluster;

      nvt_sq(string pfx) : bottle(pfx), con(in), pot(in) {
        cPtr=&con;
        pPtr=&pot;
        mrPol  = new macrorot(nvt, con, pot); 
        mtPol  = new translate(nvt, con, pot, in);
        csPol  = new crankShaft(nvt, con, pot, in);
        repPol = new reptation(nvt, con, pot, in);
        swPol  = new switch_type(nvt, con , pot, in);
        cltPol = new clustertrans(nvt, con, pot, mpolymers);
        agg    = new aggregation(con, polymers, in.getflt("cluster_def"));

        //xtc    = new ioxtc(con.len);
        swPol->full_ham=true;
        point uv;                              //'shadow' vect
#ifdef BABEL
        for(int i=0;i<in.getint("N_polymer"); i++){
          polymer pol;
          pol.babeladd( con, in );
          uv.ranunit(slp);
          pol.move(con, -pol.cm);              // Translate polymer to origo (0,0,0)
          pol.move(con, uv*con.len);           // .. and a random position
          pol.accept(con);                     // .. accept translation
          polymers.push_back(pol);
          mpolymers.push_back(pol);
        }
#endif

        // Load stored configuration?
        if(  aam.load(con, prefix+"confout.aam"))
          cout << "# Previus configuration loaded"<<endl;

        //Sloppy sync of radii and mw... blame it on openbabel
        for (int i=0; i<con.p.size(); i++){
          con.p[i].radius=atom[con.p[i].id].radius;
          con.p[i].mw=atom[con.p[i].id].mw;
          con.trial[i].radius=atom[con.p[i].id].radius;
          con.trial[i].mw=atom[con.p[i].id].mw;
        }
        usys.initialize( systemEnergy() );
       
        pqr.save(prefix+"confout.pqr",con.p);
       
        cout << con.info() << atom.info()
             << pot.info()
             << in.info();
       
        int N=in.getint("N_polymer");
       
        //Parameters for input control of Markov chain
        rot_f=in.getflt("rot_f"), crank_f=in.getflt("crank_f"), rep_f=in.getflt("rep_f"), clt_f=in.getflt("clt_f");
        trans_f=in.getflt("trans_f"), tit_f=in.getflt("tit_f"),sum=crank_f+rep_f+rot_f+trans_f+tit_f;
        rot_f/=sum, crank_f/=sum, rep_f/=sum, trans_f/=sum, tit_f/=sum, clt_f/=sum;

      }
      ~nvt_sq() {
        delete mrPol;
        delete mtPol;
        delete repPol;
        delete swPol;
        delete cltPol;
        //delete xtc;
        delete agg;
      }

      double systemEnergy() {
        return pot.energy(con.p);
      }

      void microloop(int ntimes=1) {
        int N=polymers.size();
        for (int i=0; i<N; ++i)
          polymers[i].masscenter(con);
        while ((ntimes--)>0) {
          microcnt++;
          randy=slp.random_one();
          picky=slp.random_one();
          if(trans_f>randy){
            usys+=mtPol->move(polymers[int(picky*N)]);         // translate polymers
          }
          if(trans_f+rot_f>randy && trans_f<randy){
            polymers[int(picky*N)].masscenter(con);
            usys+=mrPol->move(polymers[int(picky*N)]);         // rotate polymers
          }
          if(trans_f+rot_f+crank_f>randy && trans_f+rot_f<randy){
            usys+=csPol->move(polymers[int(picky*N)]);         // crankshaft
          }
          if(trans_f+rot_f+crank_f+rep_f>randy && trans_f+rot_f+crank_f<randy){
            usys+=repPol->move(polymers[int(picky*N)]);        // reptation
          }
          if(1.0-tit_f-clt_f<randy && 1.0-clt_f>randy){
            usys+=swPol->titrateall();                        // titrate titratable sites
          }
          if(1.0-clt_f<randy){
            usys+=cltPol->move(mpolymers);                     // cluster translation of 'all' polymers
          }
          if ( slp.random_one()>0.99 ) {
            // ANALYSIS
            agg->count();
          }
        }
      }

      void macroloop() {
        macrocnt++;
        usys.update( systemEnergy() );
        save();
        fout << "# Macroloop " << macrocnt << " completed." << endl
             << "#   Energy drift = " << usys.drift() << " kT" << endl
             << std::flush;

      }

      void save() {
        std::ostringstream o;
        for (int i=0; i<polymers.size(); i++)
          o << polymers[i].getVMDBondScript() <<endl;
        fio.writefile(prefix+"vmdbonds.tcl", o.str());
        pqr.save(prefix+"confout.pqr",con.p);
        aam.save(prefix+"confout.aam", con.p);
        agg->write(prefix+"aggregation.dat");
      }

      void finish() {
        fout << usys.info() << mtPol->info() << mrPol->info()
             << csPol->info() << swPol->info() << cltPol->info();
      }

  };

} // end of namespace
#endif
