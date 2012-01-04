#ifndef FAU_NVT_TWOBODY_H
#define FAU_NVT_TWOBODY_H

#include "faunus/bottles/base.h"

namespace Faunus {

  /*!
   * \brief Bottle for two rigid macromolecules w. explicit salt in NVT ensemble
   * \author Mikael Lund
   * \date Lund 2011
   */
  template<typename Tcon, typename Tpot> class nvt_twobody : public bottle {
    private:
      canonical nvt;
      saltmove *sm;
      dualmove *dm;
      macrorot *mr;

      string dumpfile;  //!< Filename of container dumpfile
      unsigned int microcnt, macrocnt;
      ioxtc xtc;
      salt saltgroup;
      vector<macromolecule> g;

    public:
      Tcon con; //!< Container
      Tpot pot; //!< Interaction scheme (based on energybase)

      bool movie; //!< Set to true if an xtc tracjectory file should be saved

      nvt_twobody(string pfx) : bottle(pfx), xtc(1000.), con(in), pot(in) {
        cPtr=&con;
        pPtr=&pot;
        microcnt=0;
        macrocnt=0;
        movie=in.getboo("movie", false);
        dumpfile=prefix+".dump";

        mr = new macrorot(nvt,con,pot);
        dm = new dualmove(nvt,con,pot);
        dm->setup(in);
        dm->load(in, g, 80.);
        saltgroup.add(con,in);
        sm = new saltmove(nvt,con,pot,in);

        aam.load(con,prefix+".aam");                                                   
        g[0].masscenter(con);           
        g[1].masscenter(con); 

        usys.initialize( systemEnergy() );

        fout << in.info() << con.info() << pot.info() << std::flush;
      }

      ~nvt_twobody() {
        delete sm;
        delete dm;
        delete mr;
      }

      double systemEnergy() {
        for (int i=0; i<g.size(); i++)
          g[i].masscenter(con);
        double u=0;
        u += pot.energy(con.p, g[0], g[1]) +
          pot.energy(con.p, saltgroup) +
          pot.internalElectrostatic(con.p, g[0]) +
          pot.internalElectrostatic(con.p, g[1]) +
          pot.internal(con.p, saltgroup);                // System energy analysis
        return u;
      }

      void microloop(int ntimes=1) {
        short i,n;                                                                  
        switch (rand() % 3) { 
          case 0:
            usys+=sm->move(saltgroup); 
            break;                     
          case 1:                      
            for (n=0; n<2; n++) {      
              i = rand() % g.size();   
              usys+=mr->move(g[i]);     
            }                          
            break;                     
          case 2:                      
            usys+=dm->move(g[0], g[1]); 
            break;
        }
        if (slp.random_one()>0.995) {
          //if (movie==true)
          //  xtc.save(prefix+".xtc", con);
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
        pqr.save(prefix+".pqr", con.p);
        aam.save(prefix+".aam", con.p);
        dm->gofr.write(prefix+".rdf"); 
        con.saveToDisk(dumpfile);
      }

      void finish() {
        fout << usys.info() << sm->info() << dm->info() << mr->info();
      }
  };

} // end of namespace
#endif
