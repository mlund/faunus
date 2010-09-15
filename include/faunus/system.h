#ifndef FAU_SYSTEM_H
#define FAU_SYSTEM_H

namespace Faunus {

  /*!
   * \brief Base class for entire simulation systems
   * \author Mikael Lund
   * \date Malmo, 2010
   */
  class simBottle {
    protected:
      iopqr pqr;
      ioaam aam;
      inputfile in;
      systemenergy sys;
      mcloop loop;

    public:
      container* cPtr;                  //!< Pointer to the container for this system
      string name;                      //!< Name of the system (arbitrary)
      string prefix;                    //!< Prefix for input files, output etc.
      virtual void prepare()=0;         //!< Prepare and setup systems
      virtual void microloop()=0;       //!< Things to do in each micro move
      virtual void macroloop()=0;       //!< Things to do in each macro move
      virtual void save()=0;            //!< Save data to disk
      virtual double systemEnergy()=0;  //!< Calculate system energy
      virtual string preinfo()=0;       //!< Give initial information string
      virtual string postinfo()=0;      //!< Give final information string
      simBottle(string pfx) : in(pfx+".conf"), loop(in), sys(0.) {
        prefix=pfx;
      }
  };

  /*!
   * \brief Simulation system for many molecules in the Npt ensemble
   * \author Mikael Lund
   * \date Lund 2010
   */
  template<typename Tcon, typename Tpot> class NpTmolecular : public simBottle {
    private:
      string dumpfile;  //!< Filename of container dumpfile
    protected:
      canonical nvt;
      vector<macromolecule> g;
      io fileio;

      macrorot *mmRot;
      translate *mmTrans;
      transrot *mmTransrot;
      isobaric<Tpot> *mmVol;

    public:
      Tcon con; //!< Container
      Tpot pot; //!< Interaction scheme (based on energybase)

      bool boolTranslate, boolVolume, boolRotate, boolCluster;

      NpTmolecular(string pfx) : simBottle(pfx), con(in), pot(in) {
        cPtr=&con;
        mmVol = new isobaric<Tpot>(nvt,con,pot,in);
        mmRot = new macrorot(nvt,con,pot);
        mmTrans = new translate(nvt,con,pot,in);
        mmTransrot = new transrot(nvt,con,pot);
        dumpfile=prefix+"confout.dump";
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

        sys.update( systemEnergy() ); 
      }

      double systemEnergy() {
        int n=g.size();
        double u=0;
        for (int i=0; i<n-1; i++)
          for (int j=i+1; j<n; j++)
            u+=pot.energy(con.p, g[i], g[j]);
        return u;
      }

      void microloop() {
        // Combined translation and rotation N times
        int N=g.size();
        for (int i=0; i<N; i++) {
          mmTransrot->dpt = 0.008 * pow(con.len,2);
          sys+=mmTransrot->move( g, N*slp.random_one() );
        }
      }

      void macroloop() {}

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
      }
  };

  /*!
   * \brief NVT simulation system for atomic species
   * \author Mikael Lund
   * \date Lund 2010
   */
  template<typename Tcon, typename Tpot> class NVTatomic : public simBottle {
    private:
      string dumpfile;  //!< Filename of container dumpfile
    protected:
      canonical nvt;
      io fileio;
      saltmove *sm;
    public:
      Tcon con; //!< Container
      Tpot pot; //!< Interaction scheme (based on energybase)

      NVTatomic(string pfx) : simBottle(pfx), con(in), pot(in) {
        cPtr=&con;
        dumpfile=prefix+"confout.dump";
      }

      void prepare() {}
      double systemEnergy() {}
      string preinfo() {}
      string postinfo() {}
      void microloop() {}
      void macroloop() {}
      void save() {}
  };

} // end of namespace
#endif
