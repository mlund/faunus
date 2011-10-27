#include <faunus/widom.h>
#include <faunus/species.h>
#include <faunus/energy.h>
#include <faunus/space.h>
#include <faunus/inputfile.h>
#include <faunus/geometry.h>
#include <faunus/textio.h>

namespace Faunus {

  namespace Analysis {
    double PolymerShape::gyrationRadiusSquared(const Group &pol, const Space &spc) {
      assert( spc.geo->dist(pol.cm, pol.massCenter(spc))<1e-9 && "Mass center must be in sync.");
      double x=0, r2=0, sum=0;
      Point t, o;
      for (int i=pol.beg; i<=pol.end; i++) {
        t = spc.p[i]-pol.cm;                // vector to center of mass
        spc.geo->boundary(t);               // periodic boundary (if any)
        r2 = spc.geo->sqdist(t,o);          // squared distance to cm
        x += r2 * spc.p[i].mw;              // m*r^2
        sum += spc.p[i].mw;                 // total mass
      }
      assert(sum>0 && "Zero molecular weight not allowed.");
      return x*(1./sum);
    }

    void PolymerShape::sample(const Group &pol, const Space &spc) {
      assert( pol.beg!=pol.end && "Polymer must have at least two particles.");
      double r2=gyrationRadiusSquared(pol,spc);
      Rg2[pol.name]+=r2;
      Rg[pol.name]+=sqrt(r2);
      Re2[pol.name] += spc.geo->sqdist( spc.p[pol.beg], spc.p[pol.end] ); //end-2-end squared
    }

    string PolymerShape::info() {
      char w=13;
      using namespace textio;
      std::ostringstream o;
      o << header("Polymer Shape");
      o << indent(SUBSUB) << std::left << setw(w) << "Polymer"
        << setw(w+5) << bracket("Rg"+squared)
        << setw(w+5) << bracket("Rg")+squared
        << setw(w+7) << rootof+bracket("Re"+squared)
        << setw(w) << bracket("Re"+squared)+"/"+bracket("Rg"+squared) << endl;
      for (auto &m : Rg2)
        o << indent(SUBSUB) << std::left << setw(w) << m.first << setw(w) << m.second.avg()
          << setw(w) << pow( Rg[m.first].avg(),2 )
          << setw(w) << sqrt(Re2[m.first].avg())
          << setw(w) << Re2[m.first].avg() / Rg2[m.first].avg() << endl;
      return o.str();
    }

    Widom::Widom(Space &spc, Energy::Energybase &pot) {
      spcPtr=&spc;
      potPtr=&pot;
    }

    void Widom::sample(int ghostin) {
      assert(spcPtr->geo!=NULL);
      int n=g.size();
      for (int k=0; k<ghostin; k++) {     // insert ghostin times
        double du=0;
        for (int i=0; i<n; i++)
          spcPtr->geo->randompos( g[i] ); // random ghost positions
        for (int i=0; i<n; i++)
          du+=potPtr->all2p( spcPtr->p, g[i] );    // energy with all particles in space
        for (int i=0; i<n-1; i++)
          for (int j=i+1; j<n; j++)
            du+=potPtr->p2p( g[i], g[j] );   // energy between ghost particles
        expsum += exp(-du);
      }
    }

    void Widom::addGhost(particle p) {
      g.push_back(p);
    }

    void Widom::addGhost(Space &c) {
      std::map<short,bool> map;
      for (auto p : c.p)
        map[ p.id ] = true;
      for (auto &m : map) {
        particle a;
        a=atom[m.first];
        addGhost(a);
      }
    }

    void Widom::check(UnitTest &test) {
      test.check("widomExcessChemicalPotential", muex() );
    }

    string Widom::info() {
      using namespace Faunus::textio;
      char w=30;
      std::ostringstream o;
      o << header("Multi Particle Widom Analysis")
        << pad(SUB,w, "Number of insertions") << expsum.cnt << endl
        << pad(SUB,w, "Excess chemical pot.") << muex() << kT << endl
        << pad(SUB,w, "Mean activity coefficient") << gamma() << endl
        << pad(SUB,w, "Ghost particles");
      for (auto p : g)
        o << atom[p.id].name << " ";
      o << endl;
      return o.str();
    }

    double Widom::gamma() {
      return exp(muex());
    }

    double Widom::muex() {
      return -log(expsum.avg())/g.size();
    }

    //--------------------------------------------------------------------------------

    WidomScaled::WidomScaled(int insertions){
      cnt=0;
      ghostin = insertions;
    }

    void WidomScaled::add(particle p) {
      g.push_back(p);
      init();
    }

    void WidomScaled::add(Space &c) {
      std::map<short,bool> map;
      for (auto p : c.p)
        map[ p.id ] = true;
      for (auto &m : map) {
        particle a;
        a=atom[m.first];
        add(a);
      }
    }

    bool WidomScaled::overlap(particle &a, particle &b, Space &c) {
      double s=a.radius+b.radius;
      return (c.geo->sqdist(a,b)<s*s) ? true : false;
    }

    void WidomScaled::init() {
      int gspec=g.size();
      chel.resize(gspec);
      chhc.resize(gspec);
      chex.resize(gspec);
      chtot.resize(gspec);
      ewden.resize(gspec);
      ewnom.resize(gspec);
      chint.resize(gspec);
      expuw.resize(gspec);
      chexw.resize(gspec);
      ihc.resize(gspec);
      irej.resize(gspec);

      for(int i=0;i<gspec;i++){
        chel[i]=0.;
        chhc[i]=0.;
        chex[i]=0.;
        chtot[i]=0.;
        ihc[i]=0;
        ewden[i].resize(11);
        ewnom[i].resize(11);
        chint[i].resize(11);
        for(int j=0; j<11; j++) {
          ewden[i][j]=0.;
          ewnom[i][j]=0.;
          chint[i][j]=0.;
        }
      }
    }

    void WidomScaled::insert(Space &c, double lB) {
      particle ghost;
      double u,cu;
      cnt++;
      for(int i=0; i < ghostin; i++) {
        c.geo->randompos(ghost);
        int goverlap=0;
        for(size_t k=0; k < g.size(); k++) {
          ghost.radius = g[k].radius;
          irej[k]=0;
          int j=0;
          while(!overlap(ghost,c.p[j],c) && j<(int)c.p.size())
            j++;
          if( j != (int)c.p.size()) {
            ihc[k]++;
            irej[k]=1;
            goverlap++;
          }
        }

        if ( goverlap != (int)g.size() ) {
          cu=0.;
          u=0.;  //el. potential
          for (size_t l=0; l < c.p.size(); l++) {
            double invdi=1./c.geo->dist(ghost,c.p[l]);
            cu+=invdi;
            u+=invdi*c.p[l].charge;
          } 
          cu=cu*lB;
          u=u*lB;
          double ew,ewla,ewd;
          for(size_t k=0; k < g.size(); k++) {
            if(irej[k]==0) {
              expuw[k]+=exp(-u*g[k].charge);
              for(int cint=0; cint<11; cint++) {
                ew=g[k].charge*(u-double(cint)*0.1*g[k].charge*cu/double(c.p.size()));
                ewla = ew*double(cint)*0.1;
                ewd=exp(-ewla);
                ewden[k][cint]+=ewd;
                ewnom[k][cint]+=ew*ewd;
              }
            }
          }
        }
      }
    }

    string WidomScaled::info() {
      std::ostringstream o;
      double aint4,aint2,aint1;
      for(size_t i=0; i<g.size(); i++) {
        for(int cint=0; cint<11; cint++) {
          if(ewden[i][cint]==0)
            std::cout << "# WARNING: Widom denominator equals zero" << endl;
          else
            chint[i][cint]=ewnom[i][cint]/ewden[i][cint];
        }
        aint4=chint[i][1]+chint[i][3]+chint[i][5]+chint[i][7]+chint[i][9];
        aint2=chint[i][2]+chint[i][4]+chint[i][6]+chint[i][8];
        aint1=chint[i][0]+chint[i][10];  
        chel[i]=1./30.*(aint1+2*aint2+4*aint4);
      }

      double cnttot;
      cnttot=cnt*ghostin;
      o << endl
        << "# SINGLE PARTICLE WIDOM ANALYSIS: (w. charge scaling)" << endl
        << "#   Reference:             doi:10.1080/00268978800100203" << endl
        << "#   Number of Insertions = " << cnttot << endl
        << "#   Excess chemical potentials (kT):" << endl
        << "#           total   elec.  hs            z       r"<< endl;
      for (size_t i=0; i < g.size(); i++) {
        chhc[i]=-log(double(cnttot-ihc[i])/cnttot);
        chexw[i]=-log(expuw[i]);
        chex[i]=chhc[i]+chel[i];
        o.unsetf( std::ios_base::floatfield );
        o << "#   [" << i << "] "
          << std::setprecision(4)
          << std::setw(9) << chex[i]
          << std::setw(9) << chel[i]
          << std::setw(9) << chhc[i]
          << std::setprecision(2) << std::fixed
          << std::setw(9) << g[i].charge
          << std::setw(9) << g[i].radius << endl;
      }
      return o.str();
    }
  }//namespace
}//namespace
