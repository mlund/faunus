#include "faunus/analysis.h"

namespace Faunus {

  bool analysis::runtest() { return (slp.runtest(runfraction)); };

  //! \param u_init Initial total system energy
  systemenergy::systemenergy(double u_init) {
    u0=u_init;
    sum=u0;
    cur=u0;
    confu.clear();
  }
  void systemenergy::update(double energy) {
    cur=energy;
    uavg+=cur;
  }
  void systemenergy::track() {
    confu.push_back(sum);
  }
  void systemenergy::write() {
    if (confu.size())
      fio.writefile("systemenergy.dat", confuout() );
  }
  void systemenergy::operator+=(double du) { sum+=du; }
  string systemenergy::info() {
    write();                  //!< Print dynamics of system energy 
    std::ostringstream o;
    o << endl << "# SYSTEM ENERGY (kT):" << endl;
    if (uavg.cnt>0)
      o << "#   Average <U> s      = " << uavg.avg() << " " << uavg.stdev() << "\n";
    o << "#   Initial energy     = " << u0 << endl
      << "#   Initial + changes  = " << sum << endl
      << "#   Current energy     = " << cur << endl
      << "#   Absolute drift     = " << std::abs(cur-sum) << endl;
    return o.str();
  }
  string systemenergy::confuout() {
    int j=confu.size();
    std::ostringstream o;
    o << endl << "# SYSTEM ENERGY (kT):"<< endl;
    for (int i=0;i<j;i++)
      o << i+1 << " " << confu[i] << endl;
    return o.str();
  }

  //---------------- ANGULAR CORRELATIONS ---------------------------
  angularcorr::angularcorr() {}
  void angularcorr::update(macromolecule &g1, macromolecule &g2, distributions &d) {
    double l1=g1.mu.len(), l2=g2.mu.len(), z;
    if (l1>0 && l2>0 )  {
      m1=g1.mu*(1./l1);
      m2=g1.mu*(1./l2); // unit vectors at origo
      z=g1.cm.dist(g2.cm);
      d.add("dip_z*z", z, m1.z*m2.z);
      d.add("dip_x*x", z, m1.x*m2.x);
      d.add("dip_y*y", z, m1.y*m2.y);
    }
  }

  //--------------- KIRKWOOD FACTOR ----------------------------------
  gfactor::gfactor(container &c, molecules &m) //: h(float(0.1), float(0.), float(40.))
  {
    N=m.size()/m.numatom;
    V=c.getvolume();
    mol.beg=m[0].beg;
    mol.end=m[0].end;
    mu2=pow(mol.dipole(c.p),2);  //Identical dipoles
    gamma=pow(1.602*10e-19,2)*10e10*N/V*mu2/(8.854*10e-12)/9; //in gaussian units times kT
    gscale=1./N/mu2;
  } //: h(0.1,0.1,40.); 
  double gfactor::add(vector<particle> &p, molecules &m) {
    mol.beg=0;
    mol.end=p.size()-1;
    mol.cm.x=0;
    mol.cm.y=0;
    mol.cm.z=0;
    MU2=pow( mol.dipole(p) , 2);
    g.add(MU2);
  //  h.add(MU2/N/mu2);
    return MU2/N/mu2;
  }
  string gfactor::info() {
  //  h.write("kirkwood.dat");
    std::ostringstream o;
    o << endl << "# KIRKWOOD FACTOR: "<<endl;
            o << "#"<<endl
              << "# "<<N<<" solvent molecules in a "<<V<<" A^3 cavity, rho ="<<N/V<<endl
              << "# Molecular dipole = "<<sqrt(mu2)<<" (mu^2 = "<<mu2<<" )"<<endl     
              << "# MU^2    = " << g.avg() << "(per particle), stdev =  "<<g.stdev()<<endl
              << "# g       = " << g.avg()*gscale<<endl
  //            << "# gamma = " << 4*acos(-1)*N/V*mu2/9*10e10 <<" (times kT)" <<endl
              << "# 3*gamma'= " << 3*gamma<<" (gaussian times kT)"<<endl<<endl;
    return o.str();
  }
  //---------------- DISTRIBUTIONS ----------------------------------
  /*!
   * By default the maximum and minimum x values are dynamically
   * resized to match any given dataset (for x>=0).
   *
   * \note If x<0, "min" should be specified
   * \param deltax x value resolution
   * \param min Minimum x value (optional, auto-resizes)
   * \param max Maximum x value (optional, auto-resizes)
   */
  distributions::distributions(float deltax, float min, float max) {
    dx=deltax;
    xmax_set=max;
    xmin_set=min;
    xmax=-1e7;
    xmin=-xmax;
  }

  /*! 
   * Adds a value to the distribution. The y value is automatically
   * averaged. If "name" is not found in the current distribution,
   * it will be created and the data point will be added.
   * 
   * \param name Name of distribution to add to.
   * \param x x-value
   * \param y value (add to average)
   */
  bool distributions::add( string name, float x, float y )
  {
    unsigned short i = find(name);
    d[i](x)+=y; // add to average
    if (x>xmax) // track max x value
      xmax=x;
    if (x<xmin) // track min x value
      xmin=x;
    return true;
  }

  unsigned short distributions::find( string name ) {
    unsigned short i;
    for (i=0; i<s.size(); i++) // search name vector
      if (s[i]==name)
        return i;
    s.push_back(name); // if not found, add it!
    d.resize( d.size()+1,
        xytable<float, average<float> >(dx,xmin_set,xmax_set));
    return s.size()-1;
  }

  string distributions::info() {
    unsigned char i;
    std::ostringstream o;
    o << "# DISTRIBUTION FUNCTIONS:\n"
      << "# 1 = distance" << endl;
    for (i=0; i<s.size(); i++)
      o << "# " << i+2 << " = " << s[i] << endl;
    for (float x=xmin; x<=xmax; x+=dx) {
      o << x;
      for (i=0; i<d.size(); i++)
        o << " " << d[i](x).avg();
      o << endl;
    }
    return o.str();
  }

  bool distributions::write(string name)
  {
    return fio.writefile(name, info() );
  }

  //aggregation
  aggregation::aggregation(container &C, vector<macromolecule> &G, double s) {
    con=&C;
    g.clear();
    for (int iter=0;iter<G.size();iter++)
      g.push_back(&(G[iter]));
    //g=&G;
    CNT=0;
    dist.clear();
    dist.resize(g.size(), 0);
    RG2.clear();
    average<double> start;
    for (int i=0;i<g.size();i++)
      RG2.push_back(start);
    sep=s;  
  }
  void aggregation::count() {
    ++CNT;
    for (int k=0; k<g.size();k++) {  //start looking for aggregates around k
      //prepare the vectors
      agg.clear();
      unagg.clear();
      for (int iter=0;iter<g.size();iter++)
        if (iter!=k)
          unagg.push_back(g[iter]);
      agg.push_back(g[k]);
      //start looking
      for (int i=0; i<agg.size(); i++){
        for (int j=0; j<unagg.size(); j++)
          if (coll.overlap(con->p, *agg[i], *unagg[j], sep) ) // Find the aggregated
            agg.push_back(unagg[j]);                              // pairs
        for (int i=0; i<agg.size(); i++)   
          for (int j=0; j<unagg.size(); j++)
            if (agg[i]==unagg[j]) {
              unagg.erase(unagg.begin()+j);                   // Remove any found from
              j--;                                            // unagg
            } 
      } 
      //analysis
      dist[agg.size()-1]++;                                     // One aggregate of size
      if (agg.size()+unagg.size()!=g.size()) {
        std::cout <<" Error in aggregation::count()"<<endl
          <<" g.size() = agg.size()+unagg.size()"<<endl
          <<"    "<< g.size()<<" = "<<agg.size()<<" +  "<<unagg.size()<<endl;
      }
      point CMagg, fix;
      for (int i=1; i<agg.size(); i++) {
        fix=agg[i]->cm-agg[0]->cm;
        CMagg+=fix;
      }
      CMagg+=agg[0]->cm;
      con->boundary(CMagg);
      for (int i=0; i<agg.size();i++) {
        for (int j=agg[i]->beg; j<agg[i]->end+1; j++) {
          RG2[agg.size()-1]+=con->sqdist(con->p[j], CMagg);
        }
      }
    }
  }
  void aggregation::write(string file) {
    std::ofstream f(file.c_str());
    if (f) {
      for (int i=0;i<g.size();i++) {
        if (dist[i]!=0)
          f << i+1 << "  "<<double(dist[i])/CNT/g.size() <<"  " <<RG2[i].avg() <<endl;
      }
      f.close();
    }
  }

  virial::virial(container &c) {
    runfraction=1.0;
    dr=0.1;
    conc = c.p.size()/c.getvolume();
  }
  virial::virial(container &c, vector<macromolecule> &g) {
    runfraction=1.0;
    dr=0.01;
    conc = g.size()/c.getvolume();
  }

  void virial::sample(container &c, energybase &pot) {
    if (runtest()) {
      double p=0;
      int n=c.p.size();
#pragma omp parallel for reduction (+:p)
      for (int i=0; i<n-1; i++) {
        point f,rij;
        double r;
        for (int j=i+1; j<n; j++) {
          rij=c.vdist(c.p[i],c.p[j]);
          r=rij.len();
//          f = rij * ( pot.force(c,c.p[i],c.p[j],r,dr) / r );
//          p+= rij.x*f.x + rij.y*f.y + rij.z*f.z;
          p+=pot.force(c, c.p[i], c.p[j],rij,r,dr)*r;
        }
      }
      pex+=p*pot.tokT / (3*c.getvolume());  // kT/AA^3
    }
  }
  void virial::sample(container &c, energybase &pot, vector<macromolecule> &g) {
    if (runtest()) {
      double p=0;
      int n=g.size();
//#pragma omp parallel for reduction (+:p)
      for (int i=0; i<n-1; i++) 
        for (int j=i+1; j<n; j++)
          for (int s=g[i].beg; s<=g[i].end; s++) {
            point f,rij;
            double r;
            for (int t=g[j].beg; t<=g[j].end; t++) {
              rij=c.vdist(c.p[s],c.p[t]);
              r=rij.len();
  //            f = rij * ( pot.force(c,c.p[s],c.p[t],r,dr) / r );
  //            p+= rij.x*f.x + rij.y*f.y + rij.z*f.z;
            }
          }
      pex+=p*pot.tokT / (3*c.getvolume());  // kT/AA^3
    }
  }

  string virial::info() {
    std::ostringstream o;
    char w=10;
    double tokPa=2.487*1e30/6.022e23,
           toM=2.487*1e30/(8.314*298.15*6.022e23),
           pid=conc*1e27/6.022e23,
           ex=pex.avg()*toM;

    o.unsetf( std::ios_base::floatfield );
    o << "\n# VIRIAL ANALYSIS: (Frenkel & Smith 2nd Ed. p.84 / McQuarrie p.338)\n"
      << "#   Number of force calculations = " << pex.cnt << "\n"
      << "#   Differentiation step (AA)    = " << dr << "\n"
      << "#   Run fraction                 = " << runfraction << "\n"
      << std::setprecision(4)
      << "#   Osmotic coefficient          = " << (pid+ex)/pid << "\n"
      << "#                     "
      << setw(w) << "ideal" << setw(w) << "excess" << setw(w) << "total" << setw(w+3) << "ex. stddev\n"
      << "#   Pressure (mol/l): "
      << setw(w) << pid << setw(w) << ex << setw(w) << pid+ex << setw(w) << pex.stdev()*toM
      << endl;
    return o.str();
  }
  //----------------------------------------------

  void pointpotential::add(point p, string name) {
    average<double> dummy;
    data d = {p, name, dummy , dummy };
    d.p.x+=slp.random_half()*1e-9; // avoid overlapping particles.
    list.push_back(d);
  }
  void pointpotential::sample(container &c, energybase &pot) {
    int n=list.size();
#pragma omp for
    for (int i=0; i<n; ++i) {
      double phi=pot.potential(c.p, list[i].p);
      list[i].phi+=phi;
      list[i].expphi+=exp(-phi);
    }
  }
  string pointpotential::info() {
    std::ostringstream o;
    o << endl
      << "# Point Potential Analysis:\n"
      << "#   Number of points         = " << list.size() << endl
      << "#   Number of samples        = "
      << ((list.size()==0) ? 0 : list[0].phi.cnt) << endl;
    for (int i=0; i<list.size(); ++i)
      o << "#    " << list[i].name << " " << list[i].p << " "
        << list[i].phi.avg() << " "
        << list[i].phi.stdev() << "              "
        << -log(list[i].expphi.avg()) << " "
        << -log(list[i].expphi.stdev()) << endl;
    return o.str();
  }

  twostatebinding::twostatebinding(double radius) {
    r2=radius*radius;
    vol=4./3.*acos(-1.)*pow(radius,3);
  }
  void twostatebinding::update(container &con, point &site, group &g) {
    unsigned short i,cnt=0;
    for (i=g.beg; i<=g.end; i++)
      if (con.sqdist(site,con.p[i])<r2) {
        cnt=1;
        break;
      }
    conc+=float(cnt)/vol;
  }
  void twostatebinding::update(container &con, point &site, vector<macromolecule> &g, int first) {
    unsigned short i,j,cnt=0;
    for (j=first; j<g.size(); j++)
      for (i=g[j].beg; i<=g[j].end; i++)
        if (con.sqdist(site, con.p[i])<r2) {
          cnt+=1;
          break;
        }
    conc+=float(cnt)/vol;
  }
  void twostatebinding::update(container &con, point &site, unsigned char id) {
    unsigned short i,n=con.p.size(),cnt=0;
    for (i=0; i<n; i++)
      if (con.p[i].id==id)
        if (con.sqdist(site,con.p[i])<r2) cnt++;
    conc+=float(cnt)/vol;
  }
  string twostatebinding::info() { return string(0); }
  string twostatebinding::info(double bulkconc) {
    std::ostringstream o;
    o << endl << "# TWOSTATE BINDING ANALYSIS:" << endl
      << "#   More information:  J. Phys. Chem. 1995, 99, 10412" << endl
      << "#   Site radius      = " << sqrt(r2) << endl
      << "#   Avg. site conc.  = " << conc.avg() << endl;
    if (bulkconc>0)
      o << "#   Site excess      = " << conc.avg()/bulkconc
        << " (" << -log(conc.avg()/bulkconc) << " kT)" << endl;
    return o.str();
  }
}//namespace

