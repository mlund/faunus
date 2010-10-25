#include "faunus/analysis.h"
#include "faunus/physconst.h"
#include "faunus/species.h"

namespace Faunus {

  bool analysis::runtest() { return (slp.runtest(runfraction)); };

  void systemenergy::initialize(double u_init) {
    u0=u_init;
    sum=u0;
    cur=u0;
    confu.clear();
    uavg.reset();
    update(u0);
  }

  //! \param u_init Initial total system energy
  systemenergy::systemenergy(double u_init) {
    initialize(u_init);
  }

  systemenergy::systemenergy() {};

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
  
  void systemenergy::operator+=(double du) {
    sum+=du;
  }
  
  void systemenergy::check(checkValue &test) {
    test.check("systemEnergyAverage", uavg.avg(), 0.3);
    test.check("systemEnergyDrift", std::abs(drift()), 30.);
  }
  
  string systemenergy::info() {
    write();                  //!< Print dynamics of system energy 
    std::ostringstream o;
    o << endl << "# SYSTEM ENERGY (kT):" << endl;
    if (uavg.cnt>0) {
      o << "#   Average <U> s      = " << uavg.avg() << " " << uavg.stdev() << "\n"
        << "#   Initial energy     = " << u0 << endl
        << "#   Initial + changes  = " << sum << endl
        << "#   Current energy     = " << cur << endl;
      o.precision(4);
      o << "#   Total energy drift = " << drift() << " (" << drift()/u0*100. << "\%)" << endl;
    }
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
  
  double systemenergy::drift() {
    return cur-sum;
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
  
  /*! 
   * cntinfo and cntwrite print the counter of each averaged 
   * value in the distribution to a string and write
   * a distribution function of this
   */
  string distributions::cntinfo() {
    unsigned char i;
    std::ostringstream o;
    o << "# DISTRIBUTION FUNCTIONS:\n"
      << "# 1 = distance" << endl;
    for (i=0; i<s.size(); i++)
      o << "# " << i+2 << " = " << s[i] << endl;
    for (float x=xmin; x<=xmax; x+=dx) {
      o << x;
      for (i=0; i<d.size(); i++)
        o << " " << d[i](x).cnt;
      o << endl;
    }
    return o.str();
  }

  bool distributions::cntwrite(string name)
  {
    return fio.writefile(name, cntinfo() );
  }

  //AGGREGATION
  aggregation::aggregation(container &C, vector<macromolecule> &G, double s) : avgintdist(2.0,0.0,600) {
    con=&C;
    g.clear();
    for (int iter=0;iter<G.size();iter++)
      g.push_back(&(G[iter]));
    CNT=0;
    dist.clear();
    dist.resize(g.size(), 0);
    RG2.clear();
    average<double> start;
    histogram a(2.0,0.0,600);
  //  avgintdist=a;
    for (int i=0;i<g.size();i++) {
      RG2.push_back(start);
      intdist.push_back(a);
    }
    sep=s;  
  }

  aggregation::aggregation(container &C, vector<polymer> &G, double s) : avgintdist(2.0,0.0,600) {
    con=&C;
    g.clear();
    for (int iter=0;iter<G.size();iter++)
      g.push_back(&(G[iter]));
    CNT=0;
    dist.clear();
    dist.resize(g.size(), 0);
    RG2.clear();
    contact.clear();
    average<double> start;
    histogram a(2.0,0.0,600);
    for (int i=0;i<g.size();i++) {
      RG2.push_back(start);
      contact.push_back(start);
      intdist.push_back(a);
    }
    sep=s;  
  }
/*
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
          for (int k=i; k<agg.size(); k++) 
            for (int l=agg[k]->beg; l<agg[k]->end+1; l++) {
              intdist[agg.size()-1].add(con->dist(con->p[j], con->p[l]));
              avgintdist.add(con->dist(con->p[j], con->p[l]));
            }
        }
      }
    }
  }*/
  void aggregation::count() {
    ++CNT;
    remaining=g;
    int found=0;
    double contcnt;
    while (remaining.size()>0) { 
      agg.clear();
      unagg.clear();
      if(remaining.size()>1){
        for (int iter=1;iter<remaining.size();iter++)
          unagg.push_back(remaining[iter]);
      }
      agg.push_back(remaining[0]);
      for (int i=0; i<agg.size(); i++){
        for (int j=0; j<unagg.size(); j++)
          if (coll.overlap(con->p, *agg[i], *unagg[j], sep) ) {   // Find the aggregated
            agg.push_back(unagg[j]);                              // pairs
            unagg.erase(unagg.begin()+j);                   // Remove any found from
          }
      } 
      dist[agg.size()-1]+=agg.size();                                     // One aggregate of size
      found+=agg.size();
      point CMagg, fix;
      CMagg.x=0, CMagg.y=0, CMagg.z=0;
      for (int i=1; i<agg.size(); i++) {
        fix=con->vdist(agg[i]->cm,agg[0]->cm);
        CMagg+=fix;
      }
      CMagg+=agg[0]->cm;
      con->boundary(CMagg);
      for (int i=0; i<agg.size();i++) {
        for (int j=agg[i]->beg; j<agg[i]->end+1; j++) {
          RG2[agg.size()-1]+=con->sqdist(con->p[j], CMagg);
          for (int k=i; k<agg.size(); k++) { 
            contcnt=0.0;
            for (int l=agg[k]->beg; l<agg[k]->end+1; l++) {
              if( (l!=j && k==i) || k!=i) {
                intdist[agg.size()-1].add(con->dist(con->p[j], con->p[l]));
                avgintdist.add(con->dist(con->p[j], con->p[l]));
              }
              if( k!=i)
                if( con->dist(con->p[j], con->p[l]) < con->p[j].radius+con->p[l].radius+sep)
                  contcnt++;
            }
            contact[agg.size()-1].add(contcnt);
          }
        }
      }
      for(int j=remaining.size()-1; j>=0; j--)
        for(int i=0; i<agg.size(); i++)
          if(agg[i]==remaining[j])
            remaining.erase(remaining.begin()+j);                   // Remove any found from
    }
    if(found!=g.size())
      std::cerr<<"Didn't find all! Malfunction in aggregation::count()!!!"<<std::endl;
  }
  void aggregation::write(string file) {
    std::ofstream f(file.c_str());
    avgintdist.write("aggpofr-avg.dat");
    if (f) {
      for (int i=0;i<g.size();i++) {
        if (dist[i]!=0){
          std::ostringstream o;
          o << "aggpofr-"<<i+1<<".dat";
          f << i+1 << "  "<<double(dist[i])/CNT/g.size() <<"  " <<pow(RG2[i].avg(),0.5) <<"  "<<contact[i].avg()<<std::endl;
          intdist[i].write(o.str());
        }
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
  
  void virial::check(checkValue &test) {
    test.check("virialExcessPressure", pex.avg() );
  }

  string virial::info() {
    std::ostringstream o;
    if (runfraction>1e-4) {
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
    }
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
  
  string twostatebinding::info() {
    return string(0);
  }
  
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
  
  diskoverlap::diskoverlap(vector<point> &p) {
    s1=19, s2=10, s3=p.size(), cnt=0;
    scnt.resize(19), sscale.resize(19);
    acnt.resize(10), ascale.resize(10);
    size.resize(19), asym.resize(10);
    for (i=0; i<19; i++) {
      sscale[i]=pow(2.,-9.+double(i));
      scnt[i]=0;
    }
    for (i=0; i<10; i++) {
      ascale[i]=sqrt(0.1*double(i+1));
      acnt[i]=0;
    }
    origin.x=0, origin.y=0, origin.z=0, dummy.z=0;
  }
  void diskoverlap::check(vector<point> &p) {
    cnt+=1.0;
    //Convention, major axis x, minor axis y
    for (i=0; i<s1; i++) //Check size 
      for (j=0; j<s3; j++) {
        dummy.x=p[j].x/sscale[i];
        dummy.y=p[j].y/sscale[i];
        if( origin.sqdist(dummy) < 1.0) {
          scnt[i]+=1.0;
          j=s3; // i=s1; // Break loop
        }
      }
    for (i=0; i<s2; i++) //Check assymetry
      for (j=0; j<s3; j++) {
        dummy.x=p[j].x/ascale[i];
        dummy.y=p[j].y*ascale[i];
        if( origin.sqdist(dummy) < 1.0) {
          acnt[i]+=1.0;
          j=s3; // i=s2; // Break loop
        }
      }
  }
  void diskoverlap::blockavg() {
    for (i=0; i<s1; i++){
      size[i]+=scnt[i]/cnt;
      scnt[i]=0;
    }
    for (i=0; i<s2; i++){
      asym[i]+=acnt[i]/cnt;
      acnt[i]=0;
    }
    cnt=0;
  }
  string diskoverlap::info() {
    ostringstream o;
    o << endl << "# HOLE ANALYSIS OF ELLIPSES IN A POINT LATTICE:" << endl
              << "#   Info: A testprogram for Bernard Carbane. "<< endl
              << "#   Data is presented as probability of overlaps"<<endl
              << "#____________________________________________________________"<<endl
              << "#   Scaled circle radius (Rel. sigma), P_rej, stdev   "<<endl
              << "#------------------------------------------------------------"<<endl
              << "#   "<<sscale[ 0]<<"  "<<size[ 0].avg()<<"  "<<size[ 0].stdev()<<" "<<endl
              << "#   "<<sscale[ 1]<<"  "<<size[ 1].avg()<<"  "<<size[ 1].stdev()<<" "<<endl
              << "#   "<<sscale[ 2]<<"  "<<size[ 2].avg()<<"  "<<size[ 2].stdev()<<" "<<endl
              << "#   "<<sscale[ 3]<<"  "<<size[ 3].avg()<<"  "<<size[ 3].stdev()<<" "<<endl
              << "#   "<<sscale[ 4]<<"  "<<size[ 4].avg()<<"  "<<size[ 4].stdev()<<" "<<endl
              << "#   "<<sscale[ 5]<<"  "<<size[ 5].avg()<<"  "<<size[ 5].stdev()<<" "<<endl
              << "#   "<<sscale[ 6]<<"  "<<size[ 6].avg()<<"  "<<size[ 6].stdev()<<" "<<endl
              << "#   "<<sscale[ 7]<<"  "<<size[ 7].avg()<<"  "<<size[ 7].stdev()<<" "<<endl
              << "#   "<<sscale[ 8]<<"  "<<size[ 8].avg()<<"  "<<size[ 8].stdev()<<" "<<endl
              << "#   "<<sscale[ 9]<<"  "<<size[ 9].avg()<<"  "<<size[ 9].stdev()<<" "<<endl
              << "#   "<<sscale[10]<<"  "<<size[10].avg()<<"  "<<size[10].stdev()<<" "<<endl
              << "#   "<<sscale[11]<<"  "<<size[11].avg()<<"  "<<size[11].stdev()<<" "<<endl
              << "#   "<<sscale[12]<<"  "<<size[12].avg()<<"  "<<size[12].stdev()<<" "<<endl
              << "#   "<<sscale[13]<<"  "<<size[13].avg()<<"  "<<size[13].stdev()<<" "<<endl
              << "#   "<<sscale[14]<<"  "<<size[14].avg()<<"  "<<size[14].stdev()<<" "<<endl
              << "#   "<<sscale[15]<<"  "<<size[15].avg()<<"  "<<size[15].stdev()<<" "<<endl
              << "#   "<<sscale[16]<<"  "<<size[16].avg()<<"  "<<size[16].stdev()<<" "<<endl
              << "#   "<<sscale[17]<<"  "<<size[17].avg()<<"  "<<size[17].stdev()<<" "<<endl
              << "#   "<<sscale[18]<<"  "<<size[18].avg()<<"  "<<size[18].stdev()<<" "<<endl
              << "#____________________________________________________________"<<endl
              << "#   Assymetry analysis in terms of half-axis quota a/b of a unit ellips"<<endl
              << "#   a/b, P-rej, stdev "<<endl
              << "#------------------------------------------------------------"<<endl
              << "#   0.1   "<<asym[ 0].avg()<<"  "<<asym[ 0].stdev()<<" "<<endl
              << "#   0.2   "<<asym[ 1].avg()<<"  "<<asym[ 1].stdev()<<" "<<endl
              << "#   0.3   "<<asym[ 2].avg()<<"  "<<asym[ 2].stdev()<<" "<<endl
              << "#   0.4   "<<asym[ 3].avg()<<"  "<<asym[ 3].stdev()<<" "<<endl
              << "#   0.5   "<<asym[ 4].avg()<<"  "<<asym[ 4].stdev()<<" "<<endl
              << "#   0.6   "<<asym[ 5].avg()<<"  "<<asym[ 5].stdev()<<" "<<endl
              << "#   0.7   "<<asym[ 6].avg()<<"  "<<asym[ 6].stdev()<<" "<<endl
              << "#   0.8   "<<asym[ 7].avg()<<"  "<<asym[ 7].stdev()<<" "<<endl
              << "#   0.9   "<<asym[ 8].avg()<<"  "<<asym[ 8].stdev()<<" "<<endl
              << "#   1.0   "<<asym[ 9].avg()<<"  "<<asym[ 9].stdev()<<" "<<endl
              << "#   "<<endl;
    return o.str();
  }

  /*!
   * \param name1 Name of particle 1
   * \param name2 Name of particle 2
   * \param threshold Threshold in Angstroms that defines the pair
   *
   * The particles names and properties must be defined in the species class (i.e. faunatoms.dat file).
   */
  pairing::pairing(string name1, string name2, double threshold) {
    pid1=atom[name1].id;
    pid2=atom[name2].id;
    r2=threshold*threshold;
    assert(pid1!=pid2);
  }

  /*!
   * \param in Inputfile object to search for input keywords
   *
   * The following keywords are searched for: "pair_type1", "pair_type2" and "pair_threshold".
   */
  pairing::pairing(inputfile &in) {
    pid1=atom[ in.getstr("pair_type1") ].id;
    pid2=atom[ in.getstr("pair_type2") ].id;
    r2=pow( in.getflt("pair_threshold", 10 ), 2);
    assert(pid1!=pid2);
  }

  void pairing::sample(container &c) {
    double v=c.getvolume();
    int pairs=0,        // number of pairs
        n1=0,           // number of species 1
        n2=0,           // number of species 2
        n=c.p.size();   // number of particles in container
    
    // count total number of pairing species (which may fluctuate)
    for (int i=0; i<n; ++i) {
      if (c.p[i].id==pid1) n1++;
      else if (c.p[i].id==pid2) n2++;
    }
    // count pairs within threshold
    for (int i=0; i<n-1; ++i) {
      for (int j=i; j<n; ++j) {
        if (c.p[i].id==pid1)
          if (c.p[j].id==pid2)
            if (c.sqdist(c.p[i],c.p[j])<r2)
              pairs++;
        if (c.p[i].id==pid2)
          if (c.p[j].id==pid1)
            if (c.sqdist(c.p[i],c.p[j])<r2)
              pairs++;
      }
    }
    // average concentrations
    cpair+=pairs/v;    // pairs
    c1+=(n1-pairs)/v;  // free species 1
    c2+=(n2-pairs)/v;  // free species 2
  }
  
  string pairing::info() {
    std::ostringstream o;
    if (c1.cnt>0 && c2.cnt>0) {
      double tomM=1e30/pyc.Nav; // A^-3 -> mmol/l.
      double AB=tomM*cpair.avg();
      double A=tomM*c1.avg();
      double B=tomM*c2.avg();
      o << endl
      << "# PAIRING ANALYSIS:\n"
      << "#   Number of samples           = " << cpair.cnt << endl
      << "#   Pairing particles           = " << atom[pid1].name << " " << atom[pid2].name << endl
      << "#   Pair threshold (A)          = " << sqrt(r2) << endl
      << "#   Concentrations: A B AB (mM) = " << A << " " << B << " " << AB << endl
      << "#   Association constant (1/mM) = " << AB/(A*B) << endl
      << "#   Dissociation constant (mM)  = " << A*B/AB << endl;
    }
    return o.str();
  }
  
  
  /*!
   * \param dR Width (in Angstroms) of the shell to analyse near the boundary. The smaller the more accurate.
   */
  osmoticpressure::osmoticpressure(cell &c) {
    cnt=0;
    cPtr=&c;
    hist.resize( int(c.r+.5) );
  }

  /*!
   * \param c Spherical simulation container
   * \param g Group containing all mobile ions
   */
  void osmoticpressure::sample(group &g) {
    int k;
    cnt++;
    for (int i=g.beg; i<=g.end; i++) {
      k=cPtr->p.at(i).len()+.5;
      assert(k<hist.size());
      hist[k]++;
    }
    rhoid += g.size() / cPtr->getvolume();
  }

  double osmoticpressure::getConc(double r) {
    assert(int(r+.5) < hist.size() );
    double v=4/3.*pyc.pi*( pow(r+.5,3) - pow(r-.5,3) );
    return hist.at( int(r+.5) ) / v / cnt * 1e27/pyc.Nav;
  }

  string osmoticpressure::info() {
    char w=11;
    double toM=1e27/pyc.Nav; // atoms/aa3 to mol/l
    double toPa=pyc.kB * pyc.T * pyc.Nav; // mol/l to kPa
    double id=rhoid.avg()*toM,
           tot=getConc( hist.size()-1 ), // get last point in histogram
           ex=tot-id;
    std::ostringstream o;
    o << endl
      << "# CELL MODEL PRESSURE ANALYSIS:" << endl
      << "#   More information:        doi:10.1063/1.443547" << endl
      << "#   Number of samples      = " << rhoid.cnt << endl;
    o.precision(4);
    if (rhoid.cnt>0) {
      o << "#   Osmotic coefficient    = " << tot/id << endl << std::left
        << "#                            " << setw(w) << "ideal" << setw(w) << "excess" << setw(w) << "total"  << endl
        << "#   Pressure (mM)          = " << setw(w) << id*1000 << setw(w) << ex*1000  << setw(w) << tot*1000 << endl
        << "#    - / / - (kPa)         = " << setw(w) << id*toPa << setw(w) << ex*toPa  << setw(w) << tot*toPa << endl;
    }
    return o.str();
  }

  bool osmoticpressure::write(string file) {
    std::ofstream f(file.c_str());
    if (f) {
      f << "# Radial oncentration profile for mobile ions (cell model osm. pressure)\n"
        << "# x - distance from cell origin\n"
        << "# y - concentration in mol/l\n";
      for (int r=0; r<hist.size(); r++)
        f << r << " " << getConc(r)*1000 << endl;
      f.close();
      return true;
    }
    return false;
  }

  void osmoticpressure::check(checkValue &test) {
    double tot=getConc( hist.size()-1 );
    test.check("osmotic_pressure", tot, 0.1);
    test.check("osmotic_coefficient", tot/rhoid.avg(), 0.1);
  }

}//namespace

