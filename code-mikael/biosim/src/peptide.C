#include "peptide.h"

aminoacid::aminoacid(double bjerrum) {
  int i;
  lB=bjerrum;
  d.resize(LAST);
  //ID   RADIUS           DEPROT.CHG.    PROT.CONST.    HYDROPHOBIC?      NAME
  i=ALA; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=false; d[i].name="ALA";
  i=ARG; d[i].radius=3.5; d[i].charge=0; d[i].pka=12.0; d[i].hydrp=false; d[i].name="ARG";
  i=ASN; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=false; d[i].name="ASN";
  i=ASP; d[i].radius=3.5; d[i].charge=-1;d[i].pka=4.0;  d[i].hydrp=false; d[i].name="ASP";
  i=CYS; d[i].radius=3.5; d[i].charge=-1;d[i].pka=10.8; d[i].hydrp=true;  d[i].name="CYS";
  i=GLN; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=false; d[i].name="GLN";
  i=GLU; d[i].radius=3.5; d[i].charge=-1;d[i].pka=4.4;  d[i].hydrp=false; d[i].name="GLU";
  i=GLY; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=false; d[i].name="GLY";
  i=HIS; d[i].radius=3.5; d[i].charge=0; d[i].pka=6.3;  d[i].hydrp=false; d[i].name="HIS";
  i=ILE; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=true;  d[i].name="ILE";
  i=LEU; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=true;  d[i].name="LEU";
  i=LYS; d[i].radius=3.5; d[i].charge=0; d[i].pka=10.4; d[i].hydrp=false; d[i].name="LYS";
  i=MET; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=true;  d[i].name="MET";
  i=PHE; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=true;  d[i].name="PHE";
  i=PRO; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=false; d[i].name="PRO";
  i=SER; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=false; d[i].name="SER";
  i=THR; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=false; d[i].name="THR";
  i=TRP; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=true;  d[i].name="TRP";
  i=TYR; d[i].radius=3.5; d[i].charge=-1;d[i].pka=9.6;  d[i].hydrp=true;  d[i].name="TYR";
  i=VAL; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=false; d[i].name="VAL";
  i=CTR; d[i].radius=3.5; d[i].charge=-1;d[i].pka=3.8;  d[i].hydrp=false; d[i].name="CTR";
  i=NTR; d[i].radius=3.5; d[i].charge=0; d[i].pka=7.5;  d[i].hydrp=false; d[i].name="NTR";
  i=UNK; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=false; d[i].name="UNK";
  i=CATION;d[i].radius=1.6; d[i].charge=+1; d[i].pka=0; d[i].hydrp=false; d[i].name="CATION";
  i=NA;d[i].radius=1.6; d[i].charge=+1; d[i].pka=0;     d[i].hydrp=false; d[i].name="NA";
  i=K;d[i].radius=1.8;  d[i].charge=+1; d[i].pka=0;     d[i].hydrp=false; d[i].name="K";
  i=CL;d[i].radius=1.7; d[i].charge=-1; d[i].pka=0;     d[i].hydrp=false; d[i].name="CL";
  i=BR;d[i].radius=1.7; d[i].charge=-1; d[i].pka=0;     d[i].hydrp=false; d[i].name="BR";
  i=I;d[i].radius=1.7; d[i].charge=-1; d[i].pka=0;     d[i].hydrp=false; d[i].name="I";
  i=RCOO;d[i].radius=1.7; d[i].charge=-1; d[i].pka=0; d[i].hydrp=false; d[i].name="RCOO";
  i=RCOOH;d[i].radius=1.7; d[i].charge=0; d[i].pka=0; d[i].hydrp=false; d[i].name="RCOOH";
  i=RNH3;d[i].radius=1.7; d[i].charge=0; d[i].pka=0; d[i].hydrp=false; d[i].name="RNH3";
  i=RNH4;d[i].radius=1.7; d[i].charge=1; d[i].pka=0; d[i].hydrp=false; d[i].name="RNH4";
  i=GHOST;d[i].radius=0; d[i].charge=0; d[i].pka=0; d[i].hydrp=false; d[i].name="GHOST";
};

aminoacid::id aminoacid::getId(string name) {
  for (int i=0; i<d.size(); i++)
    if (d[i].name==name)
      return id(i);
  return UNK;
};

double aminoacid::volume(double weight, double density=1.) {
  return 1.6606*density*weight;
};

double aminoacid::radius(double weight, double density=1.) {
  double vol = volume(weight, density);
  return pow(3.*vol/4./3.14, 0.33333333);
};

double aminoacid::vdW(aminoacid::id iden) {
  return 25000.;
};

bool aminoacid::loadpmf(string filename) {
  string s,a_str,b_str;
  int a,b,len,tmp;
  vector<double> x,y;
  ifstream f(filename.c_str());
  if (f) {
    //scan file, word-by-word
    while (!f.eof()) {
      f >> s;
      if (s.find("#$")==0) {
        s.clear();
        f >> a_str >> b_str >> len;
        x.resize(len);
        y.resize(len);
        a = getId(a_str);
        b = getId(b_str);
        if (a>b)
          swap(a,b);
        for (int i=0; i<len; i++)
          f >> x[i] >> y[i];
        if (pairpot[a][b].xmax==0) {
          pairpot[a][b].add( x, y );
          pairpot[a][b].comment=filename;
        };
      };
    };
    f.close();
    return true;
  }
  return false;
};

void aminoacid::loadpmf(string pmfdir, string type) {
  string n1,n2;
  unsigned int i,j=getId(type);
  for (i=FIRST; i<LAST; i++) {
    n1=d[i].name+"-"+d[j].name;
    n2=d[j].name+"-"+d[i].name;
    if (loadpmf( pmfdir + n1 + ".dat")==false)
      loadpmf( pmfdir + n2 + ".dat");
  }
};

void aminoacid::showpmf() {
  cout << "# --- LOADED PMF's ----------------------------------\n";
  cout << "# (a,b,resolution,r_max,file)\n";
  for (unsigned int i=FIRST; i<LAST; i++)
    for (unsigned int j=FIRST; j<LAST; j++)
      if (pairpot[i][j].xmax>0)
        cout << "# " << d[i].name << " " << d[j].name << " "
          << pairpot[i][j].res << " "
          << pairpot[i][j].xmax << " "
          << pairpot[i][j].comment << endl;
  cout << endl;
}

// internal energy
double aminoacid::energy(vector<particle> &p) {
  double u=0;
  unsigned int i,j,n=p.size();
  for (i=0; i<n-1; i++)
    for (j=i+1; j<n; j++)
      u+=energy(p[i], p[j]);
  return u; 
};

double aminoacid::energy(vector<particle> &p, unsigned int j) {
  if (p[j].id==GHOST) return 0;
  double u=0;
  unsigned int i,n=p.size();
  for (i=0; i<j; i++) u+=energy(p[i],p[j]);
  for (i=j+1; i<n; i++) u+=energy(p[i],p[j]);
  return u;
}

double aminoacid::energy(vector<particle> &p, group &g) {
  unsigned int glen=g.end+1, n=p.size();
  double u=0;
  if (g.beg!=-1)
    for (unsigned int i=g.beg; i<glen; i++) {
      for (unsigned int j=0; j<g.beg; j++)
        u+=energy(p[i],p[j]);
      for (unsigned int j=glen; j<n; j++)
        u+=energy(p[i],p[j]);
    };
  return u;
};

double aminoacid::energy(vector<particle> &p, group &g1, group &g2) {
  if (g1.beg==-1 || g1.beg==-1)
    return 0.;
  double u=0;
  unsigned int ilen=g1.end+1;
  unsigned int jlen=g2.end+1;
  if (g1.cm!=-1 && ilen>g1.cm) ilen=g1.cm;  //avoid ghost particles
  if (g2.cm!=-1 && jlen>g2.cm) jlen=g2.cm;  //such as "cm" and "dipv"
  for (unsigned int i=g1.beg; i<ilen; i++)
    for (unsigned int j=g2.beg; j<jlen; j++)
      u += energy(p[i],p[j]);
  return u;
}


/* P E P T I D E  C L A S S */

//total charge of particle vector
double peptide::charge(vector<particle> &p) {
  double q=0;
  for (unsigned int i=0; i<p.size(); i++)
    q+=p[i].charge;
  return q;
};

//center of mass of particle vector
point peptide::center_of_mass(vector<particle> &p) {
  particle cm;
  double x=0,y=0,z=0,m=0;
  for (int i=0; i<p.size(); i++) {
    x+=p[i].x*p[i].mw;
    y+=p[i].y*p[i].mw;
    z+=p[i].z*p[i].mw;
    m+=p[i].mw;
  };
  cm.x = x/m;
  cm.y = y/m;
  cm.z = z/m;
  return cm;
};

/*
 * Loads an Amino Acid Model coordinate file.
 * (where residues are represented by a spheres)
 */
vector<particle> peptide::loadAAModel(string filename, keys k) {
  vector<particle> p(0);
  point c;
  int n, num;
  string s,res;
  ifstream f(filename.c_str());
  if (f) {
    do getline(f,s);
    while (s.find("#")!=string::npos);
    n=atoi(s.c_str());
    p.resize(n);
    c.y=-6.*n/2.;
    for (int i=0; i<n; i++) {
      f >> res >> num >> p[i].x >> p[i].y >> p[i].z >> p[i].charge >> p[i].mw >> p[i].radius;
      p[i].id = getId(res);
      p[i].hydr = d[p[i].id].hydrp;
      if (k==CHAIN) {
        p[i] = c;
        p[i].radius = 2.0;
        c.y += 6.;
      };
    };
    f.close();
  } else
    cout << "# Error! File '"<<filename<<"' not found.";

  return p;
};

vector<particle> peptide::loadstructure(string filename, keys id) {
  vector<filefmt> l;
  vector<particle> out,aa;
  particle p;
  point cm;
  int cnt,i;

  //load structure file
  ifstream f(filename.c_str());
  if (f) {
    f >> cnt ;
    l.resize(cnt);
    for (int i=0; i<cnt; i++) {
      string unk;
      f >> l[i].atomname >> l[i].mw >> unk >> l[i].charge
        >> l[i].x >> l[i].y >> l[i].z >> l[i].atomnr >> l[i].aminoname
        >> l[i].aminonr;
    };
    f.close();
  } else {
    cout << "*** WARNING: Structure NOT loaded! ***\n";
  };
  //showel data into particle vector
  i=0;
  cnt=l[i].aminonr;
  while (i<l.size()) {
    //go through each acid
    if (cnt==l[i].aminonr) {
      p.charge=l[i].charge;
      p.mw=l[i].mw;
      p.x=l[i].x;
      p.y=l[i].y;
      p.z=l[i].z;
      p.id=getId(l[i].aminoname);
      aa.push_back(p);
      i++;
    };
    if (cnt!=l[i].aminonr) {
      cnt++;
      p.radius=3.5;
      p.id=aa[0].id;
      p.charge=charge(aa);
      cm=center_of_mass(aa);
      p.x=cm.x;
      p.y=cm.y;
      p.z=cm.z;
      out.push_back(p);
      aa.clear();
    };
  };
  cout << "# Structural information: " << endl
    << "#   Filename      = " << filename << endl
    << "#   Residues      = " << out.size() << endl
    << "#   Total charge  = " << charge(out) << endl
    << "#   Protein model = AMINOACID" << endl<<endl;
  return out;
};

bool peptide::saveAAModel(string file, vector<particle> &p, group &g) {
  ofstream f( file.c_str() );
  if (f) {
    f << g.size() << endl;
    for (int i=g.beg; i<=g.end; i++) {
      if (p[i].radius>0) {
        if (p[i].id<d.size())
          f << d[p[i].id].name <<" ";
        else 
          f << "UNK" << " ";
        f << i+1 << " ";
        f << p[i].x << " " << p[i].y << " " << p[i].z << " ";
        f << p[i].charge << " " << p[i].mw << " " << p[i].radius << endl;
      };
    };
    f.close();
    return true;
  };
  return false;
};
