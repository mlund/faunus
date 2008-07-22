#include "peptide.h"

aminoacid::aminoacid() {
  int i;
  d.resize(LAST);
  //ID   RADIUS           DEPROT.CHG.    PROT.CONST.    HYDROPHOBIC?      NAME
  i=ALA; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=false; d[i].name="ALA";
  i=ARG; d[i].radius=3.5; d[i].charge=0; d[i].pka=12.0; d[i].hydrp=false; d[i].name="ARG";
  i=ASN; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=false; d[i].name="ASN";
  i=ASP; d[i].radius=3.5; d[i].charge=-1;d[i].pka=4.8;  d[i].hydrp=false; d[i].name="ASP";
  i=CYS; d[i].radius=3.5; d[i].charge=-1;d[i].pka=10.8; d[i].hydrp=true;  d[i].name="CYS";
  i=GLN; d[i].radius=3.5; d[i].charge=0; d[i].pka=0;    d[i].hydrp=false; d[i].name="GLN";
  i=GLU; d[i].radius=3.5; d[i].charge=-1;d[i].pka=4.8;  d[i].hydrp=false; d[i].name="GLU";
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
  
  i=A37; d[i].radius=3.5; d[i].charge=-1; d[i].pka=6.5;    d[i].hydrp=false; d[i].name="UNK";
  i=CATION;d[i].radius=1.6; d[i].charge=+1; d[i].pka=0;    d[i].hydrp=false; d[i].name="CATION";
};

aminoacid::id aminoacid::getId(string name) {
  //titratable acids
  if (name=="TYR") return TYR;
  if (name=="HIS") return HIS;
  if (name=="CYS") return CYS;
  if (name=="ASP") return ASP;
  if (name=="GLU") return GLU;
  if (name=="LYS") return LYS;
  if (name=="ARG") return ARG;
  if (name=="NTR") return NTR;
  if (name=="CTR") return CTR;

  //hydrophobic
  if (name=="VAL") return VAL;
  if (name=="ILE") return ILE;
  if (name=="LEU") return LEU;
  if (name=="MET") return MET;
  if (name=="PHE") return PHE;
  if (name=="TRP") return TRP;
  if (name=="CYS") return CYS;

  //misc
  if (name=="ALA") return ALA;
  if (name=="ASN") return ASN;
  if (name=="GLN") return GLN;
  if (name=="GLY") return GLY;
  if (name=="PRO") return PRO;
  if (name=="SER") return SER;
  if (name=="THR") return THR;
  if (name=="CATION") return CATION;
  if (name=="A37") return A37;

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


/* P E P T I D E  C L A S S */

//total charge of particle vector
double peptide::charge(vector<particle> &p) {
  double q=0;
  for (int i=0; i<p.size(); i++)
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
