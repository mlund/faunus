#include "inputfile.h"

//! \param filename Input file to scan
inputfile::inputfile(string filename) {
  matrix.resize(50);
  int i=0;
  ifstream f( filename.c_str() );
  if (f) {
    while (!f.eof()) {
      f >> matrix[i].name >> matrix[i].val;
      i++;
    };
    matrix.resize(i-1);
    f.close();
    cout << "# Configuration read from: " << filename << endl;
  } else cout << "*** Failed to open inputfile ***" << endl;
};

int inputfile::findKey(string &key) {
   for (int i=0; i<matrix.size(); i++)
      if (matrix[i].name.compare(key)==0) return i;
   return -1;
};

//! \param key Keyword to look for
//! \param def Default value if keyword is not found
string inputfile::getstr(string key, string def) {
  int i = findKey(key);
  if (i!=-1)
    return matrix[i].val;
  else
    cout << "# Warning: keyword '" << key << "' not found - using "
         << def << endl;
  return def;
};

//! \param key Keyword to look for
//! \param def Default value if keyword is not found
double inputfile::getflt(string key, double def) {
  int i = findKey(key);
  if (i!=-1) return atof(matrix[i].val.c_str());
  else {
    cout << "Keyword '" << key << "' not found! - using "<<def<<"\n";
    return def;
  };
};

//! \param key Keyword to look for
//! \param def Default value if keyword is not found
int inputfile::getint(string key, int def) {
  int i = findKey(key);
  if (i!=-1) return atoi(matrix[i].val.c_str());
  else {
    cout << "Keyword '" << key << "' not found! -using "<<def<<"\n";
    return def;
  };
};

//! \param key Keyword to look for
//! \param def Default value if keyword is not found
bool inputfile::getboo(string key, bool def) {
  int i = findKey(key);
  if (i!=-1) {
    if (matrix[i].val.compare("yes")==0)
      return true;
    else
      return false;
  } else {
    cout << "Keyword '" << key << "' not found! - using "<<def<<"\n";
    return def;
  };
};

// config constructor
config::config(string filename) : inputfile(filename) {

  // integers
  seed  = getint("randomseed", 13);
  macro = getint("macro", 10);
  micro = getint("micro", 1000);
  nion1 = getint("nion1", 0);  //number of ion 1                   
  nion2 = getint("nion2", 0);  //number of ion 2             
  nion3 = getint("nion3", 0);  //number            

  // floats
  cell_r    = getflt("cell_r");
  maxsep    = getflt("maxsep",cell_r/1.5);//restrict protein separation (max)    
  minsep    = getflt("minsep",0.);         //restrict protein separation (min)     
  temp      = getflt("temp", 298);         //temperature (K)
  pH        = getflt("pH", 7);
  vdw       = getflt("hamaker", 0);
  springk   = getflt("springconst",0.1);
  springeq  = getflt("springeqdist",10.);
  u_penalty = getflt("penalty",0.);
  dielec    = getflt("dielec",78);            
  zion1     = getflt("chion1", 1);         //ionic charges
  zion2     = getflt("chion2", -1);               
  zion3     = getflt("chion3");               
  rion3     = getflt("rion3");
  prot_dp   = getflt("prot_dp");              
  prot_rot  = getflt("prot_rot");  
  dp_monomer=getflt("monomer_dp",2.); //monomer displacement factor
  clust_dp  = 10.;
  ion_dp   = getflt("ion_dp"); //ion displacement

  // strings
  jobid    = getstr("jobid", ".a");    //arbitrary job name
  protein1 = getstr("protein1"); //protein 1 coords         
  protein2 = getstr("protein2"); //protein 1 coords    
  tion1    = getstr("tion1", "NA");                
  tion2    = getstr("tion2", "CL");                
  tion3    = getstr("tion3");                
  pmfdir   = getstr("pmfdir");               

  // bools
  rotate      = getboo("rotate", true);         
  adjust_dp   = getboo("adjust_dp", true);
  hairy       = getboo("hairy", true);
  titrateBool = getboo("titrate", true);
  smear       = getboo("smear", false);
  minsnapshot = getboo("minsnapshot",false);
  adjust_dp   = getboo("adjust_dp", false);
  imdBool     = getboo("imdsupport",false);
};

