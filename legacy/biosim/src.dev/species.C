#include "species.h"

string species::info(string name) {
  int i=id(name);
  return name;
}

particle species::get(string name) {
  particle a;
  a.id = id(name);
  a.charge = d[a.id].charge;
  a.radius = d[a.id].radius;
  return a;
  
};

double species::vol(double mw, double rho=1.)
{ return 1.6606*rho*mw; };

double species::radius(double mw, double rho=1.) {
  double v = vol(mw, rho);
  return pow(3.*v/4./3.14, 0.3333);
};

void species::set(int i, string name, float r, float z,
		  float pka, bool hyd) {
  d[i].name = name;
  d[i].radius = r;
  d[i].charge = z;
  d[i].pka = pka;
  d[i].hydrp = hyd;
};

particle::type species::id(string name) {
  for (int i=0; i<d.size(); i++)
    if (d[i].name==name)
      return particle::type(i);
  return particle::UNK;
};

species::species() {
  d.resize(particle::LAST);

  //  id    name  rad   Z   pKa  hydrp
  set(particle::ALA, "ALA", 3.5,  0,  0.0, false);
  set(particle::ARG, "ARG", 3.5,  0, 12.0, false);
  set(particle::ASN, "ASN", 3.5,  0,  0.0, false);
  set(particle::ASP, "ASP", 3.5, -1,  4.0, false);
  set(particle::CYS, "CYS", 3.5, -1, 10.8, true );
  set(particle::GLN, "GLN", 3.5,  0,  0.0, false);
  set(particle::GLU, "GLU", 3.5, -1,  4.4, false);
  set(particle::GLY, "GLY", 3.5,  0,  0.0, false);
  set(particle::HIS, "HIS", 3.5,  0,  6.3, false);
  set(particle::ILE, "ILE", 3.5,  0,  0.0, true );
  set(particle::LEU, "LEU", 3.5,  0,  0.0, true );
  set(particle::LYS, "LYS", 3.5,  0, 10.4, false);
  set(particle::MET, "MET", 3.5,  0,  0.0, true );
  set(particle::PHE, "PHE", 3.5,  0,  0.0, true );
  set(particle::PRO, "PRO", 3.5,  0,  0.0, false);
  set(particle::SER, "SER", 3.5,  0,  0.0, false);
  set(particle::THR, "THR", 3.5,  0,  0.0, false);
  set(particle::TRP, "TRP", 3.5,  0,  0.0, true );
  set(particle::TYR, "TYR", 3.5, -1,  9.6, true );
  set(particle::VAL, "VAL", 3.5,  0,  0.0, false);
  set(particle::CTR, "CTR", 3.5, -1,  3.8, false);
  set(particle::NTR, "NTR", 3.5,  0,  7.5, false);
  set(particle::UNK, "UNK", 3.5,  0,  0.0, false);
  set(particle::NA,  "NA",  1.5, +1,  0.0, false);
  set(particle::K,   "K",   1.7, +1,  0.0, false);
  set(particle::CL,  "CL",  1.6, -1,  0.0, false);
  set(particle::BR,  "BR",  1.9, -1,  0.0, false);
  set(particle::I,   "I",   3.0, -1,  0.0, false);
  set(particle::GHOST,"GHOST",0,  0,  0.0, false);
};

