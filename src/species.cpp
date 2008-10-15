#include "faunus/species.h"
namespace Faunus {
  //float species::getpka(particle::type id) { return d[id].pka; }
  particle species::get(particle::type id) const { return d[id].p; }
  particle species::get(string name) const { return d[id(name)].p; }

  void species::set(
      particle::type id, string name, float r, float z, float pka, bool hydr)
  {
    d[id].name = name;
    d[id].p.id=id;
    d[id].pka = pka;
    d[id].p.radius = r;
    d[id].p.charge = z;
    d[id].p.hydrophobic = hydr;
  }

  particle::type species::id(string name) const {
    for (int i=0; i<d.size(); i++)
      if (d[i].name==name)
        return particle::type(i);
    return particle::UNK;
  }

  species::species() {
    d.resize(particle::LAST);
    //  id             name   rad   Z   pKa  hydrp
    set(particle::ALA, "ALA", 3.5,  0,  0.0, true);
    set(particle::ARG, "ARG", 3.5,  0, 12.0, false);
    set(particle::ASN, "ASN", 3.5,  0,  0.0, false);
    set(particle::ASP, "ASP", 3.5, -1,  4.8, false);
    set(particle::CYS, "CYS", 3.5, -1, 10.8, false);
    set(particle::GLN, "GLN", 3.5,  0,  0.0, false);
    set(particle::GLU, "GLU", 3.5, -1,  4.8, false);
    set(particle::GLY, "GLY", 3.5,  0,  0.0, false);
    set(particle::HIS, "HIS", 3.5,  0,  6.3, false);
    set(particle::ILE, "ILE", 3.5,  0,  0.0, true );
    set(particle::LEU, "LEU", 3.5,  0,  0.0, true );
    set(particle::LYS, "LYS", 3.5,  0, 10.4, false);
    set(particle::MET, "MET", 3.5,  0,  0.0, true );
    set(particle::PHE, "PHE", 3.5,  0,  0.0, true );
    set(particle::PRO, "PRO", 3.5,  0,  0.0, true);
    set(particle::SER, "SER", 3.5,  0,  0.0, false);
    set(particle::THR, "THR", 3.5,  0,  0.0, false);
    set(particle::TRP, "TRP", 3.5,  0,  0.0, true );
    set(particle::TYR, "TYR", 3.5, -1,  9.6, false );
    set(particle::VAL, "VAL", 3.5,  0,  0.0, true);
    set(particle::CTR, "CTR", 3.5, -1,  3.8, false);
    set(particle::NTR, "NTR", 3.5,  0,  7.5, false);
    set(particle::UNK, "UNK", 3.5,  0,  0.0, false);
    set(particle::NA,  "NA",  1.8, +1,  0.0, false);
    set(particle::K,   "K",   1.5, +1,  0.0, false);
    set(particle::F,   "F",   1.4, -1,  0.0, false);
    set(particle::CL,  "CL",  1.7, -1,  0.0, false);
    set(particle::BR,  "BR",  1.9, -1,  0.0, false);
    set(particle::I,   "I",   2.0, -1,  0.0, false);
    set(particle::SO4, "SO4", 1.6, -2,  0.0, false);
    set(particle::LA,  "LA",  2.0, +3,  0.0, false);
    set(particle::GHOST,"GHOST",0,  0,  0.0, false);
    set(particle::HYDROPHOBIC,"HYDR",3.5,0,0, true);
  }

  string species::particleinfo(particle::type id) const {
    std::ostringstream o;
    o << "# Particle specs (name,z,r,pka): "
      << d[id].name << " " << d[id].p.charge << " " << d[id].p.radius
      << " " << d[id].pka << std::endl;
    return o.str();
  }
}
