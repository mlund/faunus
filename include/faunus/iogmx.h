//beginning of a gromacs topology class...

namespace Faunus {
  class gmx_atomtype {
    public:
      void getheader(ostringstream&);
      void getio(ostringstream&);
      void set(const particle&);
      string name;
      double charge;
      double sigma;
      double epsilon;
      double mw;
      string residue;
  };

  void gmx_atomtype::getio(ostringstream &o) {
    o << name << " " << "...";
  }

  void gmx_atomtype::getheader(ostringstream &o) {
    o << "[ atomtype ]" << endl
      << "; name  bond_type    mass    charge   ptype          sigma      epsilon" << endl;
  }

  class gmx_atom : public gmx_atomtype {
    public:
      void getio(ostringstream&);
      void set(const particle&);
      double x,y,z;
  };

  class gmx_molecule {
    public:
      void getio(ostringstream&);
      void set(const vector<particle>&, group&, string);
      string name;
      vector<gmx_atom> atom;
  };

  class gmx_topology {
    public:
      vector<gmx_atomtype> atomtype;
      vector<gmx_molecule> molecule;
      float version;
      int LJmixingrule;
      write(string);
  };
}//namespace
