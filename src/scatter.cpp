#include <faunus/scatter.h>
#include <faunus/species.h>

namespace Faunus {

  namespace Scatter {

    /*!
     * This will top up the tables with related form factors
     * not loaded from disk. I.e. if you loaded data for "HIS",
     * this will be duplicated for "HHIS" if not already loaded
     * and if it is defined in the Atoms class.
     */
    void FormFactorTable::addVariants() {
      for (auto &m : F) {            // loop over loaded F(q)
        string mname = atom[m.first].name;
        for (auto &a : atom.list)    // loop over species
          if (a.id!=m.first)         // if non-tabulated species
            if (a.name.find(mname)!=std::string::npos) // and it is a mutant
              F[a.id] = F[m.first];  // duplicate!
      }
    }

    /*!
     * Example of the file format, where the first line gives
     * the q-resolution in inverse angstroms:
     * \verbatim
     * l1: 0.015
     * l2: q     ALA   ARG
     * l3: 0.000 8.983 23.527
     * l4: 0.015 8.980 23.516
     * l5: ...
     * \endverbatim
     * \param filename Multi-column file with F(q) for different species
     * \param variants True (default) if atomic variants to loaded date should be generated
     */
    bool FormFactorTable::load(string filename, bool variants) {
      std::ifstream f(filename.c_str());
      if (f) {
        float dq;
        f >> dq; // read resolution from line 1

        // read atom names from line 2
        std::vector<string> namevec;
        string line, name;
        std::getline(f, line);
        std::stringstream ss(line);
        while (ss>>name)
          if (name!="q") {
            namevec.push_back(name);
            F[ atom[name].id ] = Ttable(dq,Ttable::XYDATA);
          }

        // read F(q) starting from line 3 until eof
        while (!f.eof()) {
          int rescnt=0;
          float _f, _q;
          std::getline(f, line);
          std::stringstream ss(line);
          ss >> _q;
          while (ss >> _f) {
            particle::Tid id = atom[ namevec.at(rescnt++) ].id;
            F[id](_q) = _f;
          }
        }
        if (variants)
          addVariants();
        return true;
      }
      return false;
    }

  }//namespace
}//namespace
