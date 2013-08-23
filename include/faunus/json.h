#ifndef FAU_JSON
#define FAU_JSON
#include <faunus/common.h>
#include <faunus/auxiliary.h>
#include <faunus/species.h>
#include <faunus/picojson.h>

namespace Faunus {
  /*!
   * \brief Function for handling JSON input files (www.json.org)
   *
   *
   * \details
   * The functions in this namespace are wrappers for the PicoJSON library:
   * http://github.com/kazuho/picojson
   *
   * Example code:
   * \code
   * auto v = json::open("atoms.json");
   *  for (auto &atom : json::object("atomlist", v)) {
   *    string name = atom.first;
   *    double radius = json::value<double>(atom.second, "r", 2.0);
   *    double charge = json::value<double>(atom.second, "q", 0.0);
   *    bool hydrophobic = json::value<bool>(atom.second, "hydrophobic", false);
   *  }
   *  \endcode
   * where the \c atoms.json file may look like the following:
   \code{.js}
   { "atomlist" :
     {
       "Na" : { "q" : 1.0, "r" : 3.0 },
       "Cl" : { "q" : -1.0, "r" : 1.7 },
      "LEU" : { "r" : 3.1, "hydrophobic" : true }
     }
   } 
   \endcode
   */
  namespace json {
    typedef picojson::value Tval;  //!< JSOB value object
    typedef picojson::object Tobj; //!< JSON object object (merely a std::map)

    const Tobj emptyObject;

    /* Note in C++11 stringstream is "moveable" and can be return from a function
     * Requires gcc 4.7?
     */
    inline void stripComments(std::ifstream &in, std::stringstream &out, string id="//") {
      string line;
      while (std::getline(in, line)) {
        auto pos=line.find_first_of(id); // search for "id"
        if (pos!=string::npos)           // if present
          line.erase(pos);               // delete rest of line
        if (!line.empty())
          out << line << endl;           // send to stream
      }
    }

    /*!
     * \brief Load entire JSON file - "//" style comments are allowed.
     */
    inline Tval open(const std::string& file) {
      std::ifstream in(file.c_str());
      if (in) {
        Tval v;
        std::stringstream stripped;
        stripComments(in,stripped);
        stripped >> v;
        if (v.is<Tobj>())
          return v;
      }
      std::cerr << "Error loading JSON file, '" << file << "', or invalid content.\n";
      return Tval();
    }

    /*!
     * \brief Returns object associated with a keyword
     * \param key Keyword to look for - associated value must be an object
     * \param v The object in which to search for key
     */
    inline const Tobj& object(const std::string& key, const Tval& v) {
      if (v.is<Tobj>())
        if (v.get(key).is<Tobj>() )
          return v.get(key).get<Tobj>();
      return emptyObject;
    }

    /*!
     * \brief Search object for key and return associated value
     * \param v Object to be seached
     * \param key Keyword to search for
     * \param fallback Default value if the keyword is not found
     */
    template<typename T>
      static T value(const Tval& v, const std::string& key, const T& fallback) {
        if (v.is<Tobj>())
          if (v.contains(key))
            if (v.get(key).is<T>())
              return v.get(key).get<T>();
        return fallback;
      }

    /**
     * @brief Loads a JSON file, reads atom pair properties and returns a vector map
     *
     * Example:
     * ~~~~
     * auto map = atomPairMap("input.json", "pairproperties", "nemorep");
     * for (auto &m : map)
     *   cout << m.second.transpose() << endl; // -> 12 23 0.2 -2 3 4 5 ...
     * ~~~~
     * where the input json file could look like this:
     * ~~~~
     * {
     *   "pairproperties" : {
     *      "OW OW"  : { "nemorep":"12. 23. 0.2  -2   3   4   5" },
     *      "HW HW"  : { "nemorep":"-2. 23. 0.2   2  99   4  -5 " },
     *      "HW OW"  : { "nemorep":"112. 23. 0.2 129 391 238  23" }
     *   }
     * }
     * ~~~~
     */
    inline std::map<opair<int>,Eigen::VectorXd>
      atomPairMap(const string &file, const string &section, const string &key) {
        assert(!section.empty() && !key.empty());
        typedef Eigen::VectorXd Tvec;
        typedef opair<int> Tpair;
        std::map<Tpair,Tvec> map;
        string atom1, atom2;
        auto j=json::open(file);
        for (auto &a : json::object(section, j)) {
          std::istringstream is(a.first);
          is >> atom1 >> atom2;
          Tpair pair( atom[atom1].id, atom[atom2].id );
          string str = json::value<string>(a.second, key,"");
          std::istringstream is2(str), tmp(str);
          int size = std::distance(std::istream_iterator<string>(tmp), std::istream_iterator<string>());
          Tvec v(size);
          for (int i=0; i<size; i++)
            is2 >> v[i];
          map[pair] = v;
        }
        return map;
      }
  }//namespace
}//namespace
#endif

