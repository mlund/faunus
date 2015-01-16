#ifndef FAU_JSON
#define FAU_JSON
#include <faunus/common.h>
#include <faunus/auxiliary.h>
#include <faunus/picojson.h>

namespace Faunus {
  /**
   * @brief Functions for handling JSON input files (www.json.org)
   *
   *
   * @details
   * The functions in this namespace are wrappers for the
   * [PicoJSON library](http://github.com/kazuho/picojson/):
   *
   * Example code:
   * @code
   * auto v = json::open("atoms.json");
   *  for (auto &atom : json::object("atomlist", v)) {
   *    string name = atom.first;
   *    double radius = json::value<double>(atom.second, "r", 2.0);
   *    double charge = json::value<double>(atom.second, "q", 0.0);
   *    bool hydrophobic = json::value<bool>(atom.second, "hydrophobic", false);
   *  }
   * @endcode
   * where the `atoms.json` file may look like the following:
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
    typedef picojson::object Tobj; //!< JSON object object (`std::map<string,value>`)
    typedef Tobj::value_type Tobjitem; //!< Type if item in Tobj (std::pair<string,value>)

    const Tobj emptyObject;

    /* Note in C++11 stringstream is "moveable" and can be returned from a function
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

    /**
     * @brief Load entire JSON file - "//" style comments are allowed.
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
      } else {
        throw std::runtime_error("JSON file " + file + " not loaded.");
      }
      return Tval();
    }

    /**
     * @brief Returns object associated with a keyword
     * @param key Keyword to look for - associated value must be an object
     * @param v The object in which to search for key
     */
    inline const Tobj& object(const std::string& key, const Tval& v) {
      if (v.is<Tobj>())
        if (v.get(key).is<Tobj>() )
          return v.get(key).get<Tobj>();
      return emptyObject;
    }

    /**
     * @brief Search object for key and return associated value
     * @param v Object to be seached
     * @param key Keyword to search for
     * @param fallback Default value if the keyword is not found
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
     * @brief CD (change dir) like functionality
     * 
     * Example:
     *
     * ~~~~
     *   auto a = json::open("myfile.json");
     *   auto b = json::cd ( a, {"general","system"} );
     *   double T = value( b, "temperature", 298.15 );
     * ~~~~
     *
     * where the json file may look like this:
     *
     * ~~~~
     *   {
     *     "general" : {
     *        "system" : { "temperature":341. },
     *        "other subsection" : { ... }
     *     }
     *   }
     * ~~~~
     */
    static Tval cd( const Tval &v, std::vector<string> dir) {
      if ( dir.empty() )
        return v;
      string s = dir.front();
      dir.erase( dir.begin() );
      if ( v.contains(s) )
        return cd( v.get(s), dir );
      std::cerr << "Error: Unable to find JSON section '" << s << "'" << endl;
      exit(1);
      return Tval(); // suppresses compiler warning
    }

  }//namespace
}//namespace
#endif

