#ifndef FAU_JSON
#define FAU_JSON
#include <faunus/common.h>
#include <faunus/json.hpp>

namespace Faunus
{

  typedef nlohmann::json Tmjson;

  /** @brief Merge two json objects */
  inline Tmjson merge( const Tmjson &a, const Tmjson &b )
  {
      Tmjson result = a.flatten();
      Tmjson tmp = b.flatten();
      for ( auto it = tmp.begin(); it != tmp.end(); ++it )
          result[it.key()] = it.value();
      return result.unflatten();
  }

  inline Tmjson openjson( const string &file )
  {
      Tmjson js;
      std::ifstream f (file );
      if ( f ) {
          try {
              js << f;
          }
          catch(std::exception& e)
          {
              throw std::runtime_error("Syntax error in JSON file " + file + ": " + e.what());
          }
      }
      else
          throw std::runtime_error("Cannot find or read JSON file " + file);
      return js;
  }
}//namespace
#endif

