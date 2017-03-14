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
      std::ifstream f;
      f.exceptions( std::ifstream::failbit );
      try {
          f.open( file );
          if ( f )
              js << f;
      }
      catch(std::exception& e)
      {
          throw std::runtime_error("Error loading json file " + file + ": " + e.what());
      }
      return js;
  }
}//namespace
#endif

