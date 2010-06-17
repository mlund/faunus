#ifndef _VECTORID_
#define _VECTORID_
#include <vector>
#include <string>
namespace Faunus {
  /*!
   * \brief Vector class template with acess by string id
   *
   * This class is similar to STL vector accept that it
   * provides element access via strings.
   *
   * Example:\n
   * \code
   * vectorid<int> v;
   * v["x"]=100;
   * v["y"]=1;
   * cout << v["x"]-v["y"] << " " << v.size(); //-> 99 2
   * \endcode
   *
   * \author Mikael Lund
   * \date Malmo, 2010
   */
  template<class T> class vectorid {
  private:
    std::vector<std::string> idv;
    std::vector<T> v;
  public:
    int size() { return v.size(); }              //!< vector length
    T & at(int i) { return v.at(i); }            //!< access by index 
    std::string & id(int i) { return idv.at(i); }//!< id of i'th element
    T & operator[] (std::string id) {            //!< access by string id
      for (int i=0; i<idv.size(); i++)
        if (id==idv[i])
          return v.at(i);
      T neu;
      v.push_back(neu);
      idv.push_back(id);
      return v.back();
    }
  };
}//namespace
#endif
