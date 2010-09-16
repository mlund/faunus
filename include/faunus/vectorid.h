#ifndef _VECTORID_
#define _VECTORID_
#include <vector>
#include <string>
namespace Faunus {
  /*!
   * \brief Vector class template with access by string id
   *
   * This class is similar to STL vector accept that it
   * provides element access via strings (default). The class
   * is, however, more versatile than that. Instead of string
   * id's one could use any id type that has the "=" operator
   * - see example.
   *
   * Examples:\n
   * \code
   * vectorid<int> v;
   * v["x"]=100;
   * v["y"]=1;
   * cout << v["x"]-v["y"] << " " << v.size(); //-> 99 2
   * \endcode
   *
   * \code
   * vectorid<string, particle*> v;
   * particle a;
   * v[&a]="particleA";
   * \endcode
   *
   * \author Mikael Lund
   * \date Malmo, 2010
   */
  template<class T, class Tid=std::string> class vectorid {
  private:
    std::vector<Tid> idv;
    std::vector<T> v;
  public:
    int size() { return v.size(); }              //!< vector length

    T & at(int i) { return v.at(i); }            //!< access by index 

    Tid & id(int i) { return idv.at(i); }//!< id of i'th element

    T & operator[] (Tid id) {            //!< access by string id
      int i=find(id);
      if (i>=0)
        return v.at(i);
      T neu;
      v.push_back(neu);
      idv.push_back(id);
      return v.back();
    }

    int find(Tid id) {
      for (int i=0; i<idv.size(); i++)
        if (id==idv[i])
          return i;
      return -1;
    }
  };
}//namespace
#endif
