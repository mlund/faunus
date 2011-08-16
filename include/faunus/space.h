#include <vector>
namespace Faunus {
  class particle;
  class group;

  class bondlist {
    public:
      bool add(int i, int j);

  };

  class space {
    public:
    vector<particle> p;
    vector<particle> trial;
    vector<group*> group;

    bool insert(particle, int);
    bool insert(const vector<particle> &, int);
    bool delete(int, int);

  };
} //namespace
