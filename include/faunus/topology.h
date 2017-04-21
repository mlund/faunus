
namespace Faunus
{

  class System
  {
      Space p, trial;

      vector <Group> g;

  }

  class topology
  {
  protected:

      class atomtype
      {
      public:
          string name;
          unsigned short id;
          double charge;
          double sigma;
          double epsilon;
          double mw;
          double dp;
      };

      class bond
      {
      public:
          unsigned short id;
          int i;
          int j;
      };

  public:

      vector <atomtype> atomlist;
      vector <bond> bondlist;

  };
}//namespace
