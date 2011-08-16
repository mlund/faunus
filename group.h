#ifndef FAUNUS_GROUP_H
#define FAUNUS_GROUP_H

namespace Faunus {
  class point;
  class mcmove;
  class container;

  class group {
    public:
      int beg;                  //!< index of first particle
      int end;                  //!< index of last particle
      point cm;                 //!< mass center vector
      point cm_trial;           //!< mass center vector for trial position
      int len();                //!< number of particle in group
      vector<mcmove*> moves;    //!< pointers to move functions

      virtual double move(container&, int); //!< do a MC move
      virtual void undomove(container&);
  };

}//namespace

#endif

