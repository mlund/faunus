#ifndef FAUNUS_GROUP_H
#define FAUNUS_GROUP_H

#include <faunus/average.h>

namespace Faunus {
  class point;
  class mcmove;
  class particle;
  class container;

//  namespace Group {

    class group {
      public:
        group();
        string info();
        string name;
        int beg;                  //!< index of first particle
        int end;                  //!< index of last particle
        int len() const;          //!< number of particles in group.
        vector<mcmove*> moves;    //!< pointers to move functions

        virtual double charge(const vector<particle>&) const;      //!< Calculate total charge
        virtual point masscenter(const container&) const;          //!< Calculate mass center
        virtual point dipolemoment(const container&) const;        //!< Calculate dipole moment

        virtual void rotate(container&, const point&, double);     //!< Rotate around a vector
        virtual void translate(container&, const point&);          //!< Translate along a vector
        virtual void scale(container&, double);                    //!< Volume scaling
        virtual void undo(space&);                                 //!< Undo move operation
        virtual void accept(space&);

        bool operator==(const group&) const;
        group& operator+=(const group&);
        const group operator+(const group&) const;
    };

    class molecular  {
      protected:
        point cm_trial;           //!< mass center vector for trial position
      public:
        point cm;                 //!< mass center vector
        average<double> Q;        //!< average net charge
        average<double> mu;       //!< average dipole moment
    };

    class cigar : public molecular {
      public:
        point cap1, cap2, patch;
    };

    //  }//namespace

}//namespace

#endif

