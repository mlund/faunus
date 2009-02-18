#ifndef FAU_METROPOLIS_H
#define FAU_METROPOLIS_H

#include "faunus/slump.h"
#include "faunus/group.h"

namespace Faunus {
  /*!
   * Polymorph ensemble class
   * \todo Maybe not really needed...
   */
  class ensemble {
    protected:
      slump slp;
    public:
      bool metropolis(double du)
      {
        if (du > 0)
          if ( slp.random_one()>exp(-du) )
            return false;
        return true;
      }
  };

  /*! \brief NVT ensemble
  */
  class canonical : public ensemble {};

  /*! \brief NPT ensemble
   *  \todo Silly name due to name collision...
   */
  class isobarical : public ensemble {};

  /*!
   * Class for handling Grand Canonical steps. This class can dynamically insert and
   * erase particles while at the same time update the groups in the systems.
   * \author Mikael Lund
   * \date Lund, 2009
   * \todo A lot.
   */
  class grandcanonical : public ensemble {
    private:
      int size();  //!< Calculates the number of particles in the specifies groups
      short findgroup(unsigned int i) {
        for (unsigned int n=0; n<g.size(); n++)
          if (g[n].find(i)==true) return n;
        return -1;
      }

    public:
      vector<group> g; //!< Groups with fluctuating species
      void addgroup(group &group) { g.push_back(group); }

      short findgroup(string name) {
        for (int k=0; k<g.size(); k++)
          if (g[k].name==name) return k;
        return -1;
      }
      void searchsalt(container &c, salt &salt) {
        particle ref = c.p[salt.beg];
        group w;
        w.beg=salt.beg;
        w.end=salt.beg-1;
        for (unsigned int i=salt.beg; i<=salt.end; i++)
          if (c.p[i].id==ref.id) 
            w.end++;
          else {
            w.name=c.atom[ref.id].name;
            g.push_back(w);
            w.beg=w.end=i;
            ref.id=c.p[i].id;
          }
        w.name=c.atom[ref.id].name;
        g.push_back(w);
      }

      //\param p - particle vector to insert into.
      //\param i - We will insert at this position - the rest will be moved forward.
      //\param p - Particle to insert
      //\param num - Number of particles to insert.
      bool insert(vector<particle> &p, unsigned int i, particle a, unsigned char num=1) {
        short n=findgroup(i);
        if (n>0) {
          for (unsigned short k=0; k<g.size(); k++)
            if (g[k].beg>i) {
              g[k].beg += num;
              g[k].end += num;
            }
          p.insert( p.begin()+i, num, a );
          g[n].end+=num;
          return true;
        }
        return false;
      }
      bool erase(vector<particle> &p, unsigned int i, unsigned char num=1) {
        short n=findgroup(i);
        if (n>0) {
          for (unsigned short k=0; k<g.size(); k++)
            if (g[k].beg>i) {
              g[k].beg -= num;
              g[k].end -= num;
            }
          p.erase( p.begin()+i, p.begin()+i+num );
          g[n].end -= num;
          return true;
        }
        return false;
      }

      string info() {
        std::ostringstream o;
        o << "Grand Canonical Master Class:\n"
          << "  Number of groups     = " << g.size() << endl
          << "  Number of particles  = " << size() << endl;
      }
  };

}//namespace
#endif
