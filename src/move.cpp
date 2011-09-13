#include <faunus/move.h>

namespace Faunus {

  namespace Move {

    movebase::movebase(string pfx) {
      prefix=pfx;
      cnt=cnt_accepted=0;
      dusum=0;
      iw=22;
    }

    //void move::unittest(unittest&) {
    //}

    void movebase::pad(std::ostringstream& o) {
      o << "#   " << setw(iw) << std::left;
    }

    string movebase::info() {
      std::ostringstream o;
      o << "# " << title << endl;
      pad(o); o << "More information:" << cite << endl;
      pad(o); o << "Runfraction" << runfraction << endl;
      if (cnt>0) {
        pad(o); o << "Number of trials" << cnt << endl;
        pad(o); o << "Acceptance" << double(cnt_accepted)/cnt << endl;
        pad(o); o << "Total energy change" << dusum << " kT" << endl;
      }
      return o.str();
    }

    // TRANSLATE

    translate::translate(std::string pfx) : movebase(pfx) {
      title="Molecular Translation";
    }

    void translate::trialmove() {
    }

    void translate::acceptmove() {
      cnt_accepted++;
    }

    void translate::rejectmove() {
    }

    double translate::energychange() {return 0;}

    double translate::move() {
      return 0;
    }

    string translate::info() {
      std::ostringstream o;
      o << movebase::info();
      pad(o); o << "Displacement vector" << dp << endl; 
      return o.str();
    }


  }//namespace
}//namespace
