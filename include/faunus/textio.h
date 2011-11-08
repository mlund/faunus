#ifndef FAUNUS_TEXTIO_H
#define FAUNUS_TEXTIO_H

#include <faunus/common.h>

namespace Faunus {
  namespace textio {
    enum indentlevel {TITLE=0,SUB=2,SUBSUB=4};

    const string angstrom="\u00c5";
    const string _angstrom=" \u00c5";
    const string cubed="\u00b3";
    const string degrees="\u00b0";
    const string epsilon="\u03f5";
    const string gamma="\u0263";
    const string kT=" kT";
    const string percent="\ufe6a";
    const string pm="\u00b1";
    const string rho="\u03C1";
    const string rootof="\u221a";
    const string squared="\u00b2";
    const string sigma="\u03c3";
    const string theta="\u03b8";

    string splash();
    string bracket(const string&);
    string header(const string&);
    string indent(indentlevel);
    string pad(indentlevel, char, const string&);
    string trim(string); //!< Remove white space from string
  }//namespace
}//namespace
#endif

