#ifndef FAUNUS_TEXTIO_H
#define FAUNUS_TEXTIO_H

#include <faunus/common.h>

namespace Faunus {
  namespace textio {
    enum indentlevel {TITLE=0,SUB=2,SUBSUB=4};

    const string percent="\ufe6a";
    const string angstrom="\u00c5";
    const string _angstrom=" \u00c5";
    const string squared="\u00b2";
    const string cubed="\u00b3";
    const string rootof="\u221a";
    const string sigma="\u03c3";
    const string kT=" kT";

    string splash();
    string bracket(const string&);
    string header(const string&);
    string indent(indentlevel);
    string pad(indentlevel, char, const string&);
  }//namespace
}//namespace
#endif

