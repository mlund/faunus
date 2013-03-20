#ifndef FAUNUS_TEXTIO_H
#define FAUNUS_TEXTIO_H

#ifndef SWIG
#include <faunus/common.h>
#endif

namespace Faunus {
  /*!
   * \brief Namespace for text related operations: formatting, special characters etc.
   *
   * Unless using the macro definition AVOID_UNICODE, Faunus will use UTF-16 unicode
   * output to print fancy output with mathematical symbols, greek letters etc.
   */
  namespace textio {
    enum indentlevel {TITLE=0,SUB=2,SUBSUB=4};

#ifdef AVOID_UNICODE
    const string angstrom="AA";
    const string _angstrom=" AA";
    const string cubed="^3";
    const string cuberoot="3root";
    const string degrees="deg";
    const string epsilon="eps";
    const string gamma="gamma";
    const string Gamma="Gamma";
    const string kappa="k";
    const string kT=" kT";
    const string mu="mu";
    const string percent="%";
    const string partial="d";
    const string pm="+-";
    const string rho="rho";
    const string rootof="sqrt";
    const string squared="^2";
    const string sigma="sigma";
    const string superminus="-";
    const string subr="r";
    const string theta="theta";
#else
    const string angstrom="\u00c5";   //!< Angstrom symbol
    const string _angstrom=" \u00c5"; //!< Angstrom symbol with space in front
    const string cubed="\u00b3";      //!< Superscript 3
    const string cuberoot="\u221b";   //!< Cubic root
    const string degrees="\u00b0";    //!< Degrees
    const string epsilon="\u03f5";    //!< Greek epsilon
    const string gamma="\u0263";      //!< Greek gamma
    const string Gamma="\u0393";      //!< Greek capital gamma
    const string kappa="\u03ba";      //!< Greek kappa
    const string kT=" kT";            //!< kT (energy) with space in front
    const string mu="\u03bc";         //!< Greek mu
    const string partial="\u2202";    //!< Partial derivative
    const string percent="\ufe6a";    //!< Percent sign
    const string pm="\u00b1";         //!< Plus minus sign
    const string rho="\u03C1";        //!< Greek rho
    const string rootof="\u221a";     //!< Square root sign
    const string squared="\u00b2";    //!< Superscript 2
    const string sigma="\u03c3";      //!< Greek sigma
    const string superminus="\u207b"; //!< Superscript minus (-)
    const string subr="\u1D63";       //!< Subscript "r"
    const string theta="\u03b8";      //!< Greek theta
#endif

    string splash();                              //!< Show Faunus welcome text, version etc.
    string bracket(const string&);                //!< Put angular brackets around string
    string header(const string&);                 //!< Print header for info() functions
    string indent(indentlevel);                   //!< Indent text
    string pad(indentlevel, char, const string&); //!< Pad and indent text
    string trim(string);                          //!< Remove white space from string

    extern std::string prefix;                    //!< Unique prefix for current job. Use for file I/O.
#ifndef SWIG
    extern std::ostream &fcout;                   //!< Alias for standard output (can be redirected)
    extern std::ostream &fcerr;                   //!< Alias for standard error (can be redirected)
#endif

  }//end of textio namespace
}// end of Faunus namespace
#endif

