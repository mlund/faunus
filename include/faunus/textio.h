#ifndef FAUNUS_TEXTIO_H
#define FAUNUS_TEXTIO_H

#include <faunus/common.h>

namespace Faunus {
  /*!
   * \brief Namespace for text related operations: formatting, special characters etc.
   */
  namespace textio {
    enum indentlevel {TITLE=0,SUB=2,SUBSUB=4};

    const string angstrom="\u00c5";   //!< Angstrom symbol
    const string _angstrom=" \u00c5"; //!< Angstrom symbol with space in front
    const string cubed="\u00b3";      //!< Superscript 3
    const string degrees="\u00b0";    //!< Degrees
    const string epsilon="\u03f5";    //!< Greek epsilon
    const string gamma="\u0263";      //!< Greek gamma
    const string kT=" kT";            //!< kT (energy) with space in front
    const string mu="\u03bc";         //!< Greek mu
    const string percent="\ufe6a";    //!< Percent sign
    const string pm="\u00b1";         //!< Plus minus sign
    const string rho="\u03C1";        //!< Greek rho
    const string rootof="\u221a";     //!< Square root sign
    const string squared="\u00b2";    //!< Superscript 2
    const string sigma="\u03c3";      //!< Greek sigma
    const string theta="\u03b8";      //!< Greek theta

    string splash();                              //!< Show Faunus welcome text, version etc.
    string bracket(const string&);                //!< Put angular brackets around string
    string header(const string&);                 //!< Print header for info() functions
    string indent(indentlevel);                   //!< Indent text
    string pad(indentlevel, char, const string&); //!< Pad and indent text
    string trim(string);                          //!< Remove white space from string
  }//namespace
}//namespace
#endif

