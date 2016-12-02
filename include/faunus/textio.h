#ifndef FAUNUS_TEXTIO_H
#define FAUNUS_TEXTIO_H

#include <faunus/common.h>
#include <iterator>

namespace Faunus
{
  /**
   * @brief Namespace for text related operations: formatting, special characters etc.
   *
   * Unless using the macro definition AVOID_UNICODE, Faunus will use UTF-16 unicode
   * output to print fancy output with mathematical symbols, greek letters etc.
   */
  namespace textio
  {
    enum indentlevel { TITLE = 0, SUB = 2, SUBSUB = 4 };
    static std::string prefix = "";   //!< Unique prefix for current job. Use for file I/O.

#ifdef AVOID_UNICODE
    const string angstrom="AA";
    const string _angstrom=" AA";
    const string beta="B";
    const string cubed="^3";
    const string cuberoot="3root";
    const string degrees="deg";
    const string epsilon="eps";
    const string epsilon_m="eps_m";
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
    const string angstrom = "\u00c5";   //!< Angstrom symbol
    const string _angstrom = " \u00c5"; //!< Angstrom symbol with space in front
    const string beta = "\u03b2";       //!< Greek beta
    const string cubed = "\u00b3";      //!< Superscript 3
    const string cuberoot = "\u221b";   //!< Cubic root
    const string degrees = "\u00b0";    //!< Degrees
    const string epsilon = "\u03f5";    //!< Greek epsilon
    const string epsilon_m = "\u03b5";  //!< Greek epsilon (minuscule)
    const string gamma = "\u0263";      //!< Greek gamma
    const string Gamma = "\u0393";      //!< Greek capital gamma
    const string kappa = "\u03ba";      //!< Greek kappa
    const string kT = " kT";            //!< kT (energy) with space in front
    const string mu = "\u03bc";         //!< Greek mu
    const string partial = "\u2202";    //!< Partial derivative
    const string percent = "\ufe6a";    //!< Percent sign
    const string pm = "\u00b1";         //!< Plus minus sign
    const string rho = "\u03C1";        //!< Greek rho
    const string rootof = "\u221a";     //!< Square root sign
    const string squared = "\u00b2";    //!< Superscript 2
    const string sigma = "\u03c3";      //!< Greek sigma
    const string superminus = "\u207b"; //!< Superscript minus (-)
    const string subr = "\u1D63";       //!< Subscript "r"
    const string theta = "\u03b8";      //!< Greek theta
#endif

    /** @brief Remove white space from string */
    inline std::string trim( std::string s )
    {
        s.erase(std::remove_if(s.begin(), s.end(), isspace), s.end());
        return s;
    }

    /** @brief  Put angular brackets around string */
    inline string bracket( const string &s )
    {
#ifdef AVOID_UNICODE
        return "<"+s+">";
#else
        return "\u27e8" + s + "\u27e9";
#endif
    }

    /** @brief Print header for info() functions */
    inline string header( const string &s )
    {
        int l = s.size() + 2;
        return "\n " + std::string(l, '.') + "\n  " + s + "  \n " + string(l, '*') + "\n";
    }

    /** @brief Indent text */
    inline string indent( indentlevel level )
    {
        return std::string(level, ' ');
    }

    /** @brief Pad and indent text */
    inline string pad( indentlevel level, char width, const string &s )
    {
        std::stringstream o;
        o << indent(level) << std::left << std::setw(width) << s;
        return o.str();
    }

    /** @brief Count number of white-space separated words in a string */
    inline size_t numWords( const std::string &s )
    {
        return std::distance(std::istream_iterator<std::string>(
            std::istringstream(s) >> std::ws), std::istream_iterator<std::string>());
    }

    /**
     * @brief Convert whitespace separated words into vector of given type
     *
     * Example:
     *
     * ~~~~
     * auto v = textio::words2vec<double>( "0.2 1 100" );
     * for (auto i : v)
     *   cout << 2*i << " "; // -> 0.4 2 200
     * ~~~~
     */

    template<class T>
    std::vector<T> words2vec( const std::string &w )
    {
        std::vector<T> v(numWords(w));
        std::stringstream s(w);
        size_t i = 0;
        while ( i < v.size())
        {
            s >> v[i++];
        }
        return v;
    }

    /** @brief Convert string to lower case */
    inline std::string lowercase( std::string s )
    {
        std::transform(s.begin(), s.end(), s.begin(), ::tolower);
        return s;
    }

    /**
     * @brief Show Faunus welcome text, version etc.
     * @note See http://patorjk.com/software/taag for ASCII art generation
     */
    inline string splash()
    {
        std::ostringstream o;
        o << std::string(71, '.') << endl
          << "  Welcome to FAUNUS - A Framework for Molecule Simulation.\n"
          << "  Copyright (C) 2002-2016 Mikael Lund\n"
          << "\n"
          << "  This program is free software; you can redistribute it and/or modify\n"
          << "  it under the terms of the GNU General Public License as published by\n"
          << "  the Free Software Foundation; either version 2 of the License, or\n"
          << "  (at your option) any later version.\n"
          << "\n"
          << "   ___________                               _________\n"
          << "   \\_   _____/_____    __ __   ____   __ __ /   _____/\n"
          << "    |    __)  \\__  \\  |  |  \\ /    \\ |  |  \\\\_____  \\\n"
          << "    |     \\    / __ \\_|  |  /|   |  \\|  |  //        \\\n"
          << "    \\___  /   (____  /|____/ |___|  /|____//_______  /\n"
          << "        \\/         \\/             \\/               \\/\n"
          << "\n"
          << "  Developed by:\n"
          << "    Mikael Lund, Bj\u00F6rn Persson, Martin Trulsson,\n"
          << "    Ond\u0159ej Mar\u0161\u00E1lek, Christophe Labbez, Andre Teixeira,\n"
          << "    An\u0131l Kurut, Chris Evers, Robert V\u00E1cha, Axel Thuresson,\n"
          << "    Bj\u00F6rn Stenqvist, Jo\u00E3o Henriques, Alexei Abrikossov,\n"
          << "    Giulio Tesei, Luk\u00E1\u0161 Suken\u00edk, Niels Kouwenhoven\n"
          << "\n"
          << "  References:\n"
          << "    http://doi.org/nvn\n"
          << "    http://doi.org/dfqgch\n"
          << "\n"
          << "  Library build details:\n"
          << "    Compiled on " << __DATE__ << " " << __TIME__
          #ifdef GIT_COMMIT_HASH
          << "\n    Git commit " << GIT_COMMIT_HASH
          #endif
          << "\n" << string(71, '*') << endl;
        return o.str();
    }

  }//end of textio namespace
}// end of Faunus namespace
#endif

