#include <faunus/textio.h>

namespace Faunus {
  /*!
   * Please direct all output to stdout here. By default this is exactly the same
   * as using std::cout but by using this alias it is possible to redirect all output as
   * needed in for example MPI code.
   */
  std::ostream& textio::fcout = std::cout;

  /*!
   * As textio::fcout but for standard error.
   */
  std::ostream& textio::fcerr = std::cerr;

  string textio::trim(string s) {
    s.erase( std::remove_if(s.begin(), s.end(), isspace), s.end());
    return s;
  }

  string textio::bracket(const string &s) {
#ifdef AVOID_UNICODE
    return "<"+s+">";
#else
    return "\u27e8"+s+"\u27e9";
#endif
  }

  string textio::header(const string &s) {
    int l=s.size()+2;
    return "\n " + std::string(l,'.') + "\n  " + s + "  \n " + string(l,'*') + "\n";
  }

  string textio::indent(indentlevel level) {
    return std::string(level, ' ');
  }

  string textio::pad(indentlevel level, char width, const string &s) {
    std::stringstream o;
    o << indent(level) << std::left << std::setw(width) << s;
    return o.str();
  }

  /*!
   * \note See http://patorjk.com/software/taag for ASCII art generation
   */
  string textio::splash() {
    std::ostringstream o;
    o << std::string(71,'.') << endl
      << "  Welcome to FAUNUS - A Framework for Molecule Simulation.\n"
      << "  Copyright (C) 2002-2012 Mikael Lund\n"
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
      << "    An\u0131l Kurut, Chris Evers, Robert V\u00E1cha\n"
      << "\n"
      << "  Reference:\n"
      << "    Source Code Biol. Med. (2008) 3:1\n"
      << "    doi:10.1186/1751-0473-3-1\n"
      << "\n"
      << "  Library build details:\n"
      << "    Compiled on " << __DATE__ << " " << __TIME__
#ifdef __VERSION__
      << "\n    " << __VERSION__
#endif
#ifdef __SVN_REV__
      << "\n    SVN revision: " << __SVN_REV__ << "."
#endif
      << "\n" << string(71,'*') << endl;
    return o.str();
  }

}//namespace
