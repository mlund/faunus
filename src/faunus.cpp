#include "faunus/faunus.h"
#include <faunus/common.h>

namespace Faunus {
  
  string faunus_splash() {
    std::ostringstream o;
    o << "-----------------------------------------------------------------------\n"
      << "  Welcome to FAUNUS - A Framework for Molecule Simulation.\n"
      << "  Copyright (C) 2002-2011 Mikael Lund\n"
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
      << "  Programmed by:\n"
      << "    Mikael Lund, Bj\u00F6rn Persson, Martin Trulsson,\n"
      << "    Ond\u0159ej Mar\u0161\u00E1lek, Christophe Labbez, Andre Teixeira,\n"
      << "    An\u0131l Kurut, Chris Evers, Robert V\u00E1cha\n"
      << "\n"
      << "  Reference:\n"
      << "    Source Code Biol. Med. (2008) 3:1\n"
      << "\n"
      << "  Library build details:\n"
      << "    Compiled on " << __DATE__ << " " << __TIME__
#ifdef __VERSION__
      << "\n    " << __VERSION__
#endif
#ifdef __SVN_REV__
      << "\n    SVN revision: " << __SVN_REV__ << "."
#endif
      << "\n"
      << "----------------------------------------------------------------------\n" << std::endl;
    return o.str();
  }

  string header(const string &s) {
    int l=s.size()+4;
    return "\n" + std::string(l,'=') + "\n  " + s + "  \n" + string(l,'=') + "\n";
  }

  string indent(indentlevel level) { return std::string(level, ' '); }

  string pad(const string &s, char width, indentlevel level) {
    std::ostringstream o;
    o << indent(level) << std::left << std::setw(width) << s;
    return o.str();
  }

}//namespace
