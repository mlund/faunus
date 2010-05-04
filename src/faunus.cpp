#include <faunus/faunus.h>

namespace Faunus {
  string faunus_splash() {
    std::ostringstream o;
    o << "-----------------------------------------------------------------------\n"
      << "  Welcome to FAUNUS - A Framework for Molecule Simulation.\n"
      << "  Copyright (C) 2002-2010 Mikael Lund\n"
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
      << "    Mikael Lund, Bjorn Persson, Martin Trulsson,\n"
      << "    Ondrej Marsalek, Christophe Labbez, Andre Teixeira,\n"
      << "    Anil Kurut.\n"
      << "\n"
      << "  Reference:\n"
      << "    Source Code Biol. Med. (2008) 3:1\n"
      << "\n"
      << "  Library build details:\n"
      << "    Compiled on " << __DATE__ << " " << __TIME__ 
#ifdef __VERSION__
      << " using " << __VERSION__
#endif
#ifdef __SVN_REV__
      << "\n    SVN revision: " << __SVN_REV__ << "."
#endif
      << "\n"
      << "----------------------------------------------------------------------\n" << std::endl;
    return o.str();
  }
}//namespace
