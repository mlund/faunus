/***************************************************************************
  Faunus -- A Framework for Molecular Modelling 
  Copyright (C) 2002-2009 Mikael Lund 

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or 
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 ***************************************************************************/
#ifndef FAUNUS_H
#define FAUNUS_H

#ifndef SWIG
#include <faunus/common.h>
#endif

namespace Faunus {
  enum indentlevel {TITLE=0,SUB=2,SUBSUB=4};
  string faunus_splash();
  string header(const string&);
  string indent(indentlevel);
  string pad(const string&, char, indentlevel);

  class textOutput {
    private:
      unsigned short width;
    public:
      template<class T> void printValue(string s, T v) {
        std::cout << s << " = " << v << std::endl;
      }
  };
}//namespace

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/geometry.h>
#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/energy.h>
#include <faunus/potentials.h>
#include <faunus/average.h>
#include <faunus/space.h>
#include <faunus/move.h>
#include <faunus/mcloop.h>
#include <faunus/group.h>
#include <faunus/io.h>
#include <faunus/drift.h>
#include <faunus/xytable.h>
#endif

#endif
