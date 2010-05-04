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

#include "faunus/common.h"
#include "faunus/slump.h"
#include "faunus/physconst.h"
#include "faunus/histogram.h"
#include "faunus/analysis.h"
#include "faunus/particles.h"
#include "faunus/io.h"
#include "faunus/inputfile.h"
#include "faunus/container.h"
#include "faunus/countdown.h"
#include "faunus/mcloop.h"
#include "faunus/energy/base.h"
#include "faunus/widom.h"
#include "faunus/moves/translational.h"
#include "faunus/moves/rotational.h"
#include "faunus/moves/charge.h"
#include "faunus/moves/volume.h"
#include "faunus/moves/miscmove.h"
#include "faunus/moves/rosenbluth.h"
#include "faunus/moves/clustermove.h"
#include "faunus/moves/crankshaft.h"
#include "faunus/moves/saltbath.h"

namespace Faunus {
  string faunus_splash();

  class textOutput {
    private:
      unsigned short width;
    public:
      void printTitle(string);
      template<class T> void printValue(string s, T v) {
        std::cout << s << " = " << v << std::endl;
      }
  };
}//namespace
#endif
