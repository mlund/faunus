/***************************************************************************
  Faunus -- A Framework for Molecular Modelling 
  Copyright (C) 2002-2012 Mikael Lund 

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
#include <faunus/json.hpp> // lohmann modern json
#include <faunus/common.h>
#include <faunus/textio.h>
#include <faunus/point.h>
#include <faunus/geometry.h>
#include <faunus/json.h>
#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/energy.h>
#include <faunus/potentials.h>
#include <faunus/multipole.h>
#include <faunus/externalpotential.h>
#include <faunus/average.h>
#include <faunus/space.h>
#include <faunus/move.h>
#include <faunus/titrate.h>
#include <faunus/mcloop.h>
#include <faunus/group.h>
#include <faunus/io.h>
#include <faunus/tabulate.h>
#include <faunus/analysis.h>
#include <faunus/scatter.h>
#include <faunus/spherocylinder.h>

#endif

#endif
