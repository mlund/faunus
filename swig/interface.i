%module (docstring="A Molecular Framework for Moleculer Simulation") pyfaunus
%{
#include "faunus/common.h"
#include "faunus/point.h"
#include "faunus/average.h"
#include "faunus/inputfile.h"
#include "faunus/textio.h"
#include "faunus/potentials.h"
#include "faunus/geometry.h"
#include "faunus/slump.h"
#include "faunus/species.h"
#include "faunus/space.h"
#include "faunus/group.h"
#include "faunus/energy.h"
#include "faunus/move.h"
#include "faunus/drift.h"
#include "faunus/analysis.h"
#include "faunus/physconst.h"
#include "faunus/titrate.h"
#include "faunus/mcloop.h"
#include "faunus/histogram.h"
#include "faunus/io.h"
//#include "faunus/mpi.h"
%}

/* -----------------------------------
 *          Do some renaming
 * ----------------------------------- */
%rename(assign) Faunus::Point::operator=;
%rename(assign) Faunus::PointParticle::operator=;
%rename(__int__) Faunus::Average<int>::operator int;
%rename(__float__) Faunus::Average<double>::operator double;
%rename(__float__) Faunus::Average<float>::operator float;
%rename(__getitem__) Faunus::AtomTypes::operator[];

%ignore *::operator<<(std::istream&);
%ignore *::operator<<(std::ostream&, const Faunus::Point&);
%ignore Faunus::Average::operator=(T);
%ignore Faunus::BlockCorrelation::operator[] (unsigned int);
%ignore Faunus::Move::AtomTracker::operator[] (particle::Tid);

%include "std_vector.i"
%include "std_string.i"
//%include "std_iostream.i"
//%include "std_streambuf.i"
//%include "std_common.i"
//%include "std_map.i"
%include "faunus/common.h"
%include "faunus/point.h"
%include "faunus/average.h"
%include "faunus/inputfile.h"
%include "faunus/textio.h"
%include "faunus/potentials.h"
%include "faunus/geometry.h"
%include "faunus/slump.h"
%include "faunus/species.h"
%include "faunus/space.h"
%include "faunus/group.h"
%include "faunus/energy.h"
%include "faunus/move.h"
%include "faunus/drift.h"
%include "faunus/analysis.h"
%include "faunus/physconst.h"
%include "faunus/titrate.h"
%include "faunus/mcloop.h"
%include "faunus/histogram.h"
%include "faunus/io.h"
//%include "faunus/mpi.h"

/* -----------------------------------
 *              Templates
 * ----------------------------------- */
%template(Average_int) Faunus::Average<int>;
%template(Average_dbl) Faunus::Average<double>;
%template(Average_flt) Faunus::Average<float>;
%template(vecAverage_int) std::vector< Faunus::Average<int> >;
%template(vecAverage_dbl) std::vector< Faunus::Average<double> >;
%template(vecAverage_flt) std::vector< Faunus::Average<float> >;
%template(vecParticle) std::vector< Faunus::particle >;
%template(vecGroup) std::vector<Faunus::Group>;
%template(vecGroupMoleculer) std::vector<Faunus::GroupMolecular>;
%template(vecGroupAtomic) std::vector<Faunus::GroupAtomic>;

/* Energy classes */
%template(Nonbonded_DebyeHuckelLJ_Cuboid)
Faunus::Energy::Nonbonded<Faunus::Potential::DebyeHuckelLJ, Faunus::Geometry::Cuboid>;
%template(Nonbonded_DebyeHuckelLJ_Sphere)
Faunus::Energy::Nonbonded<Faunus::Potential::DebyeHuckelLJ, Faunus::Geometry::Sphere>;

%template(Nonbonded_CoulombLJ_Cuboid)
Faunus::Energy::Nonbonded<Faunus::Potential::CoulombLJ, Faunus::Geometry::Cuboid>;
%template(Nonbonded_CoulombLJ_Sphere)
Faunus::Energy::Nonbonded<Faunus::Potential::CoulombLJ, Faunus::Geometry::Sphere>;

%template(Nonbonded_CoulombHS_Sphere)
Faunus::Energy::Nonbonded<Faunus::Potential::CoulombHS, Faunus::Geometry::Sphere>;

/* Simulate typedefs */
//%pythoncode { particle = PointParticle }
//%pythoncode { p_vec = vecParticle }


