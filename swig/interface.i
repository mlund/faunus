%module pyfaunus
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
#include "faunus/mpi.h"
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
%rename(myinsert) Faunus::Space::insert;
%rename(myerase) Faunus::Space::erase;

%ignore *::operator<<(std::istream&);
%ignore *::operator<<(std::ostream&, const Faunus::Point&);
%ignore Faunus::Average::operator=(T);
%ignore Faunus::BlockCorrelation::operator[] (unsigned int);
%ignore Faunus::Move::AtomTracker::operator[] (particle::Tid);

%include "std_vector.i"
%include "std_string.i"
%include "std_iostream.i"
%include "std_streambuf.i"
%include "std_common.i"
%include "std_map.i"
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
%include "faunus/mpi.h"

%template(average_int) Faunus::Average<int>;
%template(average_dbl) Faunus::Average<double>;
%template(average_flt) Faunus::Average<float>;
%template(vector_int) std::vector< Faunus::Average<int> >;
%template(vector_dblavg) std::vector< Faunus::Average<double> >;
%template(vector_fltavg) std::vector< Faunus::Average<float> >;
%template(vector_particle) std::vector< Faunus::particle >;
%template(vector_group) std::vector<Faunus::Group>;
%template(vector_groupmoleculer) std::vector<Faunus::GroupMolecular>;

