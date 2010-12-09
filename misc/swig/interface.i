%module pyfaunus
%{
#include "faunus/faunus.h"
#include "faunus/potentials/pot_coulomb.h"
#include "faunus/potentials/pot_hscoulomb.h"
#include "faunus/potentials/pot_minimage.h"
#include "faunus/potentials/pot_test_rf.h"
#include "faunus/potentials/pot_silkspider.h"
%}

using namespace std;

/* -----------------------------------
 *          Do some renaming
 * ----------------------------------- */
%rename(assign) Faunus::point::operator=;
%rename(assign) Faunus::particle::operator=;
%rename(assign) Faunus::spherical::operator=;
%rename(assign) Faunus::macromolecule::operator=;
%rename(_print) Faunus::grandcanonical::print;
%rename(_print) Faunus::clusterrotate::print;
%rename(_print) Faunus::inputfile::print;
%rename(__getitem__) Faunus::species::operator[];
%ignore Faunus::point::operator << (std::ostream &, Faunus::point &);

/* -----------------------------------
 *            Header import
 * ----------------------------------- */
%include "std_vector.i"
%include "std_string.i"
%include "std_iostream.i"
//%include "std_sstream.i"
//%include "std_basic_string.i"
#include "faunus/physconst.h"
#include "faunus/faunus.h"
#include "faunus/slump.h"
#include "faunus/point.h"
#include "faunus/average.h"
#include "faunus/particles.h"
#include "faunus/group.h"
#include "faunus/inputfile.h"
#include "faunus/species.h"
#include "faunus/container.h"
#include "faunus/xytable.h"
#include "faunus/ensemble.h"
#include "faunus/io.h"
#include "faunus/analysis.h"
#include "faunus/widom.h"
#include "faunus/mcloop.h"
#include "faunus/titrate.h"
#include "faunus/histogram.h"
#include "faunus/notification.h"
#include "faunus/hardsphere.h"
#include "faunus/moves/base.h"
#include "faunus/moves/translational.h"
#include "faunus/moves/rotational.h"
#include "faunus/moves/charge.h"
#include "faunus/moves/volume.h"
#include "faunus/moves/miscmove.h"
//#include "faunus/moves/clustermove.h"
#include "faunus/moves/crankshaft.h"
#include "faunus/moves/saltbath.h"
#include "faunus/moves/eqtitrate.h"
#include "faunus/moves/replicaexchange.h"
#include "faunus/energy/base.h"
#include "faunus/energy/springinteraction.h"
#include "faunus/bottles/base.h"
#include "faunus/bottles/npt_molecular.h"

/* -----------------------------------
 *        Define some templates
 * ----------------------------------- */
%template(average_int) Faunus::average<int>;
%template(average_dbl) Faunus::average<double>;
%template(vector_group) std::vector<Faunus::group>;
%template(vector_polymer) std::vector<Faunus::polymer>;
%template(vector_particle) std::vector<Faunus::particle>;
%template(vector_macromolecule) std::vector<Faunus::macromolecule>;
namespace std {
  %template(vector_int) std::vector<int>;
  %template(vector_dbl) std::vector<double>;
}

/* -----------------------------------
 *           Energy functions
 * ----------------------------------- */
%feature("notabstract") interaction;
%feature("notabstract") springinteraction;

#include "faunus/potentials/pot_coulomb.h"
%template(interaction_coulomb) Faunus::interaction<Faunus::pot_coulomb>;

#include "faunus/potentials/pot_hscoulomb.h"
%template(interaction_hscoulomb) Faunus::interaction<Faunus::pot_hscoulomb>;

#include "faunus/potentials/pot_minimage.h"
%template(interaction_r12minimage) Faunus::interaction<Faunus::pot_r12minimage>;
%template(springinteraction_r12minimage) Faunus::springinteraction<Faunus::pot_r12minimage>;
%template(interaction_minimage) Faunus::interaction<Faunus::pot_minimage>;

#include "faunus/potentials/pot_silkspider.h"
%template(interaction_silkspider) Faunus::interaction<Faunus::pot_silkspider>;
%template(springinteraction_silkspider) Faunus::springinteraction<Faunus::pot_silkspider>;

#include "faunus/potentials/pot_test_rf.h"
//%template(interaction_testrf) Faunus::interaction<Faunus::pot_test_rf>;


