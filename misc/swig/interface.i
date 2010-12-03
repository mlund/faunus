%module pyfaunus
%{
#include "faunus/faunus.h"
#include "faunus/potentials/pot_hscoulomb.h"
#include "faunus/potentials/pot_minimage.h"
%}

using namespace std;

/* -----------------------------------
 *           Plain classes
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
#include "faunus/hardsphere.h"
#include "faunus/moves/base.h"
#include "faunus/moves/translational.h"
#include "faunus/moves/rotational.h"
#include "faunus/moves/charge.h"
#include "faunus/moves/volume.h"
#include "faunus/moves/miscmove.h"
#include "faunus/moves/clustermove.h"
#include "faunus/moves/crankshaft.h"
#include "faunus/moves/saltbath.h"
#include "faunus/moves/eqtitrate.h"
//#include "faunus/bottles/base.h"
//#include "faunus/bottles/npt_molecular.h"

/* -----------------------------------
 *        Define some templates
 * ----------------------------------- */
%template(average_int) Faunus::average<int>;
%template(average_dbl) Faunus::average<double>;
%template(vector_group) std::vector<Faunus::group>;
%template(vector_particle) std::vector<Faunus::particle>;
namespace std {
  %template(vector_int) std::vector<int>;
  %template(vector_dbl) std::vector<double>;
}

/* -----------------------------------
 *           Energy functions
 * ----------------------------------- */
%feature("notabstract") interaction;
%feature("notabstract") springinteraction;
#include "faunus/energy/base.h"
#include "faunus/potentials/pot_hscoulomb.h"
%template(interaction_hscoulomb) Faunus::interaction<Faunus::pot_hscoulomb>;
#include "faunus/potentials/pot_minimage.h"
%template(springinteraction_r12minimage) Faunus::springinteraction<Faunus::pot_r12minimage>;

