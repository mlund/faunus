#include <faunus/faunus.h>
#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/particles.h>
#include <faunus/slump.h>
#include <faunus/inputfile.h>
#include <faunus/species.h>
#include <faunus/io.h>
#include <faunus/container.h>
#include <faunus/potentials/pot_coulomb.h>
#include <faunus/moves/base.h>
#include <faunus/moves/translational.h>
#include <faunus/moves/charge.h>
#include <faunus/group.h>
#include <faunus/analysis.h>
#include <faunus/widom.h>
#include <faunus/energy.h>
#include <faunus/ensemble.h>
#include <faunus/histogram.h>
#include <faunus/xytable.h>
#include <faunus/mcloop.h>
#include <faunus/average.h>

namespace pyplusplus{ namespace aliases{

    typedef Faunus::interaction< Faunus::pot_coulomb > interaction_coulomb;
    typedef Faunus::interaction< Faunus::pot_hscoulomb > interaction_hscoulomb;
    typedef Faunus::average< float > average_float;
    typedef Faunus::average< double > average_double;

    inline void instantiate() {
        sizeof(Faunus::interaction< Faunus::pot_coulomb >);
        sizeof(Faunus::interaction< Faunus::pot_hscoulomb >);
        sizeof(Faunus::average< float >);
        sizeof(Faunus::average< double >);
    }
} }
