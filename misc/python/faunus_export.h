#ifndef __faunus_export_h__
#define __faunus_export_h__

#include "faunus/point.h"
#include "faunus/particles.h"
#include "faunus/slump.h"
#include "faunus/inputfile.h"
#include "faunus/species.h"
#include "faunus/io.h"
#include "faunus/container.h"
#include "faunus/potentials/pot_coulomb.h"
#include "faunus/moves/base.h"
#include "faunus/moves/translational.h"
#include "faunus/group.h"
#include "faunus/analysis.h"
#include "faunus/widom.h"
#include "faunus/energy.h"
#include "faunus/ensemble.h"


namespace pyplusplus{ namespace aliases{

    typedef Faunus::interaction<Faunus::pot_coulomb> interaction_coulomb;

    inline void instantiate() {
        sizeof(Faunus::interaction<Faunus::pot_coulomb>);
    }

} }


#endif // __faunus_export_hpp__
