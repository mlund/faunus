#include <boost/python.hpp>
#include "faunus/faunus.h"
#include "faunus/potentials/pot_coulomb.h"

namespace Faunus {
    typedef interaction<pot_coulomb> interaction_coulomb;
}

BOOST_PYTHON_MODULE( pyfaunus ) {

    using namespace boost::python;

    class_<Faunus::point>("point")
        .def_readwrite("x", &Faunus::point::x)
        .def_readwrite("y", &Faunus::point::y)
        .def_readwrite("z", &Faunus::point::z)
        .def("clear", &Faunus::point::clear)
    ;

    class_<Faunus::inputfile>("inputfile", init<std::string>())
    ;
 
    class_<Faunus::interaction_coulomb>("interaction_coulomb", init<Faunus::inputfile&>())
        .def("dipdip", &Faunus::interaction_coulomb::dipdip)
    ;
}

