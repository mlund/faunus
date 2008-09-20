#include <boost/python.hpp>
#include "faunus/point.h"
#include "faunus/energy.h"
#include "faunus/inputfile.h"

namespace Faunus {
    typedef Faunus::interaction<Faunus::pot_coulomb> interaction_coulomb;
}

BOOST_PYTHON_MODULE( pyfaunus ) {

    using namespace boost::python;

/*    def("greet", greet);

    def("fn", fn);

    def("mod", mod);

    class_<VerySimple>("VerySimple")
        .def_readwrite("a", &VerySimple::a)
    ;

    class_<Simple>("Simple")
        .def_readwrite("b", &Simple::b)
        .def_readwrite("v", &Simple::v)
    ;

    class_< Test<Simple> >("TestSimple")
        .def_readwrite("data", &Test<Simple>::data)
    ;
*/

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

