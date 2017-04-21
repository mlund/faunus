#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/json.h>
#include <faunus/inputfile.h>
#include <faunus/mcloop.h>
#include <faunus/space.h>
#include <faunus/move.h>
#include <faunus/energy.h>

namespace py = pybind11;
using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;
typedef Space<Tgeometry> Tspace;
typedef typename Tspace::ParticleVector Tpvec;

PYBIND11_PLUGIN(pyfaunus)
{
    py::module m("pyfaunus", "pybind11 faunus plugin");

    // Point
    py::class_ <Point> _point(m, "Point");

    _point.def(py::init<>());

    _point.def_property("x", []( Point &b ) { return &b.x(); }, []( Point &b, double x ) { b.x() = x; });
    _point.def_property("y", []( Point &b ) { return &b.y(); }, []( Point &b, double y ) { b.y() = y; });
    _point.def_property("z", []( Point &b ) { return &b.z(); }, []( Point &b, double z ) { b.z() = z; });

    _point.def_buffer([]( Point &m ) -> py::buffer_info
                      {
                          return py::buffer_info(
                              m.data(),                              /* Pointer to buffer */
                              sizeof(double),                         /* Size of one scalar */
                              py::format_descriptor<double>::value(), /* Python struct-style format descriptor */
                              1,                                     /* Number of dimensions */
                              {(unsigned long) (m.rows()),
                               (unsigned long) (m.cols())},                /* Buffer dimensions */
                              {sizeof(double) * m.cols(),            /* Strides (in bytes) for each index */
                               sizeof(double)}
                          );
                      });

    // PointParticle
    py::class_ <PointParticle> _pointparticle(m, "PointParticle", _point);
    _pointparticle.def(py::init<>());
    _pointparticle.def_readwrite("charge", &PointParticle::charge);
    _pointparticle.def_readwrite("radius", &PointParticle::radius);
    _pointparticle.def_readwrite("mw", &PointParticle::mw);
    _pointparticle.def("volume", &PointParticle::volume);

    py::implicitly_convertible<PointParticle, Point>();

    // Tmjson and InputMap
    py::class_ <Tmjson> js(m, "Tmjson");

    js.def(py::init<>());

    py::class_<InputMap>(m, "InputMap", js)
        .def(py::init<>())
        .def(py::init<string>());

    // MCLoop
    py::class_<MCLoop>(m, "MCLoop")
        .def("info", &MCLoop::info)
        .def(py::init<Tmjson &>());

    // Group
    py::class_<Group>(m, "Group")
        .def(py::init<int, int>())
        .def("__len__", []( Group &g ) { return g.size(); })
        .def_readwrite("name", &Group::name)
        .def_readwrite("molId", &Group::molId)
        .def_readwrite("cm", &Group::cm)
        .def_readwrite("cm_trial", &Group::cm_trial)
        .def("range", []( Group &g )
        {
            py::list l;
            for ( auto i : g )
                l.append(py::cast(i));
            return l;
        })
        .def("front", []( Group &g ) { return g.front(); })
        .def("back", []( Group &g ) { return g.back(); });

    // Geometries
    py::class_ <Geometry::Geometrybase> _Geometrybase(m, "Geometrybase");
    _Geometrybase.def("info", &Geometry::Geometrybase::info);
    _Geometrybase.def("getVolume", &Geometry::Geometrybase::getVolume);
    _Geometrybase.def("setVolume", &Geometry::Geometrybase::setVolume);
    _Geometrybase.def("dist", &Geometry::Geometrybase::dist);
    _Geometrybase.def("vdist", &Geometry::Geometrybase::vdist);
    _Geometrybase.def("sqdist", &Geometry::Geometrybase::sqdist);
    _Geometrybase.def("boundary", &Geometry::Geometrybase::boundary);
    _Geometrybase.def("randompos", &Geometry::Geometrybase::randompos);
    _Geometrybase.def("collision", &Geometry::Geometrybase::collision);

    py::class_<Geometry::Cuboid>(m, "Cuboid", _Geometrybase).def(py::init<Tmjson &>());

    m.def("massCenter", &Geometry::massCenter<Tgeometry, Tpvec, Group>);
    m.def("cm2origo", &Geometry::cm2origo<Tgeometry, Tpvec>);
    m.def("translate", &Geometry::translate<Tgeometry, Tpvec>);
    m.def("calcVolume", &Geometry::calcVolume<Tpvec>);

    m.def("dipoleMoment",
          &Geometry::dipoleMoment<Tspace, Group>,
          "Calculates the dipole moment of a group");
    //py::arg("cutoff")=1e9, py::arg("mu") = Point());

    // Space
    py::class_<Tspace>(m, "Space")
        .def("info", &Tspace::info)
        .def("save", &Tspace::save)
        .def("load", &Tspace::load)
        .def("atomList", []( Tspace &s ) { return s.atomList(); })
        .def("groupList", []( Tspace &s )
        {
            std::vector<Group> g;
            g.reserve(s.groupList().size());
            for ( auto gPtr : s.groupList())
                g.push_back(*gPtr);
            return g;
        })
        .def_readwrite("p", &Tspace::p, py::return_value_policy::reference_internal)
        .def_readwrite("trial", &Tspace::trial)
        .def_readwrite("geo", &Tspace::geo)
        .def(py::init<Tmjson &>());

    //py::class_<Move::Propagator<Tspace>>(m, "Propagator")
    //    .def("info", &Move::Propagator<Tspace>::info);

    return m.ptr();
}
