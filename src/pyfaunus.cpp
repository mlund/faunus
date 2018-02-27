#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <src/core.h>
#include <src/space.h>
#include <src/move.h>

namespace py = pybind11;
using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;
typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;

typedef Space<Geometry::Cuboid, Tparticle> Tspace;
typedef typename Tspace::Tpvec Tpvec;
typedef typename Tspace::Tgroup Tgroup;
typedef AtomData<Tparticle> Tatomdata;
 
PYBIND11_MODULE(pyfaunus, m)
{
    using namespace pybind11::literals;

    // Geometries
    py::class_<Geometry::GeometryBase>(m, "Geometrybase")
        .def("getVolume", &Geometry::GeometryBase::getVolume, "Get container volume", "dim"_a=3)
        .def("setVolume", &Geometry::GeometryBase::setVolume, "Set container volume", "volume"_a, "directions"_a)
        .def("randompos", &Geometry::GeometryBase::randompos, "Generates a random position within container", "point"_a, "randomfunc"_a)
        .def("vdist", &Geometry::GeometryBase::vdist, "Minimum vector distance, a-b", "a"_a, "b"_a);

    py::class_<Geometry::Box, Geometry::GeometryBase>(m, "Box")
        .def("getLength", &Geometry::Box::getLength, "Get cuboid sidelengths")
        .def("setLength", &Geometry::Box::setLength, "Set cuboid sidelengths", "sides"_a);

    py::class_<Geometry::Cuboid, Geometry::Box>(m, "Cuboid")
        .def(py::init<>());

    // Particle properties
    py::class_<Radius>(m, "Radius")
        .def(py::init<>())
        .def_readwrite("radius", &Radius::radius);

    py::class_<Charge>(m, "Charge")
        .def(py::init<>())
        .def_readwrite("charge", &Charge::charge, "Particle charge (monopole)");

    py::class_<Tparticle, Radius, Charge>(m, "Particle")
        .def(py::init<>())
        .def_readwrite("pos", &Tparticle::pos, "Particle position");

    // Particle vector and it's iterator
    py::class_<typename Tpvec::iterator>(m, "ParticleVectorIterator")
        .def("__add__", [](typename Tpvec::iterator it, int i){ return it+i; } )
        .def("__sub__", [](typename Tpvec::iterator it, int i){ return it-i; } );

    auto _pvec = py::bind_vector<Tpvec>(m, "ParticleVector");
    _pvec
        .def("positions", [](Tpvec &p){ return asEigenMatrix(p.begin(), p.end(), &Tparticle::pos); })
        .def("charges", [](Tpvec &p){ return asEigenVector(p.begin(), p.end(), &Tparticle::charge); })
        .def("begin", [](Tpvec &p){ return p.begin(); })
        .def("end", [](Tpvec &p){ return p.end(); });

    // Group
    py::class_<Tgroup>(m, "Group")
        .def(py::init<typename Tpvec::iterator, typename Tpvec::iterator>())
        .def(py::init<Tpvec&>())
        .def_readwrite("id", &Tgroup::id, "Molecule id")
        .def_readwrite("cm", &Tgroup::cm, "Center of mass")
        .def_readwrite("atomic", &Tgroup::atomic)
        .def("__iter__", [](Tgroup &v) {
                return py::make_iterator(v.begin(), v.end()); },
                py::keep_alive<0, 1>());

    py::bind_vector<std::vector<Tgroup>>(m, "GroupVector");

    // AtomData
    py::class_<Tatomdata>(m, "AtomData")
        .def(py::init<>())
        .def_readwrite("name", &Tatomdata::name)
        .def_readwrite("activity", &Tatomdata::activity, "Activity = chemical potential in log scale (mol/l)")
        .def_readwrite("p", &Tatomdata::p);

    py::bind_vector<std::vector<Tatomdata>>(m, "AtomDataVector");
    m.attr("atoms") = &atoms<Tparticle>;

    // Change
    py::class_<Change>(m, "Change")
        .def(py::init<>())
        .def_readwrite("dV", &Change::dV);

    // Space
    py::class_<Tspace>(m, "Space")
        .def(py::init<>())
        .def_readwrite("geo", &Tspace::geo)
        .def_readwrite("p", &Tspace::p)
        .def_readwrite("groups", &Tspace::groups)
        .def("findMolecules", &Tspace::findMolecules); 
}
