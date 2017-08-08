#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <src/core.h>

namespace py = pybind11;
using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;
typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
typedef Space<Geometry::Cuboid, Tparticle> Tspace;
typedef AtomData<Tparticle> Tatomdata;
 
PYBIND11_PLUGIN(pyfaunus)
{
    py::module m("pyfaunus", "pybind11 faunus plugin");

    // Geometries
    py::class_<Geometry::GeometryBase>(m, "Geometrybase")
        .def("getVolume", &Geometry::GeometryBase::getVolume)
        .def("setVolume", &Geometry::GeometryBase::setVolume)
        .def("vdist", &Geometry::GeometryBase::vdist);

    py::class_<Geometry::Cuboid, Geometry::GeometryBase>(m, "Cuboid")
        .def(py::init<>());

    // Particle properties
    py::class_<Radius>(m, "Radius")
        .def(py::init<>())
        .def_readwrite("radius", &Radius::radius);

    py::class_<Charge>(m, "Charge")
        .def(py::init<>())
        .def_readwrite("charge", &Charge::charge);

    py::class_<Tparticle, Radius, Charge>(m, "Particle")
        .def(py::init<>())
        .def_readwrite("pos", &Tparticle::pos);

    // Particle vector
    py::bind_vector<std::vector<Tparticle>>(m, "ParticleVector");

    // AtomData
    py::class_<Tatomdata>(m, "AtomData")
        .def(py::init<>())
        .def_readwrite("name", &Tatomdata::name)
        .def_readwrite("activity", &Tatomdata::activity)
        .def_readwrite("p", &Tatomdata::p);

    // AtomData vector and static instance
    py::bind_vector<std::vector<Tatomdata>>(m, "AtomDataVector");

    m.attr("atoms") = &atoms<Tparticle>;

    // Space
    py::class_<Tspace>(m, "Space")
        .def_readwrite("geo", &Tspace::geo)
        .def("findMolecules", &Tspace::findMolecules); 

    return m.ptr();
}
