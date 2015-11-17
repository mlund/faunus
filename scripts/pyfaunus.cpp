#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

typedef Space<Geometry::Cuboid> Tspace;

PYBIND11_PLUGIN(pyfaunus) {
  py::module m("pyfaunus", "pybind11 faunus plugin");

  py::class_<Eigen::MatrixXd>(m, "MatrixXd")
    .def("__init__", [](Eigen::MatrixXd &m, py::buffer b) {
        /* Request a buffer descriptor from Python */
        py::buffer_info info = b.request();

        /* Some sanity checks ... */
        if (info.format != py::format_descriptor<double>::value())
        throw std::runtime_error("Incompatible format: expected a double array!");

        if (info.ndim != 2)
        throw std::runtime_error("Incompatible buffer dimension!");

        if (info.strides[0] == sizeof(double)) {
        /* Buffer has the right layout -- directly copy. */
        new (&m) Eigen::MatrixXd(info.shape[0], info.shape[1]);
        memcpy(m.data(), info.ptr, sizeof(double) * m.size());
        } else {
        /* Oops -- the buffer is transposed */
        new (&m) Eigen::MatrixXd(info.shape[1], info.shape[0]);
        memcpy(m.data(), info.ptr, sizeof(double) * m.size());
        m.transposeInPlace();
        }
    });

  py::class_<PointParticle>(m, "PointParticle")
    .def(py::init<>())
    .def_readwrite("charge", &PointParticle::charge)
    .def("x", [](PointParticle &b) { return b.x(); })
    .def("volume", &PointParticle::volume);

  py::class_<Tmjson> js(m, "Tmjson");
  js.def(py::init<>());

  py::class_<InputMap>(m, "InputMap", js)
    .def(py::init<>())
    .def(py::init<string>());

  py::class_<MCLoop>(m, "MCLoop")
    .def("info", &MCLoop::info)
    .def(py::init<Tmjson&>());

  py::class_<Tspace>(m, "Space")
    .def("info", &Tspace::info)
    .def("p", [](Tspace &s) { return s.p; } )
    .def("trial", [](Tspace &s) { return s.trial; } )
    .def(py::init<Tmjson&>());

  //py::class_<Move::Propagator<Tspace>>(m, "Propagator")
  //    .def("info", &Move::Propagator<Tspace>::info);

  return m.ptr();
}
