#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>

#include <src/core.h>
#include <src/space.h>
#include <src/move.h>
#include <src/analysis.h>

namespace py = pybind11;
using namespace Faunus;

typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;

typedef Space<Tparticle> Tspace;
typedef typename Tspace::Tpvec Tpvec;
typedef typename Tspace::Tgroup Tgroup;
typedef Energy::Hamiltonian<Tspace> Thamiltonian;
typedef MCSimulation<Tparticle> Tmcsimulation;

inline json dict2json(py::dict dict) {
    py::object dumps = py::module::import("json").attr("dumps");
    return json::parse( dumps(dict).cast<std::string>() );
} // python dict --> c++ json

inline py::dict json2dict(const json &j) {
    py::object loads = py::module::import("json").attr("loads");
    return loads( j.dump() ) ;
}

template<class T>
std::unique_ptr<T> from_dict(py::dict dict) {
    auto ptr = new T();
    *ptr = dict2json(dict);
    return std::unique_ptr<T>(ptr);
} // convert py::dict to T through Faunus::json

PYBIND11_MODULE(pyfaunus, m)
{
    using namespace pybind11::literals;

    // json
    py::class_<json>(m, "json")
        .def(py::init( [](std::string arg) {
                    return std::unique_ptr<json>(new json(json::parse(arg)));
                    } ) )
    .def(py::init([](py::dict dict) {
                py::object dumps = py::module::import("json").attr("dumps");
                std::string s = dumps(dict).cast<std::string>();
                return std::unique_ptr<json>(new json(json::parse(s)));
                } ) );

    // Random
    py::class_<Random>(m, "Random")
        .def(py::init<>())
        .def("__call__", [](Random &self){ return self(); } ); // function operator

    m.attr("random") = &Faunus::random; // global instance

    // Geometries
    py::enum_<Geometry::VolumeMethod>(m, "VolumeMethod")
        .value("ISOTROPIC", Geometry::VolumeMethod::ISOTROPIC)
        .value("ISOCHORIC", Geometry::VolumeMethod::ISOCHORIC)
        .value("XY", Geometry::VolumeMethod::XY)
        .value("Z", Geometry::VolumeMethod::Z)
        .export_values();

    py::class_<Geometry::GeometryBase>(m, "Geometrybase")
        .def("getVolume", &Geometry::GeometryBase::getVolume, "Get container volume", "dim"_a=3)
        .def("setVolume", &Geometry::GeometryBase::setVolume, "Set container volume", "volume"_a,
                "method"_a=Geometry::ISOTROPIC)
        .def("collision", &Geometry::GeometryBase::collision, "pos"_a, "Checks if point is inside container")
        .def("getLength", &Geometry::GeometryBase::getLength, "Get cuboid sidelengths")
        .def("vdist", &Geometry::GeometryBase::vdist, "Minimum vector distance, a-b", "a"_a, "b"_a)
        .def("sqdist", &Geometry::GeometryBase::sqdist, "Squared minimum distance, |a-b|^2", "a"_a, "b"_a)
        .def("randompos", [](Geometry::GeometryBase &g, Random &rnd) { 
                Point pos;
                g.randompos(pos, rnd);
                return pos;
                })
        .def("boundary", [](Geometry::GeometryBase &g, const Point &pos) {
                Point a = pos; // we cannot modify `pos` directly
                g.boundary(a); // as in c++
                return a;
                }, R"(
                    Apply periodic boundaries

                    If applicable for the geometry type, this will apply
                    periodic boundaries to a coordinate.

                    Args:
                       pos (array): coordinate

                    Returns:
                       array: position with wrapped coordinate
                )");

    py::class_<Geometry::Chameleon, Geometry::GeometryBase>(m, "Chameleon")
        .def(py::init<>())
        .def(py::init( [](py::dict dict) {
                    auto ptr = std::make_unique<Geometry::Chameleon>();
                    Faunus::Geometry::from_json( dict2json(dict), *ptr);
                    return ptr;
                    } ) );

    // Particle properties
    py::class_<Radius>(m, "Radius")
        .def(py::init<>())
        .def_readwrite("radius", &Radius::radius);

    py::class_<Charge>(m, "Charge")
        .def(py::init<>())
        .def_readwrite("charge", &Charge::charge, "Particle charge (monopole)");

    py::class_<Tparticle, Radius, Charge>(m, "Particle")
        .def(py::init<>())
        .def("traits", &Tparticle::traits)
        .def_readwrite("id", &Tparticle::id, "Particle ID")
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
        .def_readwrite("groups", &Tgroup::id, "Molecule id")
        .def_readwrite("id", &Tgroup::id, "Molecule id")
        .def_readwrite("cm", &Tgroup::cm, "Center of mass")
        .def_readwrite("atomic", &Tgroup::atomic)
        .def("__len__", [](Tgroup &self) { return self.size(); } )
        .def("__iter__", [](Tgroup &v) {
                return py::make_iterator(v.begin(), v.end()); },
                py::keep_alive<0, 1>())
        .def("traits", &Tgroup::traits)
        .def("contains", &Tgroup::contains)
        .def("capacity", &Tgroup::capacity)
        .def("deactivate", &Tgroup::deactivate)
        .def("activate", &Tgroup::activate)
        .def("begin", (Tpvec::iterator& (Tgroup::*)() ) &Tgroup::begin)
        .def("end", (Tpvec::iterator& (Tgroup::*)() ) &Tgroup::end);

    py::bind_vector<std::vector<Tgroup>>(m, "GroupVector");

    // AtomData
    py::class_<AtomData>(m, "AtomData")
        .def(py::init<>())
        .def_readwrite("eps", &AtomData::eps)
        .def_readwrite("sigma", &AtomData::sigma)
        .def_readwrite("name", &AtomData::name)
        .def_readwrite("activity", &AtomData::activity, "Activity = chemical potential in log scale (mol/l)")
        .def("id", (const int& (AtomData::*)() const) &AtomData::id); // explicit signature due to overload in c++

    auto _atomdatavec = py::bind_vector<std::vector<AtomData>>(m, "AtomDataVector");
    _atomdatavec
        .def("from_dict", [](std::vector<AtomData> &a, py::dict dict) {
                Faunus::from_json(dict2json(dict), a); } );

    m.attr("atoms") = &Faunus::atoms; // global instance

    // Temperature and other globals etc.
    m.def("getTemperature", []() { return pc::temperature; } );
    m.def("setTemperature", [](double T) { pc::temperature = T; } );

    // Potentials
    py::class_<Potential::FunctorPotential<Tparticle>>(m, "FunctorPotential")
        .def(py::init( [](py::dict dict) {
                    return from_dict<Potential::FunctorPotential<Tparticle>>(dict);
                    } ) )
    .def("energy", [](Potential::FunctorPotential<Tparticle> &pot, const Tparticle &a, const Tparticle &b, const Point &r){
            return pot(a,b,r); 
            });

    // Change::Data
    py::class_<Change::data>(m, "ChangeData")
        .def(py::init<>())
        .def_readwrite("dNatomic", &Change::data::dNatomic)
        .def_readwrite("dNswap", &Change::data::dNswap)
        .def_readwrite("index", &Change::data::index)
        .def_readwrite("internal", &Change::data::internal)
        .def_readwrite("all", &Change::data::all)
        .def_readwrite("atoms", &Change::data::atoms);

    py::bind_vector<std::vector<Change::data>>(m, "ChangeDataVec");

    // Change
    py::class_<Change>(m, "Change")
        .def(py::init<>())
        .def_readwrite("all", &Change::all)
        .def_readwrite("dV", &Change::dV)
        .def_readwrite("dN", &Change::dN)
        .def_readwrite("groups", &Change::groups);

    // Space
    py::class_<Tspace>(m, "Space")
        .def(py::init<>())
        .def_readwrite("geo", &Tspace::geo)
        .def_readwrite("p", &Tspace::p)
        .def_readwrite("groups", &Tspace::groups)
        .def("findMolecules", &Tspace::findMolecules)
        .def("from_dict", [](Tspace &spc, py::dict dict) { from_json(dict2json(dict), spc); } );

    // Hamiltonian
    py::class_<Thamiltonian>(m, "Hamiltonian")
        .def(py::init<Tspace&, const json&>())
        .def(py::init([](Tspace &spc, py::dict dict) {
                    json j = dict2json(dict);
                    return std::unique_ptr<Thamiltonian>(new Thamiltonian(spc,j));
                    }))
        .def("init", &Thamiltonian::init)
        .def("energy", &Thamiltonian::energy);

    // IdealTerm
    m.def("IdealTerm", &IdealTerm<Tspace>);

    // MCSimulation
    py::class_<Tmcsimulation>(m, "MCSimulation")
        .def(py::init([](py::dict dict) {
                    json j = dict2json(dict);
                    return std::unique_ptr<Tmcsimulation>(new Tmcsimulation(j,Faunus::MPI::mpi));
                    }))
        .def(py::init([](py::dict dict, Faunus::MPI::MPIController &mpi) {
                    json j = dict2json(dict);
                    return std::unique_ptr<Tmcsimulation>(new Tmcsimulation(j,mpi));
                    }));

    // CombinedAnalysis
    py::class_<Analysis::CombinedAnalysis>(m, "Analysis")
        .def(py::init([](Tspace &spc, Thamiltonian &pot, py::dict dict) {
                    json j = dict2json(dict);
                    return std::unique_ptr<Analysis::CombinedAnalysis>(new Analysis::CombinedAnalysis(j,spc,pot));
                    }))
        .def("to_dict", [](Analysis::CombinedAnalysis &self) {
                json j;
                Faunus::to_json(j, self);
                return json2dict(j);
                } )
        .def("sample", &Analysis::CombinedAnalysis::sample);
}
