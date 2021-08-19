#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>

#include <core.h>
#include <space.h>
#include <move.h>
#include <analysis.h>
#include <potentials.h>
#include <regions.h>
#include <montecarlo.h>
#include <energy.h>

namespace py = pybind11;
using namespace Faunus;

typedef typename Space::Tpvec Tpvec;
typedef typename Space::Tgroup Tgroup;
typedef Energy::Hamiltonian Thamiltonian;
typedef MetropolisMonteCarlo Tmcsimulation;

inline json dict2json(py::dict dict) {
    py::object dumps = py::module::import("json").attr("dumps");
    return json::parse( dumps(dict).cast<std::string>() );
} // python dict --> c++ json

inline json list2json(py::list list) {
    py::object dumps = py::module::import("json").attr("dumps");
    return json::parse(dumps(list).cast<std::string>());
} // python list --> c++ json

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
        .def("getVolume", &Geometry::GeometryBase::getVolume, "Get container volume", "dim"_a = 3)
        .def("setVolume", &Geometry::GeometryBase::setVolume, "Set container volume", "volume"_a,
             "method"_a = Geometry::ISOTROPIC)
        .def("collision", &Geometry::GeometryBase::collision, "pos"_a, "Checks if point is inside container")
        .def("getLength", &Geometry::GeometryBase::getLength, "Get cuboid sidelengths")
        .def("vdist", &Geometry::GeometryBase::vdist, "Minimum vector distance, a-b", "a"_a, "b"_a)
        .def("randompos",
             [](Geometry::GeometryBase &g, Random &rnd) {
                 Point pos;
                 g.randompos(pos, rnd);
                 return pos;
             })
        .def(
            "boundary",
            [](Geometry::GeometryBase &g, const Point &pos) {
                Point a = pos; // we cannot modify `pos` directly
                g.boundary(a); // as in c++
                return a;
            },
            R"(
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
        .def(py::init([](py::dict dict) {
            auto ptr = std::make_unique<Geometry::Chameleon>();
            Faunus::Geometry::from_json(dict2json(dict), *ptr);
            return ptr;
        }))
        .def("sqdist", &Geometry::Chameleon::sqdist, "Squared minimum distance, |a-b|^2", "a"_a, "b"_a);

    // Particle properties
    py::class_<Charge>(m, "Charge")
        .def(py::init<>())
        .def_readwrite("charge", &Charge::charge, "Particle charge (monopole)");

    py::class_<Particle>(m, "Particle")
        .def(py::init<>())
        .def("traits", &Particle::traits)
        .def_readwrite("id", &Particle::id, "Particle ID")
        .def_readwrite("pos", &Particle::pos, "Particle position")
        .def_readwrite("charge", &Particle::charge, "Particle charge (monopole)");

    // Particle vector and it's iterator
    py::class_<typename Tpvec::iterator>(m, "ParticleVectorIterator")
        .def("__add__", [](typename Tpvec::iterator it, int i){ return it+i; } )
        .def("__sub__", [](typename Tpvec::iterator it, int i){ return it-i; } );

    auto _pvec = py::bind_vector<Tpvec>(m, "ParticleVector");
    _pvec.def("positions", [](Tpvec &p) { return asEigenMatrix(p.begin(), p.end(), &Particle::pos); })
        .def("charges", [](Tpvec &p) { return asEigenVector(p.begin(), p.end(), &Particle::charge); })
        .def("begin", [](Tpvec &p) { return p.begin(); })
        .def("end", [](Tpvec &p) { return p.end(); });

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

    // Region
    py::enum_<Region::RegionBase::RegionType>(m, "RegionType")
        .value("SPHERE", Region::RegionBase::RegionType::SPHERE)
        .value("CUBOID", Region::RegionBase::RegionType::CUBOID)
        .value("WITHIN", Region::RegionBase::RegionType::WITHIN)
        .value("NONE", Region::RegionBase::RegionType::NONE)
        .export_values();

    py::class_<Region::RegionBase>(m, "RegionBase")
        .def_readwrite("name", &Region::RegionBase::name)
        .def_readwrite("type", &Region::RegionBase::type)
        .def("volume", &Region::RegionBase::isInside)
        .def("isInside", &Region::RegionBase::isInside);

    // AtomData
    py::class_<AtomData>(m, "AtomData")
        .def(py::init<>())
        .def_property(
            "eps", [](const AtomData& a) { return a.interaction.at("eps"); },
            [](AtomData& a, double val) { a.interaction.at("eps") = val; })
        .def_property(
            "sigma", [](const AtomData& a) { return a.interaction.at("sigma"); },
            [](AtomData& a, double val) { a.interaction.at("sigma") = val; })
        .def_readwrite("name", &AtomData::name)
        .def_readwrite("activity", &AtomData::activity, "Activity = chemical potential in log scale (mol/l)")
        .def("id", (const int& (AtomData::*)() const) & AtomData::id); // explicit signature due to overload in c++

    auto _atomdatavec = py::bind_vector<std::vector<AtomData>>(m, "AtomDataVector");
    _atomdatavec
        .def("from_list", [](std::vector<AtomData> &a, py::list list) {
                Faunus::from_json(list2json(list), a); } );

    m.attr("atoms") = &Faunus::atoms; // global instance

    // Temperature and other globals etc.
    m.def("getTemperature", []() { return pc::temperature; } );
    m.def("setTemperature", [](double T) { pc::temperature = T; } );

    // --------- Pair Potentials ---------

    // Base
    py::class_<Potential::PairPotentialBase>(m, "PairPotentialBase")
        .def_readwrite("name", &Potential::PairPotentialBase::name)
        .def_readwrite("cite", &Potential::PairPotentialBase::cite)
        .def_readwrite("isotropic", &Potential::PairPotentialBase::isotropic)
        .def_readwrite("selfEnergy", &Potential::PairPotentialBase::selfEnergy)
        .def("force", &Potential::PairPotentialBase::force)
        .def("energy", [](Potential::PairPotentialBase &pot, const Particle &a, const Particle &b, double r2,
                          const Point &r) { return pot(a, b, r2, r); });

    // Potentials::FunctorPotential
    py::class_<Potential::FunctorPotential, Potential::PairPotentialBase>(m, "FunctorPotential")
        .def(py::init([](py::dict dict) { return from_dict<Potential::FunctorPotential>(dict); }));

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
    py::class_<Space>(m, "Space")
        .def(py::init<>())
        .def_readwrite("geo", &Space::geo)
        .def_readwrite("p", &Space::p)
        .def_readwrite("groups", &Space::groups)
        .def("findMolecules", &Space::findMolecules)
        .def("from_dict", [](Space &spc, py::dict dict) { from_json(dict2json(dict), spc); });

    // Hamiltonian
    py::class_<Thamiltonian>(m, "Hamiltonian")
        .def(py::init<Space &, const json &>())
        .def(py::init([](Space &spc, py::list list) {
            json j = list2json(list);
            return std::unique_ptr<Thamiltonian>(new Thamiltonian(spc, j));
        }))
        .def("init", &Thamiltonian::init)
        .def("energy", &Thamiltonian::energy);

    // TranslationalEntropy
    py::class_<TranslationalEntropy>(m, "TranslationalEntropy")
        .def(py::init<Space &, Space &>())
        .def("energy", &TranslationalEntropy::energy);

    // MCSimulation
    py::class_<Tmcsimulation>(m, "MetropolisMonteCarlo")
        .def(py::init([](py::dict dict) {
            json j = dict2json(dict);
            return std::unique_ptr<Tmcsimulation>(new Tmcsimulation(j, Faunus::MPI::mpi));
        }))
        .def(py::init([](py::dict dict, Faunus::MPI::MPIController &mpi) {
            json j = dict2json(dict);
            return std::unique_ptr<Tmcsimulation>(new Tmcsimulation(j, mpi));
        }));

    // Analysisbase
    py::class_<Analysis::Analysisbase>(m, "Analysisbase")
        .def_readwrite("name", &Analysis::Analysisbase::name)
        .def_readwrite("cite", &Analysis::Analysisbase::cite)
        .def("to_disk", &Analysis::Analysisbase::to_disk)
        .def("sample", &Analysis::Analysisbase::sample)
        .def("to_dict", [](Analysis::Analysisbase &self) {
            json j;
            Analysis::to_json(j, self);
            return json2dict(j);
        });

    py::bind_vector<std::vector<std::shared_ptr<Analysis::Analysisbase>>>(m, "AnalysisVector");

    // CombinedAnalysis
    py::class_<Analysis::CombinedAnalysis>(m, "Analysis")
        .def(py::init([](Space &spc, Thamiltonian &pot, py::list list) {
            json j = list2json(list);
            return std::unique_ptr<Analysis::CombinedAnalysis>(new Analysis::CombinedAnalysis(j, spc, pot));
        }))
        .def_readwrite("vector", &Analysis::CombinedAnalysis::vec)
        .def("to_dict",
             [](Analysis::CombinedAnalysis &self) {
                 json j;
                 Faunus::to_json(j, self);
                 return json2dict(j);
             })
        .def("sample", &Analysis::CombinedAnalysis::sample);
}
