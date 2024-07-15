#include <memory>
#include <pybind11_json/pybind11_json.hpp>
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

using Tgroup = typename Space::GroupType;
using Thamiltonian = Energy::Hamiltonian;
using Tmcsimulation = MetropolisMonteCarlo;

template <class T> std::unique_ptr<T> from_dict(py::dict dict)
{
    auto ptr = std::make_unique<T>();
    *ptr = static_cast<T>(json(dict));
    return ptr;
} // convert py::dict to T through Faunus::json

PYBIND11_MODULE(pyfaunus, m)
{
    using namespace pybind11::literals;

    // Random
    py::class_<Random>(m, "Random").def(py::init<>()).def("__call__", [](Random& self) {
        return self();
    }); // function operator

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
             "method"_a = Geometry::VolumeMethod::ISOTROPIC)
        .def("collision", &Geometry::GeometryBase::collision, "pos"_a,
             "Checks if point is inside container")
        .def("getLength", &Geometry::GeometryBase::getLength, "Get cuboid sidelengths")
        .def("vdist", &Geometry::GeometryBase::vdist, "Minimum vector distance, a-b", "a"_a, "b"_a)
        .def("randompos",
             [](Geometry::GeometryBase& g, Random& rnd) {
                 Point pos;
                 g.randompos(pos, rnd);
                 return pos;
             })
        .def(
            "boundary",
            [](Geometry::GeometryBase& g, const Point& pos) {
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
            Faunus::Geometry::from_json(dict, *ptr);
            return ptr;
        }))
        .def("sqdist", &Geometry::Chameleon::sqdist, "Squared minimum distance, |a-b|^2", "a"_a,
             "b"_a);

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
    py::class_<ParticleVector::iterator>(m, "ParticleVectorIterator")
        .def("__add__", [](ParticleVector::iterator it, int i) { return it + i; })
        .def("__sub__", [](ParticleVector::iterator it, int i) { return it - i; });

    auto _pvec = py::bind_vector<ParticleVector>(m, "ParticleVector");
    _pvec
        .def("positions",
             [](ParticleVector& particles) {
                 return asEigenMatrix(particles.begin(), particles.end(), &Particle::pos);
             })
        .def("charges",
             [](ParticleVector& particles) {
                 return asEigenVector(particles.begin(), particles.end(), &Particle::charge);
             })
        .def("begin", [](ParticleVector& particles) { return particles.begin(); })
        .def("end", [](ParticleVector& particles) { return particles.end(); });

    // Group
    py::class_<Tgroup>(m, "Group")
        .def(py::init<MoleculeData::index_type, ParticleVector::iterator,
                      ParticleVector::iterator>())
        .def_readwrite("groups", &Tgroup::id, "Molecule id")
        .def_readwrite("id", &Tgroup::id, "Molecule id")
        .def_readwrite("cm", &Tgroup::mass_center, "Center of mass")
        .def("__len__", [](Tgroup& self) { return self.size(); })
        .def(
            "__iter__", [](Tgroup& v) { return py::make_iterator(v.begin(), v.end()); },
            py::keep_alive<0, 1>())
        .def("isAtomic", &Tgroup::isAtomic)
        .def("isMolecular", &Tgroup::isMolecular)
        .def("traits", &Tgroup::traits)
        .def("contains", &Tgroup::contains)
        .def("capacity", &Tgroup::capacity)
        .def("deactivate", &Tgroup::deactivate)
        .def("activate", &Tgroup::activate)
        .def("begin", (ParticleVector::iterator & (Tgroup::*)()) & Tgroup::begin)
        .def("end", (ParticleVector::iterator & (Tgroup::*)()) & Tgroup::end);

    py::bind_vector<std::vector<Tgroup>>(m, "GroupVector");

    // Region
    py::enum_<Region::RegionType>(m, "RegionType")
        .value("WITHIN_PARTICLE", Region::RegionType::WITHIN_PARTICLE)
        .value("WITHIN_MOLID", Region::RegionType::WITHIN_MOLID)
        .value("WITHIN_ELLIPSOID", Region::RegionType::WITHIN_ELLIPSOID)
        .value("NONE", Region::RegionType::INVALID)
        .export_values();

    py::class_<Region::RegionBase>(m, "RegionBase")
        .def_readonly("type", &Region::RegionBase::type)
        .def("volume", &Region::RegionBase::volume);

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
        .def_readwrite("activity", &AtomData::activity,
                       "Activity = chemical potential in log scale (mol/l)")
        .def("id", (const AtomData::index_type& (AtomData::*)() const) &
                       AtomData::id); // explicit signature due to overload in c++

    auto _atomdatavec = py::bind_vector<std::vector<AtomData>>(m, "AtomDataVector");
    _atomdatavec.def("from_list",
                     [](std::vector<AtomData>& a, py::list list) { Faunus::from_json(list, a); });

    m.attr("atoms") = &Faunus::atoms; // global instance

    // Temperature and other globals etc.
    m.def("getTemperature", []() { return pc::temperature; });
    m.def("setTemperature", [](double T) { pc::temperature = T; });

    // --------- Pair Potentials ---------

    // Base
    py::class_<pairpotential::PairPotential>(m, "PairPotentialBase")
        .def_readwrite("name", &pairpotential::PairPotential::name)
        .def_readwrite("cite", &pairpotential::PairPotential::cite)
        .def_readwrite("isotropic", &pairpotential::PairPotential::isotropic)
        .def_readwrite("selfEnergy", &pairpotential::PairPotential::selfEnergy)
        .def("force", &pairpotential::PairPotential::force)
        .def("energy", [](pairpotential::PairPotential& pot, const Particle& a, const Particle& b,
                          double r2, const Point& r) { return pot(a, b, r2, r); });

    // Potentials::FunctorPotential
    py::class_<pairpotential::FunctorPotential, pairpotential::PairPotential>(m, "FunctorPotential")
        .def(py::init([](py::dict dict) {
            auto pairpot = pairpotential::makePairPotential<pairpotential::FunctorPotential>(dict);
            return std::make_unique<pairpotential::FunctorPotential>(pairpot);
        }));

    // Change::Data
    py::class_<Change::GroupChange>(m, "ChangeData")
        .def(py::init<>())
        .def_readwrite("dNatomic", &Change::GroupChange::dNatomic)
        .def_readwrite("dNswap", &Change::GroupChange::dNswap)
        .def_readwrite("index", &Change::GroupChange::group_index)
        .def_readwrite("internal", &Change::GroupChange::internal)
        .def_readwrite("all", &Change::GroupChange::all)
        .def_readwrite("atoms", &Change::GroupChange::relative_atom_indices);

    py::bind_vector<std::vector<Change::GroupChange>>(m, "ChangeDataVec");

    // Change
    py::class_<Change>(m, "Change")
        .def(py::init<>())
        .def_readwrite("everything", &Change::everything)
        .def_readwrite("volume_change", &Change::volume_change)
        .def_readwrite("matter_change", &Change::matter_change)
        .def_readwrite("groups", &Change::groups);

    // Space
    py::class_<Space>(m, "Space")
        .def(py::init<>())
        .def_readwrite("geo", &Space::geometry)
        .def_readwrite("particles", &Space::particles)
        .def_readwrite("groups", &Space::groups)
        // https://stackoverflow.com/questions/65812046/disambiguate-non-const-and-const-access-methods-pybind11
        // .def("findMolecules", &Space::findMolecules)
        .def("from_dict", [](Space& spc, py::dict dict) { from_json(dict, spc); });

    // Hamiltonian
    py::class_<Thamiltonian>(m, "Hamiltonian")
        .def(py::init<Space&, const json&>())
        .def(py::init(
            [](Space& spc, py::list list) { return std::make_unique<Thamiltonian>(spc, list); }))
        .def("init", &Thamiltonian::init)
        .def("energy", &Thamiltonian::energy);

    // TranslationalEntropy
    py::class_<TranslationalEntropy>(m, "TranslationalEntropy")
        .def(py::init<Space&, Space&>())
        .def("energy", &TranslationalEntropy::energy);

    // MCSimulation
    py::class_<Tmcsimulation>(m, "MetropolisMonteCarlo").def(py::init([](py::dict dict) {
        return std::make_unique<Tmcsimulation>(dict);
    }));

    // Analysisbase
    py::class_<analysis::Analysis>(m, "Analysisbase")
        .def_readonly("name", &analysis::Analysis::name)
        .def_readwrite("cite", &analysis::Analysis::cite)
        .def("to_disk", &analysis::Analysis::to_disk)
        .def("sample", &analysis::Analysis::sample)
        .def("to_dict", [](analysis::Analysis& self) {
            json j;
            analysis::to_json(j, self);
            return py::dict(j);
        });

    py::bind_vector<std::vector<std::shared_ptr<analysis::Analysis>>>(m, "AnalysisVector");

    // CombinedAnalysis
    py::class_<analysis::CombinedAnalysis>(m, "Analysis")
        .def(py::init([](Space& spc, Thamiltonian& pot, py::list list) {
            return std::make_unique<analysis::CombinedAnalysis>(list, spc, pot);
        }))
        .def_readwrite("vector", &analysis::CombinedAnalysis::vec)
        .def("to_dict",
             [](analysis::CombinedAnalysis& self) {
                 json j;
                 Faunus::to_json(j, self);
                 return py::dict(j);
             })
        .def("sample", &analysis::CombinedAnalysis::sample);
}
