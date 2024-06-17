#include "energy.h"
#include "penalty.h"
#include "potentials.h"
#include "externalpotential.h"
#include <functional>
#include <range/v3/view/zip.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <numeric>
#include <utility>

#ifdef ENABLE_FREESASA
#include <freesasa.h>
struct freesasa_parameters_fwd : public freesasa_parameters {
    freesasa_parameters_fwd(const freesasa_parameters& p) : freesasa_parameters(p){};
};
#endif

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

namespace Faunus::Energy {

EwaldData::EwaldData(const json &j) {
    alpha = j.at("alpha");                          // damping-parameter
    realspace_cutoff = j.at("cutoff");              // real space cut-off
    invspace_cutoff = j.at("ncutoff");              // reciprocal space cut-off
    use_spherical_sum = j.value("spherical_sum", true); // Using spherical summation of k-vectors in reciprocal space?
    bjerrum_length = pc::bjerrumLength(j.at("epsr").get<double>());
    surface_dielectric_constant = j.value("epss", 0.0); // dielectric constant of surrounding medium
    kappa = j.value("kappa", 0.0);
    kappa_squared = kappa * kappa;

    if (j.count("kcutoff")) {
        faunus_logger->warn("`kcutoff` is deprecated, use `ncutoff` instead");
        invspace_cutoff = j.at("kcutoff").get<double>();
    } else {
        invspace_cutoff = j.at("ncutoff").get<double>();
    }
    if (j.value("ipbc", false)) { // look for legacy bool `ipbc`
        faunus_logger->warn("key `ipbc` is deprecated, use `ewaldscheme: ipbc` instead");
        policy = EwaldData::IPBC;
    } else {
        policy = j.value("ewaldscheme", EwaldData::PBC);
        if (policy == EwaldData::INVALID)
            throw std::runtime_error("invalid `ewaldpolicy`");
    }
}

void EwaldData::sync(const EwaldData& other, const Change& change) {
    // hard-coded sync; should be expanded when dipolar ewald is supported
    if (&other == this) {
        return;
    }
    if (change.everything or change.volume_change) {
        Q_dipole.resize(0); // dipoles are currently unsupported
        *this = other;
    } else {
        Q_ion = other.Q_ion;
    }
}

bool EwaldData::tinfoilSurrounding() const {
    return surface_dielectric_constant < 1.0 || std::isinf(surface_dielectric_constant);
}

double EwaldData::volume() const { return box_length.prod(); }

void to_json(json& j, const EwaldData& ewald_data) {
    j = {{"lB", ewald_data.bjerrum_length},
         {"epss", ewald_data.surface_dielectric_constant},
         {"alpha", ewald_data.alpha},
         {"cutoff", ewald_data.realspace_cutoff},
         {"ncutoff", ewald_data.invspace_cutoff},
         {"wavefunctions", ewald_data.k_vectors.cols()},
         {"spherical_sum", ewald_data.use_spherical_sum},
         {"kappa", ewald_data.kappa},
         {"ewaldscheme", ewald_data.policy}};
}

TEST_CASE("[Faunus] Ewald - EwaldData") {
    using doctest::Approx;

    Space spc;
    EwaldData data(R"({
                "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "ncutoff": 11.0, "spherical_sum": true, "cutoff": 5.0})"_json);

    CHECK_EQ(data.policy, EwaldData::PBC);
    CHECK_EQ(data.tinfoilSurrounding(), false);
    CHECK_EQ(data.alpha, 0.894427190999916);

    // Check number of wave-vectors using PBC
    PolicyIonIon ionion;
    ionion.updateBox(data, Point(10, 10, 10));
    CHECK_EQ(data.k_vectors.cols(), 2975);
    CHECK_EQ(data.Q_ion.size(), data.k_vectors.cols());

    // Check number of wave-vectors using IPBC
    data.policy = EwaldData::IPBC;
    PolicyIonIonIPBC ionionIPBC;
    ionionIPBC.updateBox(data, Point(10, 10, 10));
    CHECK_EQ(data.k_vectors.cols(), 846);
    CHECK_EQ(data.Q_ion.size(), data.k_vectors.cols());
}

//----------------- Ewald Policies -------------------

std::unique_ptr<EwaldPolicyBase> EwaldPolicyBase::makePolicy(EwaldData::Policies policy) {
    switch (policy) {
    case EwaldData::PBC:
        return std::make_unique<PolicyIonIon>();
    case EwaldData::PBCEigen:
        return std::make_unique<PolicyIonIonEigen>();
    case EwaldData::IPBC:
        return std::make_unique<PolicyIonIonIPBC>();
    case EwaldData::IPBCEigen:
        return std::make_unique<PolicyIonIonIPBCEigen>();
    case EwaldData::METALSLIT:
        return std::make_unique<PolicyIonIonMetalSlit>();
    default:
        throw std::runtime_error("invalid Ewald policy");
    }
    return nullptr;
}

PolicyIonIon::PolicyIonIon() { cite = "doi:10.1063/1.481216"; }
PolicyIonIonIPBC::PolicyIonIonIPBC() { cite = "doi:10/css8"; }

/**
 * Resize k-vectors according to current variables and box length
 */
void PolicyIonIon::updateBox(EwaldData& ewald_data, const Point& box) const {
    assert(ewald_data.policy == EwaldData::PBC or ewald_data.policy == EwaldData::PBCEigen);
    ewald_data.box_length = box;
    int n_cutoff_ceil = ceil(ewald_data.invspace_cutoff);
    ewald_data.check_k2_zero = 0.1 * std::pow(2 * pc::pi / ewald_data.box_length.maxCoeff(), 2);
    int k_vector_size = (2 * n_cutoff_ceil + 1) * (2 * n_cutoff_ceil + 1) * (2 * n_cutoff_ceil + 1) - 1;
    if (k_vector_size == 0) {
        ewald_data.k_vectors.resize(3, 1);
        ewald_data.Aks.resize(1);
        ewald_data.k_vectors.col(0) = Point(1, 0, 0); // Just so it is not the zero-vector
        ewald_data.Aks[0] = 0;
        ewald_data.num_kvectors = 1;
        ewald_data.Q_ion.resize(1);
        ewald_data.Q_dipole.resize(1);
    } else {
        double nc2 = ewald_data.invspace_cutoff * ewald_data.invspace_cutoff;
        ewald_data.k_vectors.resize(3, k_vector_size);
        ewald_data.Aks.resize(k_vector_size);
        ewald_data.num_kvectors = 0;
        ewald_data.k_vectors.setZero();
        ewald_data.Aks.setZero();
        int start_value = 1;
        for (int nx = 0; nx <= n_cutoff_ceil; nx++) {
            auto dnx2 = double(nx * nx);
            double factor = (nx > 0) ? 2.0 : 1.0; // optimization of PBC Ewald (and
                                                  // always the case for IPBC Ewald)
            for (int ny = -n_cutoff_ceil * start_value; ny <= n_cutoff_ceil; ny++) {
                auto dny2 = double(ny * ny);
                for (int nz = -n_cutoff_ceil * start_value; nz <= n_cutoff_ceil; nz++) {
                    Point kv = 2 * pc::pi * Point(nx, ny, nz).cwiseQuotient(ewald_data.box_length);
                    const auto k2 = kv.squaredNorm() + ewald_data.kappa_squared; // last term is only for Yukawa-Ewald
                    if (k2 < ewald_data.check_k2_zero) {                         // Check if k2 != 0
                        continue;
                    }
                    if (ewald_data.use_spherical_sum) {
                        const auto dnz2 = double(nz * nz);
                        if ((dnx2 + dny2 + dnz2) / nc2 > 1) {
                            continue;
                        }
                    }
                    ewald_data.k_vectors.col(ewald_data.num_kvectors) = kv;
                    ewald_data.Aks[ewald_data.num_kvectors] =
                        factor * exp(-k2 / (4 * ewald_data.alpha * ewald_data.alpha)) / k2;
                    ewald_data.num_kvectors++;
                }
            }
        }
        ewald_data.Q_ion.resize(ewald_data.num_kvectors);
        ewald_data.Q_dipole.resize(ewald_data.num_kvectors);
        ewald_data.Aks.conservativeResize(ewald_data.num_kvectors);
        ewald_data.k_vectors.conservativeResize(3, ewald_data.num_kvectors);
    }
}

/**
 * @todo Add OpenMP pragma to first loop
 */
void PolicyIonIon::updateComplex(EwaldData& ewald_data, const Space::GroupVector& groups) const {
    namespace rv = ranges::cpp20::views;
    for (int k = 0; k < ewald_data.k_vectors.cols(); k++) {
        const Point& q = ewald_data.k_vectors.col(k);
        auto positions = groups | rv::join | rv::transform(&Particle::pos);
        auto charge = groups | rv::join | rv::transform(&Particle::charge);
        ewald_data.Q_ion[k] = sumWavevector(q, ranges::views::zip(positions, charge));
    }
}

void PolicyIonIonEigen::updateComplex(EwaldData& ewald_data, const Space::GroupVector& groups) const {
    auto [pos, charge] = mapGroupsToEigen(groups);            // throws if inactive particles
    Eigen::MatrixXd kr = pos.matrix() * ewald_data.k_vectors; // ( N x 3 ) * ( 3 x K ) = N x K
    ewald_data.Q_ion.real() =
        (kr.array().cos().colwise() * charge).colwise().sum();  // real part of 'Q^q', see eq. 25 in ref.
    ewald_data.Q_ion.imag() = kr.array().sin().colwise().sum(); // imaginary part of 'Q^q', see eq. 25 in ref.
}

EwaldData::Tcomplex PolicyIonIon::sumWavevectorGroups(const Point& wavevector, const Group& group,
                                                      const std::vector<Change::index_type>& indices,
                                                      [[maybe_unused]] EwaldData& ewald_data) const {
    namespace rv = ranges::cpp20::views;
    auto new_indices = indices | rv::filter([&](auto i) { return i < group.size(); }) | ranges::to_vector;
    return sumWavevector(wavevector, ranges::views::zip(group[new_indices] | rv::transform(&Particle::pos),
                                                        group[new_indices] | rv::transform(&Particle::charge)));
}

void PolicyIonIon::updateComplexOptimized(EwaldData& ewald_data, const Change& change, const Space::GroupVector& groups,
                                          const Space::GroupVector& oldgroups) const {
    assert(groups.size() == oldgroups.size());
    namespace rv = ranges::cpp20::views;
    for (int k = 0; k < ewald_data.k_vectors.cols(); k++) {
        auto& Q = ewald_data.Q_ion[k];
        const Point& wavevector = ewald_data.k_vectors.col(k);

        for (const auto& changed_group : change.groups) {
            const auto& g_new = groups.at(changed_group.group_index);
            const auto& g_old = oldgroups.at(changed_group.group_index);
            const auto max_group_size = std::max(g_new.size(), g_old.size());
            auto indices = (changed_group.all) ? ranges::cpp20::views::iota(0U, max_group_size) |
                                                     ranges::to<std::vector<Change::index_type>>
                                               : changed_group.relative_atom_indices;
            Q += sumWavevectorGroups(wavevector, g_new, indices, ewald_data);
            Q -= sumWavevectorGroups(wavevector, g_old, indices, ewald_data);
        }
    }
}

TEST_CASE("[Faunus] Ewald - IonIonPolicy") {
    using doctest::Approx;
    Space spc;
    spc.particles.resize(2);
    spc.geometry = R"( {"type": "cuboid", "length": 10} )"_json;
    spc.particles.at(0) = R"( {"id": 0, "pos": [0,0,0], "q": 1.0} )"_json;
    spc.particles.at(1) = R"( {"id": 0, "pos": [1,0,0], "q": -1.0} )"_json;
    if (Faunus::molecules.empty()) {
        Faunus::molecules.resize(1);
    }
    Group g(0, spc.particles.begin(), spc.particles.end());
    spc.groups.push_back(g);

    auto data = static_cast<EwaldData>(R"({
                "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "ncutoff": 11.0, "spherical_sum": true, "cutoff": 5.0})"_json);
    Change c;
    c.everything = true;
    data.policy = EwaldData::PBC;

    SUBCASE("PBC") {
        PolicyIonIon ionion;
        ionion.updateBox(data, spc.geometry.getLength());
        ionion.updateComplex(data, spc.groups);
        CHECK_EQ(ionion.selfEnergy(data, c, spc.groups), Approx(-1.0092530088080642 * data.bjerrum_length));
        CHECK_EQ(ionion.surfaceEnergy(data, c, spc.groups), Approx(0.0020943951023931952 * data.bjerrum_length));
        CHECK_EQ(ionion.reciprocalEnergy(data), Approx(0.21303063979675319 * data.bjerrum_length));
    }

    SUBCASE("PBCEigen") {
        PolicyIonIonEigen ionion;
        ionion.updateBox(data, spc.geometry.getLength());
        ionion.updateComplex(data, spc.groups);
        CHECK_EQ(ionion.selfEnergy(data, c, spc.groups), Approx(-1.0092530088080642 * data.bjerrum_length));
        CHECK_EQ(ionion.surfaceEnergy(data, c, spc.groups), Approx(0.0020943951023931952 * data.bjerrum_length));
        CHECK_EQ(ionion.reciprocalEnergy(data), Approx(0.21303063979675319 * data.bjerrum_length));
    }

    SUBCASE("IPBC") {
        PolicyIonIonIPBC ionion;
        data.policy = EwaldData::IPBC;
        ionion.updateBox(data, spc.geometry.getLength());
        ionion.updateComplex(data, spc.groups);
        CHECK_EQ(ionion.selfEnergy(data, c, spc.groups), Approx(-1.0092530088080642 * data.bjerrum_length));
        CHECK_EQ(ionion.surfaceEnergy(data, c, spc.groups), Approx(0.0020943951023931952 * data.bjerrum_length));
        CHECK_EQ(ionion.reciprocalEnergy(data), Approx(0.0865107467 * data.bjerrum_length));
    }

    // IPBCEigen is under construction
    /*SUBCASE("IPBCEigen") {
        PolicyIonIonIPBCEigen ionion();
        data.type = EwaldData::IPBCEigen;
        ionion.updateBox(data, spc.geo.getLength());
        ionion.updateComplex(data, spc.groups);
        CHECK_EQ(ionion.selfEnergy(data, c, spc.groups), Approx(-1.0092530088080642 * data.lB));
        CHECK_EQ(ionion.surfaceEnergy(data, c, spc.groups), Approx(0.0020943951023931952 * data.lB));
        CHECK_EQ(ionion.reciprocalEnergy(data), Approx(0.0865107467 * data.lB));
    }*/
}

TEST_CASE("[Faunus] Ewald - IonIonPolicy Benchmarks") {
    Space spc;
    spc.geometry = R"( {"type": "cuboid", "length": 80} )"_json;
    spc.particles.resize(200);
    for (auto& particle : spc.particles) {
        particle.charge = 1.0;
        particle.pos = (random() - 0.5) * spc.geometry.getLength();
    }
    Group group(0, spc.particles.begin(), spc.particles.end());
    spc.groups.push_back(group);

    EwaldData data(R"({
                "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "ncutoff": 11.0, "spherical_sum": true, "cutoff": 9.0})"_json);
    Change change;
    change.everything = true;
    data.policy = EwaldData::PBC;

    {
        PolicyIonIon pbc;
        PolicyIonIonEigen pbc_eigen;
        pbc.updateBox(data, spc.geometry.getLength());
        pbc_eigen.updateBox(data, spc.geometry.getLength());

        ankerl::nanobench::Bench bench;
        bench.minEpochIterations(20);
        bench.run("PBC", [&] { pbc.updateComplex(data, spc.groups); });
        bench.run("PBCEigen", [&] { pbc_eigen.updateComplex(data, spc.groups); });
    }
}

//----------------- IPBC Ewald -------------------

/**
 * Resize k-vectors according to current variables and box length
 */
void PolicyIonIonIPBC::updateBox(EwaldData &data, const Point &box) const {
    assert(data.policy == EwaldData::IPBC or data.policy == EwaldData::IPBCEigen);
    data.box_length = box;
    int ncc = std::ceil(data.invspace_cutoff);
    data.check_k2_zero = 0.1 * std::pow(2 * pc::pi / data.box_length.maxCoeff(), 2);
    int k_vector_size = (2 * ncc + 1) * (2 * ncc + 1) * (2 * ncc + 1) - 1;
    if (k_vector_size == 0) {
        data.k_vectors.resize(3, 1);
        data.Aks.resize(1);
        data.k_vectors.col(0) = Point(1, 0, 0); // Just so it is not the zero-vector
        data.Aks[0] = 0;
        data.num_kvectors = 1;
        data.Q_ion.resize(1);
        data.Q_dipole.resize(1);
    } else {
        double nc2 = data.invspace_cutoff * data.invspace_cutoff;
        data.k_vectors.resize(3, k_vector_size);
        data.Aks.resize(k_vector_size);
        data.num_kvectors = 0;
        data.k_vectors.setZero();
        data.Aks.setZero();
        int start_value = 0;
        for (int nx = 0; nx <= ncc; nx++) {
            auto dnx2 = double(nx * nx);
            double xfactor = (nx > 0) ? 2.0 : 1.0; // optimization of PBC Ewald
            for (int ny = -ncc * start_value; ny <= ncc; ny++) {
                auto dny2 = double(ny * ny);
                double yfactor = (ny > 0) ? 2.0 : 1.0; // optimization of PBC Ewald
                for (int nz = -ncc * start_value; nz <= ncc; nz++) {
                    double factor = xfactor * yfactor;
                    if (nz > 0) {
                        factor *= 2.0;
                    }
                    Point kv = 2.0 * pc::pi * Point(nx, ny, nz).cwiseQuotient(data.box_length);
                    const auto k2 = kv.squaredNorm() + data.kappa_squared; // last term is only for Yukawa-Ewald
                    if (k2 < data.check_k2_zero) {                         // Check if k2 != 0
                        continue;
                    }
                    if (data.use_spherical_sum) {
                        auto dnz2 = double(nz * nz);
                        if ((dnx2 + dny2 + dnz2) / nc2 > 1)
                            continue;
                        }
                    }
                    data.k_vectors.col(data.num_kvectors) = kv;
                    data.Aks[data.num_kvectors] = factor * exp(-k2 / (4 * data.alpha * data.alpha)) / k2;
                    data.num_kvectors++;
                }
            }
        }
        data.Q_ion.resize(data.num_kvectors);
        data.Q_dipole.resize(data.num_kvectors);
        data.Aks.conservativeResize(data.num_kvectors);
        data.k_vectors.conservativeResize(3, data.num_kvectors);
    }
}

void PolicyIonIonIPBC::updateComplex(EwaldData& d, const Space::GroupVector& groups) const {
    assert(d.policy == EwaldData::IPBC or d.policy == EwaldData::IPBCEigen);
    for (int k = 0; k < d.k_vectors.cols(); k++) {
        const Point &q = d.k_vectors.col(k);
        EwaldData::Tcomplex Q(0, 0);
        for (const auto& particle : groups | ranges::cpp20::views::join) {
            Q += q.cwiseProduct(particle.pos).array().cos().prod() * particle.charge; // see eq. 2 in doi:10/css8
        }
        d.Q_ion[k] = Q;
    }
}

void PolicyIonIonIPBCEigen::updateComplex(EwaldData& d, const Space::GroupVector& groups) const {
    assert(d.policy == EwaldData::IPBC or d.policy == EwaldData::IPBCEigen);
    auto [pos, charge] = mapGroupsToEigen(groups); // throws if inactive particles
    d.Q_ion.real() = (d.k_vectors.array().cwiseProduct(pos).array().cos().prod() * charge)
                         .colwise()
                         .sum(); // see eq. 2 in doi:10/css8
}

void PolicyIonIonIPBC::updateComplexOptimized(EwaldData& d, const Change& change, const Space::GroupVector& groups,
                                              const Space::GroupVector& oldgroups) const {
    assert(d.policy == EwaldData::IPBC or d.policy == EwaldData::IPBCEigen);
    assert(groups.size() == oldgroups.size());

    for (int k = 0; k < d.k_vectors.cols(); k++) {
        auto &Q = d.Q_ion[k];
        const Point &q = d.k_vectors.col(k);
        for (const auto& changed_group : change.groups) {
            const auto& g_new = groups.at(changed_group.group_index);
            const auto& g_old = oldgroups.at(changed_group.group_index);
            for (auto i : changed_group.relative_atom_indices) {
                if (i < g_new.size()) {
                    Q += q.cwiseProduct(g_new[i].pos).array().cos().prod() * g_new[i].charge;
                }
                if (i < g_old.size()) {
                    Q -= q.cwiseProduct(g_old[i].pos).array().cos().prod() * g_old[i].charge;
                }
            }
        }
    }
}

/**
 * @note The surface energy cannot be calculated for a partial change due to the squared `qr`.
 */
double PolicyIonIon::surfaceEnergy(const EwaldData& ewald_data, const Change& change,
                                   const Space::GroupVector& groups) {
    if (change.empty() || ewald_data.tinfoilSurrounding()) {
        return 0.0;
    }
    namespace rv = ranges::cpp20::views;
    auto charge_x_position = [](const Particle& particle) -> Point { return particle.charge * particle.pos; };
    auto qr_range = groups | rv::join | rv::transform(charge_x_position);
    const auto qr_squared = ranges::accumulate(qr_range, Point(0.0, 0.0, 0.0)).squaredNorm();

    return 2.0 * pc::pi / ((2.0 * ewald_data.surface_dielectric_constant + 1.0) * ewald_data.volume()) * qr_squared *
           ewald_data.bjerrum_length;
}

double PolicyIonIon::selfEnergy(const EwaldData& ewald_data, Change& change, Space::GroupVector& groups) {
    double charges_squared = 0.0;
    double charge_total = 0.0;
    if (change.matter_change) {
        for (const auto& changed_group : change.groups) {
            const auto& group = groups.at(changed_group.group_index);
            for (auto index : changed_group.relative_atom_indices) {
                if (index < group.size()) {
                    charges_squared += std::pow(group[index].charge, 2);
                    charge_total += group[index].charge;
                }
            }
        }
    } else if (change.everything and not change.volume_change) { // @todo shouldn't this be `and not` -> `or` ??
        for (const auto& particle : groups | ranges::cpp20::views::join) {
            charges_squared += particle.charge * particle.charge;
            charge_total += particle.charge;
        }
    }
    return selfEnergyFromChargeSums(ewald_data, charges_squared, charge_total); // depends on volume
}

double PolicyIonIon::selfEnergyFromChargeSums(const EwaldData& d, double charges_squared, double charge_total) const {
    auto Vcc = -pc::pi / 2.0 / d.alpha / d.alpha / d.volume() * charge_total *
               charge_total; // compensate with neutralizing background (if non-zero total charge in system)
    const auto beta = d.kappa / (2.0 * d.alpha);
    if (beta > 1e-6) {
        Vcc *= (1.0 - exp(-beta * beta)) / beta / beta;
    } // same as above but for Yukawa-systems
    return (-d.alpha * charges_squared / sqrt(pc::pi) * (exp(-beta * beta) + sqrt(pc::pi) * beta * erf(beta)) + Vcc) *
           d.bjerrum_length;
}

/**
 * Updates the reciprocal space terms 'Q^q' and 'A_k'.
 * See eqs. 24 and 25 in ref. for PBC Ewald, and eq. 2 in doi:10/css8 for IPBC Ewald.
 */
double PolicyIonIon::reciprocalEnergy(const EwaldData& ewald_data) {
    double energy = 0.0;
    for (int k = 0; k < ewald_data.Q_ion.size(); k++) {
        energy += ewald_data.Aks[k] * std::norm(ewald_data.Q_ion[k]);
    }
    return 2.0 * pc::pi * energy * ewald_data.bjerrum_length / ewald_data.volume();
}

double PolicyIonIonEigen::reciprocalEnergy(const EwaldData& ewald_data) {
    double energy = ewald_data.Aks.cwiseProduct(ewald_data.Q_ion.cwiseAbs2()).sum();
    return 2.0 * pc::pi * ewald_data.bjerrum_length * energy / ewald_data.volume();
}

Ewald::Ewald(const Space& spc, const EwaldData& data)
    : spc(spc)
    , data(data) {
    name = "ewald";
    policy = EwaldPolicyBase::makePolicy(data.policy);
    citation_information = policy->cite;
    init();
}

Ewald::Ewald(const json& j, const Space& spc)
    : Ewald(spc, static_cast<EwaldData>(j)) {}

void Ewald::init() {
    policy->updateBox(data, spc.geometry.getLength());
    policy->updateComplex(data, spc.groups); // brute force. todo: be selective
}

/**
 * If `old_groups` have been set and if the change object is only partial, this will attempt
 * to perform a faster, partial update of the k-vectors. Otherwise perform a full (slower) update.
 */
void Ewald::updateState(const Change& change) {
    if (change) {
        if (!change.groups.empty() && old_groups && !change.everything && !change.volume_change) {
            policy->updateComplexOptimized(data, change, spc.groups, *old_groups); // partial update (fast)
        } else {                                                          // full update (slow)
            policy->updateBox(data, spc.geometry.getLength());
            policy->updateComplex(data, spc.groups);
        }
    }
}

/**
 * the selfEnergy() is omitted as this is added as a separate term in `Hamiltonian`
 * (The pair-potential is responsible for this)
 */
double Ewald::energy(const Change& change) {
    if (change) {
        return policy->surfaceEnergy(data, change, spc.groups) + policy->reciprocalEnergy(data);
    }
    return 0.0;
}

/**
 * @param forces Destination force vector
 *
 * Calculate forces from reciprocal space. Note that
 * the destination force vector will *not* be zeroed
 * before addition.
 */
void Ewald::force(std::vector<Point> &forces) {
    assert(forces.size() == spc.particles.size());
    const double volume = spc.geometry.getVolume();

    // Surface contribution
    Point total_dipole_moment = {0.0, 0.0, 0.0};
    for (auto& particle : spc.particles) {
        auto mu = particle.hasExtension() ? particle.getExt().mu * particle.getExt().mulen : Point(0, 0, 0);
        total_dipole_moment += particle.pos * particle.charge + mu;
    }

    auto force = forces.begin(); // iterator to force vector on first particle

    assert(data.k_vectors.cols() == data.Q_ion.size());
    data.Q_dipole.resize(data.Q_ion.size());

    for (auto& particle : spc.particles) { // loop over particles
        (*force) = total_dipole_moment * particle.charge / (2.0 * data.surface_dielectric_constant + 1.0);
        double mu_scalar = particle.hasExtension() ? particle.getExt().mulen : 0.0;
        std::complex<double> qmu(mu_scalar, particle.charge);
        for (size_t i = 0; i < data.k_vectors.cols(); i++) { // loop over k vectors
            std::complex<double> Q = data.Q_ion[i] + data.Q_dipole[i];
            double k_dot_r = data.k_vectors.col(i).dot(particle.pos);
            std::complex<double> expKri(std::cos(k_dot_r), std::sin(k_dot_r));
            std::complex<double> repart = expKri * qmu * std::conj(Q);
            (*force) += std::real(repart) * data.k_vectors.col(i) * data.Aks[i];
        }
        (*force) *= -4.0 * pc::pi / volume * data.bjerrum_length; // to units of kT/Angstrom^2
        force++;                                                  // advance to next force vector
    }
}

/**
 * When updating k-vectors an optimization can be performed if the old group positions
 * are known. Use this function to set a pointer to the old groups (usually from the accepted Space)
 */
void Ewald::setOldGroups(const Space::GroupVector& old_groups) { this->old_groups = &old_groups; }

void Ewald::sync(EnergyTerm* energybase, const Change& change) {
    if (auto* other = dynamic_cast<const Ewald*>(energybase)) {
        if (!old_groups && other->state == MonteCarloState::ACCEPTED) {
            setOldGroups(other->spc.groups);
        }
        data.sync(other->data, change);
    } else {
        throw std::runtime_error("sync error");
    }
}

void Ewald::to_json(json &j) const { j = data; }

TEST_CASE("[Faunus] Energy::Ewald") {
    EwaldData data(R"({
                "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "ncutoff": 11.0, "spherical_sum": true, "cutoff": 5.0})"_json);
    data.policy = EwaldData::PBC;

    auto copy_position = [](auto& pos, auto& particle) { particle.pos = pos; };

    Space space;
    SpaceFactory::makeNaCl(space, 1, R"( {"type": "cuboid", "length": 10} )"_json);
    PointVector positions = {{0, 0, 0}, {1, 0, 0}};
    space.updateParticles(positions.begin(), positions.end(), space.particles.begin(), copy_position);
    auto ewald = Ewald(space, data);
    Change change;
    change.everything = true;
    // reference energy from IonIonPolicy:
    const auto reference_energy = (0.0020943951023931952 + 0.21303063979675319) * data.bjerrum_length;
    const auto energy_change = -16.8380445846; // move first position (0,0,0) --> (0.1, 0.1, 0.1)

    SUBCASE("energy") { CHECK_EQ(ewald.energy(change), doctest::Approx(reference_energy)); }

    SUBCASE("update and restore (full update)") {
        positions[0] = {0.1, 0.1, 0.1};
        space.updateParticles(positions.begin(), positions.end(), space.particles.begin(), copy_position);
        ewald.updateState(change);
        CHECK_EQ(ewald.energy(change), doctest::Approx(103.7300260099));
        CHECK_EQ(ewald.energy(change) - reference_energy, doctest::Approx(energy_change));

        positions[0] = {0.0, 0.0, 0.0};
        space.updateParticles(positions.begin(), positions.end(), space.particles.begin(), copy_position);
        CHECK((ewald.energy(change) != doctest::Approx(reference_energy))); // no match since k-vectors not updated
        ewald.updateState(change);
        CHECK_EQ(ewald.energy(change), doctest::Approx(reference_energy)); // that's better
    }

    Space trial_space;
    SpaceFactory::makeNaCl(trial_space, 1, R"( {"type": "cuboid", "length": 10} )"_json);
    positions = {{0, 0, 0}, {1, 0, 0}};
    trial_space.updateParticles(positions.begin(), positions.end(), trial_space.particles.begin(), copy_position);
    auto trial_ewald = Ewald(trial_space, data);
    CHECK_EQ(trial_ewald.energy(change), doctest::Approx(reference_energy));

    auto check_energy_change = [&](const Change& change) { // perturb from (0,0,0) --> (0.1, 0.1, 0.1)
        positions[0] = {0.0, 0.0, 0.0};
        trial_space.updateParticles(positions.begin(), positions.end(), trial_space.particles.begin(), copy_position);
        trial_ewald.updateState(change);
        ewald.sync(&trial_ewald, change);
        space.sync(trial_space, change);
        const auto old_energy = trial_ewald.energy(change);

        positions[0] = {0.1, 0.1, 0.1};
        trial_space.updateParticles(positions.begin(), positions.end(), trial_space.particles.begin(), copy_position);
        trial_ewald.updateState(change);
        ewald.sync(&trial_ewald, change);
        space.sync(trial_space, change);
        const auto new_energy = trial_ewald.energy(change); // position (0.1,0.1,0.1)
        CHECK_EQ(energy_change, doctest::Approx(new_energy - old_energy));
    };

    trial_ewald.setOldGroups(space.groups); // give access to positions in accepted state

    SUBCASE("update and restore (trial ewald)") {
        change.everything = true;
        check_energy_change(change);
    }
    SUBCASE("update and restore (partial using `all`)") {
        change.clear();
        auto& group_change = change.groups.emplace_back();
        group_change.group_index = 0;
        group_change.all = true;
        check_energy_change(change);
    }
    SUBCASE("update and restore (partial using `relative_atom_indices`)") {
        change.clear();
        auto& group_change = change.groups.emplace_back();
        group_change.group_index = 0;
        group_change.all = false;
        group_change.relative_atom_indices.push_back(0);
        check_energy_change(change);
    }
}

double Example2D::energy(const Change&) {
    double s =
        1 + std::sin(2.0 * pc::pi * particle.x()) + std::cos(2.0 * pc::pi * particle.y()) * static_cast<double>(use_2d);
    s *= scale_energy;
    if (particle.x() >= -2.00 && particle.x() <= -1.25)
        return 1 * s;
    if (particle.x() >= -1.25 && particle.x() <= -0.25)
        return 2 * s;
    if (particle.x() >= -0.25 && particle.x() <= 0.75)
        return 3 * s;
    if (particle.x() >= 0.75 && particle.x() <= 1.75)
        return 4 * s;
    if (particle.x() >= 1.75 && particle.x() <= 2.00)
        return 5 * s;
    return 1e10;
}

Example2D::Example2D(const json& j, Space& spc) : particle(spc.particles.at(0).pos) {
    scale_energy = j.value("scale", 1.0);
    use_2d = j.value("2D", true);
    name = "Example2D";
}
void Example2D::to_json(json &j) const {
    j["scale"] = scale_energy;
    j["2D"] = use_2d;
}

double ContainerOverlap::energy(const Change& change) {
    if (change && spc.geometry.type != Geometry::Variant::CUBOID) { // no need to check in PBC systems
        // *all* groups
        if (change.volume_change or change.everything) {
            return energyOfAllGroups();
        }
        // *subset* of groups
        for (const auto& group_change : change.groups) {
            if (groupIsOutsideContainer(group_change)) {
                return pc::infty;
            }
        }
    }
    return 0.0; // all particle insie simulation container
}

/**
 * @return infinity if any active particle is outside; zero otherwise
 */
double ContainerOverlap::energyOfAllGroups() const {
    auto positions = spc.groups | ranges::cpp20::views::join | ranges::cpp20::views::transform(&Particle::pos);
    bool outside = std::any_of(positions.begin(), positions.end(),
                               [&](const auto& position) { return spc.geometry.collision(position); });
    return outside ? pc::infty : 0.0;
}

/**
 * @brief Check a single group based on change object
 * @return true is any particle in group is outside; false otherwise
 */
bool ContainerOverlap::groupIsOutsideContainer(const Change::GroupChange& group_change) const {
    const auto& group = spc.groups.at(group_change.group_index);
    // *all* atoms
    if (group_change.all) {
        return std::any_of(group.begin(), group.end(),
                           [&](auto const& particle) { return spc.geometry.collision(particle.pos); });
    }
    // *subset* of atoms
    for (const auto particle_index : group_change.relative_atom_indices) {
        if (particle_index < group.size()) { // condition due to speciation move?
            if (spc.geometry.collision((group.begin() + particle_index)->pos)) {
                return true;
            }
        }
    }
    return false; // no overlap
}
ContainerOverlap::ContainerOverlap(const Space& spc) : spc(spc) { name = "ContainerOverlap"; }

// ------------- Isobaric ---------------

Isobaric::Isobaric(const json& j, const Space& spc) : spc(spc) {
    name = "isobaric";
    citation_information = "doi:10/dhn4v6 alt. Frenkel & Smith, 2nd Ed (Eq. 5.4.13)";

    for (const auto& [key, conversion_factor] : pressure_units) {
        if (auto it = j.find(key); it != j.end()) {
            pressure = it->get<double>() * conversion_factor; // convert to internal units
            return;
        }
    }
    throw ConfigurationError("specify pressure");
}

void Isobaric::to_json(json& j) const {
    for (const auto& [key, conversion_factor] : pressure_units) {
        j[key] = pressure / conversion_factor;
    }
    roundJSON(j, 5);
}

/**
 * @brief Calculates the energy contribution p × V / kT - (N + 1) × ln(V)
 * @return energy in kT
 */
double Isobaric::energy(const Change& change) {
    if (change.volume_change || change.everything || change.matter_change) {
        auto group_is_active = [](const auto& group) { return !group.empty(); };
        auto count_particles = [](const auto& group) { return group.isAtomic() ? group.size() : 1; };
        auto particles_per_group = spc.groups | ranges::cpp20::views::filter(group_is_active) |
                                   ranges::cpp20::views::transform(count_particles);
        auto number_of_particles = std::accumulate(particles_per_group.begin(), particles_per_group.end(), 0);
        const auto volume = spc.geometry.getVolume();
        return pressure * volume - static_cast<double>(number_of_particles + 1) * std::log(volume);
    }
    return 0.0;
}

const std::map<std::string, double> Isobaric::pressure_units = {{"P/atm", 1.0_atm},
                                                                {"P/bar", 1.0_bar},
                                                                {"P/kT", 1.0_kT},
                                                                {"P/mM", 1.0_millimolar},
                                                                {"P/Pa", 1.0_Pa}}; // add more if you fancy...

TEST_CASE("Energy::Isobaric") {
    Space spc;
    CHECK_NOTHROW(Isobaric(json({{"P/atm", 0.0}}), spc));
    CHECK_NOTHROW(Isobaric(json({{"P/bar", 0.0}}), spc));
    CHECK_NOTHROW(Isobaric(json({{"P/kT", 0.0}}), spc));
    CHECK_NOTHROW(Isobaric(json({{"P/mM", 0.0}}), spc));
    CHECK_NOTHROW(Isobaric(json({{"P/Pa", 0.0}}), spc));
    CHECK_THROWS(Isobaric(json({{"P/unknown_unit", 0.0}}), spc));

    SUBCASE("to_json") {
        json j;
        Isobaric(json({{"P/atm", 0.5}}), spc).to_json(j);
        CHECK_EQ(j.at("P/atm").get<double>(), doctest::Approx(0.5));
        Isobaric(json({{"P/bar", 0.4}}), spc).to_json(j);
        CHECK_EQ(j.at("P/bar").get<double>(), doctest::Approx(0.4));
        Isobaric(json({{"P/kT", 0.3}}), spc).to_json(j);
        CHECK_EQ(j.at("P/kT").get<double>(), doctest::Approx(0.3));
        Isobaric(json({{"P/mM", 0.2}}), spc).to_json(j);
        CHECK_EQ(j.at("P/mM").get<double>(), doctest::Approx(0.2));
        Isobaric(json({{"P/Pa", 0.1}}), spc).to_json(j);
        CHECK_EQ(j.at("P/Pa").get<double>(), doctest::Approx(0.1));
    }
}

Constrain::Constrain(const json& j, Space& spc) {
    name = "constrain";
    type = j.at("type").get<std::string>();
    if (const auto it = j.find("harmonic"); it != j.end()) {
        harmonic = std::make_optional<Faunus::pairpotential::HarmonicBond>();
        harmonic->from_json(*it);
    }
    coordinate = ReactionCoordinate::createReactionCoordinate({{type, j}}, spc);
}

double Constrain::energy(const Change& change) {
    if (change) {
        const auto value = (*coordinate)(); // calculate reaction coordinate
        if (harmonic) {
            return harmonic->half_force_constant * std::pow(harmonic->equilibrium_distance - value, 2);
        }
        if (not coordinate->inRange(value)) { // if outside allowed range ...
            return pc::infty;                 // ... return infinite energy
        }
    }
    return 0.0;
}
void Constrain::to_json(json& j) const {
    j = json(*coordinate).at(type);
    j.erase("resolution");
    j["type"] = type;
    if (harmonic) {
        harmonic->to_json(j["harmonic"]);
        j["harmonic"].erase("index");
        j.erase("range");
    }
}

void Bonded::updateGroupBonds(const Space::GroupType& group) {
    const auto first_particle_index = spc.getFirstParticleIndex(group);
    const auto group_index = spc.getGroupIndex(group);
    auto& bonds = internal_bonds[group_index];              // access or insert
    for (const auto& generic_bond : group.traits().bonds) { // generic bonds defined in topology
        const auto& bond = bonds.push_back(generic_bond->clone());
        bond->shiftIndices(first_particle_index); // shift to absolute particle index
        bond->setEnergyFunction(spc.particles);
    }
}

/**
 * Ensures that all internal bonds are updated according to bonds defined in the topology.
 */
void Bonded::updateInternalBonds() {
    internal_bonds.clear();
    std::for_each(spc.groups.begin(), spc.groups.end(), [&](auto& group) { updateGroupBonds(group); });
}

double Bonded::sumBondEnergy(const Bonded::BondVector& bonds) const {
#if (defined(__clang__) && __clang_major__ >= 10) || (defined(__GNUC__) && __GNUC__ >= 10)
    auto bond_energy = [&](const auto& bond) { return bond->energyFunc(spc.geometry.getDistanceFunc()); };
    return std::transform_reduce(bonds.begin(), bonds.end(), 0.0, std::plus<>(), bond_energy);
#else
    double energy = 0.0;
    for (const auto& bond : bonds) {
        energy += bond->energyFunc(spc.geometry.getDistanceFunc());
    }
    return energy;
#endif
}

Bonded::Bonded(const Space& spc, BondVector external_bonds = BondVector())
    : spc(spc), external_bonds(std::move(external_bonds)) {
    name = "bonded";
    updateInternalBonds();
    for (auto& bond : this->external_bonds) {
        bond->setEnergyFunction(spc.particles);
    }
}

Bonded::Bonded(const json& j, const Space& spc) : Bonded(spc, j.value("bondlist", BondVector())) {}

void Bonded::to_json(json& j) const {
    if (!external_bonds.empty()) {
        j["bondlist"] = external_bonds;
    }
    if (!internal_bonds.empty()) {
        json& array_of_bonds = j["bondlist-intramolecular"] = json::array();
        for ([[maybe_unused]] const auto& [group_index, bonds] : internal_bonds) {
            std::for_each(bonds.begin(), bonds.end(), [&](auto& bond) { array_of_bonds.push_back(bond); });
        }
    }
}

double Bonded::energy(const Change& change) {
    double energy = 0.0;
    if (change) {
        energy += sumBondEnergy(external_bonds);
        if (change.everything || change.volume_change) { // calc. for everything!
            for (const auto& [group_index, bonds] : internal_bonds) {
                if (!spc.groups.at(group_index).empty()) {
                    energy += sumBondEnergy(bonds);
                }
            }
        } else { // calc. for a subset of groups
            for (const auto& group_change : change.groups) {
                energy += internalGroupEnergy(group_change);
            }
        }
    }
    return energy;
}

double Bonded::internalGroupEnergy(const Change::GroupChange& changed) {
    using namespace ranges::cpp20::views; // @todo cpp20 --> std::ranges
    double energy = 0.0;
    const auto& group = spc.groups.at(changed.group_index);
    if (changed.internal && !group.empty()) {
        const auto& bonds = internal_bonds.at(changed.group_index);
        if (changed.all) { // all internal positions updated
            energy += sumBondEnergy(bonds);
        } else { // only partial update of affected atoms
            const auto first_particle_index = spc.getFirstParticleIndex(group);
            auto particle_indices = changed.relative_atom_indices |
                                    transform([first_particle_index](auto i) { return i + first_particle_index; });
            energy += sumEnergy(bonds, particle_indices);
        }
    }
    return energy;
}

/**
 * @param forces Target force vector for *all* particles in the system
 *
 * Each element in `force` represents the force on a particle and this
 * updates (add) the bonded force.
 *
 * - loop over groups and their internal bonds and _add_ force
 * - loop over inter-molecular bonds and _add_ force
 *
 * Force unit: kT/Å
 *
 * @warning Untested
 */
void Bonded::force(std::vector<Point>& forces) {
    auto distance_function = spc.geometry.getDistanceFunc();

    auto calculateForces = [&](const auto& bond) {
        if (!bond->hasForceFunction()) {
            throw std::runtime_error("force not implemented!");
        }
        for (const auto& [index, force] : bond->forceFunc(distance_function)) {
            forces.at(index) += force;
        }
    };

    for ([[maybe_unused]] const auto& [group_index, bonds] : internal_bonds) {
        for (const auto& bond : bonds) {
            calculateForces(bond);
        }
    }

    for (const auto& bond : external_bonds) {
        calculateForces(bond);
    }
}

//---------- Hamiltonian ------------

void Hamiltonian::to_json(json &j) const {
    std::for_each(energy_terms.cbegin(), energy_terms.cend(), [&](auto energy) { j.push_back(*energy); });
}

void Hamiltonian::addEwald(const json& j, Space& spc) {
    // note this will currently not detect multipolar energies ore
    // deeply nested "coulomb" pair-potentials
    json _j;
    if (j.count("default") == 1) { // try to detect FunctorPotential
        for (const auto& i : j["default"]) {
            if (i.count("coulomb") == 1) {
                _j = i["coulomb"];
                break;
            }
        }
    } else if (j.count("coulomb") == 1) {
        _j = j["coulomb"];
    } else {
        return;
    }

    if (_j.count("type") == 1) {
        if (_j.at("type") == "ewald") {
            if (_j.at("ewaldscheme") == "METALSLIT") {
                emplace_back<Energy::MetalSlitEwald>(_j, spc);
            } else {
                emplace_back<Energy::Ewald>(_j, spc);
            }
            faunus_logger->debug("hamiltonian expanded with {}", energy_terms.back()->name);
        }
    }
}

void Hamiltonian::force(PointVector& forces) {
    std::for_each(energy_terms.begin(), energy_terms.end(), [&](auto energy) { energy->force(forces); });
}

/**
 * @todo Move addEwald to Nonbonded as it now has access to Hamiltonian and can add
 */
Hamiltonian::Hamiltonian(Space& spc, const json& j) : energy_terms(this->vec) {
    name = "hamiltonian";
    if (!j.is_array()) {
        throw ConfigurationError("energy: json array expected");
    }

    // add container overlap energy for non-cuboidal geometries
    if (spc.geometry.type != Geometry::Variant::CUBOID) {
        emplace_back<Energy::ContainerOverlap>(spc);
        faunus_logger->debug("hamiltonian expanded with {}", energy_terms.back()->name);
    }

    for (const auto& j_energy : j) { // loop over energy list
        try {
            const auto& [key, value] = Faunus::jsonSingleItem(j_energy);
            try {
                if (key == "maxenergy") {
                    // looks like an unfortumate json scheme decision that requires special handling here
                    maximum_allowed_energy = value.get<double>();
                } else {
                    energy_terms.push_back(createEnergy(spc, key, value));
                    faunus_logger->debug("hamiltonian expanded with {}", key);
                    addEwald(value, spc); // add reciprocal Ewald terms if appropriate
                }
            } catch (std::exception& e) {
                usageTip.pick(key);
                throw ConfigurationError("{} -> {}", key, e.what());
            }
        } catch (std::exception& e) {
            throw ConfigurationError("energy -> {}", e.what()).attachJson(j_energy);
        }
    }
    latest_energies.reserve(energy_terms.size());
    checkBondedMolecules();
}

void Hamiltonian::checkBondedMolecules() const {
    if (find<Energy::Bonded>().empty()) { // no bond potential added? issue warning if molecules w. bonds
        auto molecules_with_bonds = Faunus::molecules | ranges::cpp20::views::filter([](const auto& molecule) {
                                        return !molecule.bonds.empty();
                                    });
        for (const auto& molecule : molecules_with_bonds) {
            faunus_logger->warn("{} bonds specified in topology but missing in energy", molecule.name);
        }
    }
}

double Hamiltonian::energy(const Change& change) {
    latest_energies.clear();
    for (auto& energy_ptr : energy_terms) {
        energy_ptr->state = state; // is this needed?
        energy_ptr->timer.start();
        const auto energy = energy_ptr->energy(change);
        latest_energies.push_back(energy);
        energy_ptr->timer.stop();
        if (energy >= maximum_allowed_energy || std::isnan(energy)) {
            break; // stop summing energies
        }
    }
    return std::accumulate(latest_energies.begin(), latest_energies.end(), 0.0);
}
void Hamiltonian::init() {
    std::for_each(energy_terms.begin(), energy_terms.end(), [&](auto& energy) { energy->init(); });
}

void Hamiltonian::updateState(const Change& change) {
    std::for_each(energy_terms.begin(), energy_terms.end(), [&](auto& energy) { energy->updateState(change); });
}

void Hamiltonian::sync(EnergyTerm* other_hamiltonian, const Change& change) {
    if (auto* other = dynamic_cast<Hamiltonian*>(other_hamiltonian)) {
        if (other->size() == size()) {
            latest_energies = other->latestEnergies();
            auto other_energy_term = other->energy_terms.cbegin();
            std::for_each(energy_terms.begin(), energy_terms.end(), [&](auto& energy_term) {
                energy_term->sync(other_energy_term->get(), change);
                std::advance(other_energy_term, 1);
            });
            return;
        }
    }
    throw std::runtime_error("hamiltonian mismatch");
}

/**
 * @brief Factory function to generate energy instances based on their name and json input
 * @param spc Space to use
 * @param name Name of energy
 * @param j Configuration for energy
 * @return shared pointer to energy base class
 * @throw if unknown name or configuration error
 *
 * New energy terms should be added to the if-else chain in the function
 */
std::unique_ptr<EnergyTerm> Hamiltonian::createEnergy(Space& spc, const std::string& name, const json& j) {
    using namespace pairpotential;
    using CoulombLJ = CombinedPairPotential<NewCoulombGalore, LennardJones>;
    using CoulombWCA = CombinedPairPotential<NewCoulombGalore, WeeksChandlerAndersen>;
    using PrimitiveModelWCA = CombinedPairPotential<Coulomb, WeeksChandlerAndersen>;
    using PrimitiveModel = CombinedPairPotential<Coulomb, HardSphere>;

    // only a single pairing policy and cutoff scheme so far
    using PairingPolicy = GroupPairing<GroupPairingPolicy<GroupCutoff>>;

    try {
        if (name == "nonbonded_coulomblj" || name == "nonbonded_newcoulomblj") {
            return std::make_unique<Nonbonded<PairEnergy<CoulombLJ, false>, PairingPolicy>>(j, spc, *this);
        }
        if (name == "nonbonded_coulomblj_EM") {
            return std::make_unique<NonbondedCached<PairEnergy<CoulombLJ, false>, PairingPolicy>>(j, spc, *this);
        }
        if (name == "nonbonded_splined") {
            return std::make_unique<Nonbonded<PairEnergy<pairpotential::SplinedPotential, false>, PairingPolicy>>(j, spc,
                                                                                                              *this);
        }
        if (name == "nonbonded" || name == "nonbonded_exact") {
            return std::make_unique<Nonbonded<PairEnergy<pairpotential::FunctorPotential, true>, PairingPolicy>>(j, spc,
                                                                                                             *this);
        }
        if (name == "nonbonded_cached") {
            return std::make_unique<NonbondedCached<PairEnergy<pairpotential::SplinedPotential>, PairingPolicy>>(j, spc,
                                                                                                             *this);
        }
        if (name == "nonbonded_coulombwca") {
            return std::make_unique<Nonbonded<PairEnergy<CoulombWCA, false>, PairingPolicy>>(j, spc, *this);
        }
        if (name == "nonbonded_pm" || name == "nonbonded_coulombhs") {
            return std::make_unique<Nonbonded<PairEnergy<PrimitiveModel, false>, PairingPolicy>>(j, spc, *this);
        }
        if (name == "nonbonded_pmwca") {
            return std::make_unique<Nonbonded<PairEnergy<PrimitiveModelWCA, false>, PairingPolicy>>(j, spc, *this);
        }
        if (name == "bonded") {
            return std::make_unique<Bonded>(j, spc);
        }
        if (name == "customexternal") {
            return std::make_unique<CustomExternal>(j, spc);
        }
        if (name == "akesson") {
            return std::make_unique<ExternalAkesson>(j, spc);
        }
        if (name == "confine") {
            return std::make_unique<Confine>(j, spc);
        }
        if (name == "constrain") {
            return std::make_unique<Constrain>(j, spc);
        }
        if (name == "custom-groupgroup") {
            return std::make_unique<CustomGroupGroup>(j, spc);
        }
        if (name == "example2d") {
            return std::make_unique<Example2D>(j, spc);
        }
        if (name == "isobaric") {
            return std::make_unique<Isobaric>(j, spc);
        }
        if (name == "penalty") {
#ifdef ENABLE_MPI
            return std::make_unique<PenaltyMPI>(j, spc, MPI::mpi);
#else
            return std::make_unique<Penalty>(j, spc);
#endif
        }
        if (name == "freesasa") {
#if defined ENABLE_FREESASA
            return std::make_unique<FreeSASAEnergy>(j, spc);
#else
            throw ConfigurationError("faunus not compiled with sasa support");
#endif
        }
        if (name == "sasa_reference") {
            return std::make_unique<SASAEnergyReference>(j, spc);
        }
        if (name == "sasa") {
            return std::make_unique<SASAEnergy>(j, spc);
        }
        throw ConfigurationError("'{}' unknown", name);
    } catch (std::exception& e) {
        // @todo For unknown reasons, ConfigurationError displays an *empty*
        //       what() message, which is why std::runtime_error is used instead.
        throw std::runtime_error("error creating energy -> "s + e.what());
    }
}

const std::vector<double>& Hamiltonian::latestEnergies() const { return latest_energies; }

#ifdef ENABLE_FREESASA

FreeSASAEnergy::FreeSASAEnergy(const Space& spc, const double cosolute_molarity, const double probe_radius)
    : spc(spc), cosolute_molarity(cosolute_molarity),
      parameters(std::make_unique<freesasa_parameters_fwd>(freesasa_default_parameters)) {
    name = "sasa";
    citation_information = "doi:10.12688/f1000research.7931.1";
    parameters->probe_radius = probe_radius;
    init();
    if (spc.geometry.asSimpleGeometry()->boundary_conditions.isPeriodic().count() != 0) {
        faunus_logger->error("PBC unsupported by FreeSASA; expect unphysical results");
    }
}

FreeSASAEnergy::FreeSASAEnergy(const json& j, const Space& spc)
    : FreeSASAEnergy(spc, j.value("molarity", 0.0) * 1.0_molar, j.value("radius", 1.4) * 1.0_angstrom) {}

void FreeSASAEnergy::updateSASA(const Change& change) {
    const auto particles = spc.activeParticles();
    const auto number_of_active_particles = std::distance(particles.begin(), particles.end());
    updateRadii(particles.begin(), particles.end(), change);
    updatePositions(particles.begin(), particles.end(), change);

    auto* result = freesasa_calc_coord(positions.data(), radii.data(), number_of_active_particles, parameters.get());
    if (result != nullptr && result->n_atoms == number_of_active_particles) {
        sasa.clear();
        sasa.reserve(number_of_active_particles);
        sasa.insert(sasa.begin(), result->sasa, result->sasa + result->n_atoms); // copy
        freesasa_result_free(result);
        assert(sasa.size() == number_of_active_particles);
    } else {
        throw std::runtime_error("FreeSASA failed");
    }
}

void FreeSASAEnergy::init() {
    Change change;
    change.everything = true;
    updateSASA(change);
}

double FreeSASAEnergy::energy(const Change& change) {
    double energy = 0.0;
    double surface_area = 0.0;
    updateSASA(change);
    for (const auto& [area, particle] : ranges::views::zip(sasa, spc.activeParticles())) {
        surface_area += area;
        energy += area * (particle.traits().tension + cosolute_molarity * particle.traits().tfe);
    }
    mean_surface_area += surface_area; // sample average area for accepted confs.
    return energy;
}

void FreeSASAEnergy::sync(EnergyTerm* energybase_ptr, const Change& change) {
    // since the full SASA is calculated in each energy evaluation, this is
    // currently not needed.
    if (auto* other = dynamic_cast<FreeSASAEnergy*>(energybase_ptr); other != nullptr) {
        sasa = other->sasa;
        if (change.everything || change.volume_change || change.matter_change) {
            radii = other->radii;
            positions = other->positions;
        } else {
            for (const auto& group_change : change.groups) {
                const auto& group = spc.groups.at(group_change.group_index);
                const auto offset = spc.getFirstActiveParticleIndex(group);
                auto absolute_atom_index = group_change.relative_atom_indices |
                                           ranges::cpp20::views::transform([offset](auto i) { return i + offset; });
                for (auto i : absolute_atom_index) {
                    radii.at(i) = other->radii.at(i);
                    for (size_t k = 0; k < 3; ++k) {
                        positions.at(3 * i + k) = other->positions.at(3 * i + k);
                    }
                }
            }
        }
    }
}

void FreeSASAEnergy::to_json(json& j) const {
    j["molarity"] = cosolute_molarity / 1.0_molar;
    j["radius"] = parameters->probe_radius / 1.0_angstrom;
    j[unicode::bracket("SASA") + "/" + unicode::angstrom + unicode::squared] = mean_surface_area.avg() / 1.0_angstrom;
    roundJSON(j, 5); // set json output precision
}

TEST_CASE("[Faunus] FreeSASA") {
    using doctest::Approx;
    Change change; // change object telling that a full energy calculation
    change.everything = true;
    pc::temperature = 300.0_K;
    atoms = R"([
        { "A": { "sigma": 4.0, "tfe": 1.0 } },
        { "B": { "sigma": 2.4, "tfe": 1.0 } }
    ])"_json.get<decltype(atoms)>();
    molecules = R"([
        { "M": { "atoms": ["A", "B"], "atomic": true } }
    ])"_json.get<decltype(molecules)>();
    json j = R"({
        "geometry": {"type": "sphere", "radius": 100 },
        "insertmolecules": [ { "M": { "N": 1 } } ]
    })"_json;
    Space spc = j;
    spc.particles.at(0).pos = {0.0, 0.0, 0.0};
    spc.particles.at(1).pos = {0.0, 0.0, 20.0};

    SUBCASE("Separated atoms") {
        FreeSASAEnergy sasa(spc, 1.5_molar, 1.4_angstrom);
        CHECK_EQ(sasa.energy(change), Approx(4.0 * pc::pi * (3.4 * 3.4 + 2.6 * 2.6) * 1.5 * 1.0_kJmol));
    }

    SUBCASE("Intersecting atoms") {
        FreeSASAEnergy sasa(spc, 1.5_molar, 1.4_angstrom);
        std::vector<std::pair<double, double>> distance_energy = {
            {0.0, 87.3576}, {2.5, 100.4612}, {5.0, 127.3487}, {7.5, 138.4422}, {10.0, 138.4422}};
        for (const auto& [distance, energy] : distance_energy) {
            spc.particles.at(1).pos = {0.0, 0.0, distance};
            CHECK_EQ(sasa.energy(change), Approx(energy).epsilon(0.02));
        }
    }

    SUBCASE("PBC") {}
}
#endif

/**
 * @param spc
 * @param cosolute_molarity in particles per angstrom cubed
 * @param probe_radius in angstrom
 * @param slices_per_atom number of slices of spheres in SASA calculation
 * @param dense_container flag specifying if a fast memory heavy version of cell_list container is used
 */
SASAEnergyReference::SASAEnergyReference(const Space& spc, double cosolute_molarity, double probe_radius,
                                         int slices_per_atom, bool dense_container)
    : spc(spc)
    , cosolute_molarity(cosolute_molarity) {
    using SASA::SASACellList;
    const auto periodic_dimensions = spc.geometry.asSimpleGeometry()->boundary_conditions.isPeriodic().count();
    switch (periodic_dimensions) {
    case 3: // PBC in all directions
        if (dense_container) {
            sasa = std::make_unique<SASACellList<SASA::DensePeriodicCellList>>(spc, probe_radius, slices_per_atom);
        } else {
            sasa = std::make_unique<SASACellList<SASA::SparsePeriodicCellList>>(spc, probe_radius, slices_per_atom);
        }
        break;
    case 0:
        if (dense_container) {
            sasa = std::make_unique<SASACellList<SASA::DenseFixedCellList>>(spc, probe_radius, slices_per_atom);
        } else {
            sasa = std::make_unique<SASACellList<SASA::SparseFixedCellList>>(spc, probe_radius, slices_per_atom);
        }
        break;
    default:
        sasa = std::make_unique<SASA::SASA>(spc, probe_radius, slices_per_atom);
        faunus_logger->warn("CellList neighbour search not available yet for current geometry");
        break;
    }
    name = "sasa";
    citation_information = "doi:10.12688/f1000research.7931.1";
    init();
}

SASAEnergyReference::SASAEnergyReference(const json& j, const Space& spc)
    : SASAEnergyReference(spc, j.at("molarity").get<double>() * 1.0_molar, j.value("radius", 1.4) * 1.0_angstrom,
                          j.value("slices", 25), j.value("dense", true)) {}

void SASAEnergyReference::init() {
    sasa->init(spc);
    areas.resize(spc.particles.size(), 0.0);
}

// -----------------------------------

/**
 * @param spc
 * @param cosolute_molarity in particles per angstrom cubed
 * @param probe_radius in angstrom
 * @param slices_per_atom number of slices of spheres in SASA calculation
 * @param dense_container flag specifying if a fast memory heavy version of cell_list container is used
 */
SASAEnergy::SASAEnergy(const Space& spc, double cosolute_molarity, double probe_radius, int slices_per_atom,
                       bool dense_container)
    : SASAEnergyReference(spc, cosolute_molarity, probe_radius, slices_per_atom, dense_container) {
    name = "sasa";
    citation_information = "doi:10.12688/f1000research.7931.1";
    init();
    areas = sasa->getAreas();
}

SASAEnergy::SASAEnergy(const json& j, const Space& spc)
    : SASAEnergy(spc, j.at("molarity").get<double>() * 1.0_molar, j.value("radius", 1.4) * 1.0_angstrom,
                 j.value("slices", 25), j.value("dense", true)) {}

void SASAEnergy::init() {
    areas.resize(spc.particles.size(), 0.0);
    current_neighbours.resize(spc.particles.size());
    changed_indices.reserve(spc.particles.size());
}

double SASAEnergyReference::energy(const Change& change) {
    double energy = 0.;
    if (state != MonteCarloState::ACCEPTED) {
        sasa->update(spc, change);
    }
    if (!change.everything) {
        sasa->needs_syncing = true;
    }
    const auto particles = spc.activeParticles();

    std::vector<index_type> target_indices;
    auto to_index = [this](const auto& particle) { return indexOf(particle); };
    target_indices = particles | ranges::cpp20::views::transform(to_index) | ranges::to<std::vector>;

    const auto neighbours = sasa->calcNeighbourData(spc, target_indices);
    sasa->updateSASA(neighbours, target_indices);

    const auto& new_areas = sasa->getAreas();
    ranges::cpp20::for_each(target_indices,
                            [this, &new_areas](const auto index) { areas.at(index) = new_areas.at(index); });

    auto accumulate_energy = [this, &energy](const auto& particle) {
        energy += areas.at(indexOf(particle)) * (particle.traits().tension + cosolute_molarity * particle.traits().tfe);
    };
    ranges::cpp20::for_each(particles, accumulate_energy);
    return energy;
}

void SASAEnergyReference::sync(EnergyTerm* energybase_ptr, const Change& change) {
    if (auto* other = dynamic_cast<SASAEnergyReference*>(energybase_ptr)) {
        areas = other->areas;
        if (sasa->needs_syncing) {
            sasa->update(other->spc, change);
        }
        other->sasa->needs_syncing = false;
        sasa->needs_syncing = false;
    }
}

void SASAEnergyReference::to_json(json& j) const {
    j["molarity"] = cosolute_molarity / 1.0_molar;
    roundJSON(j, 6); // set json output precision
}

const std::vector<double>& SASAEnergyReference::getAreas() const { return areas; }

/**
 * @brief Finds absolute indices of particles whose SASA has changed
 * @param change Change object
 */
void SASAEnergy::updateChangedIndices(const Change& change) {
    //!< if the state is ACCEPTED there is no need to recalculate SASAs so return empty target_indices
    if (state == MonteCarloState::ACCEPTED) {
        return;
    }

    //!< using set to get rid of repeating target_indices in touched_atoms neighbours lists
    std::set<index_type> target_indices;
    for (const auto& group_change : change.groups) {
        const auto& group = spc.groups.at(group_change.group_index);
        const auto offset = spc.getFirstParticleIndex(group);
        auto insert_changed = [this, &target_indices](const auto index) {
            target_indices.insert(index);
            insertChangedNeighboursOf(index, target_indices);
        };

        if (group_change.relative_atom_indices.empty()) {
            const auto indices = ranges::cpp20::views::iota(offset, group.size() + offset);
            ranges::cpp20::for_each(indices, insert_changed);
        } else {
            const auto indices = group_change.relative_atom_indices |
                                 ranges::cpp20::views::transform([offset](auto i) { return offset + i; });
            ranges::cpp20::for_each(indices, insert_changed);
        }
    }
    changed_indices.assign(target_indices.begin(), target_indices.end());
}
/**
 @brief
* for each touched atoms we need to recalculate SASA for the touched atom
* and for its neighbours before it moved and after it moved
*
*      * @param index_type index of particle whose neighbours sasa has changed
*      * @param target_indices placeholder to insert changed indices
**/
void SASAEnergy::insertChangedNeighboursOf(const index_type index, std::set<index_type>& target_indices) const {
    const auto& current_neighbour = sasa->calcNeighbourDataOfParticle(spc, index).indices;
    const auto& past_neighbour = current_neighbours.at(index);
    target_indices.insert(past_neighbour.begin(), past_neighbour.end());
    target_indices.insert(current_neighbour.begin(), current_neighbour.end());
}

double SASAEnergy::energy(const Change& change) {
    double energy(0.);

    //! update not needed when the state is already accepted
    if (state != MonteCarloState::ACCEPTED) {
        sasa->update(spc, change);
    }
    const auto particles = spc.activeParticles();
    changed_indices.clear();
    if (change.everything) { //! all the active particles will be used for SASA calculation
        auto to_index = [this](const auto& particle) { return indexOf(particle); };
        changed_indices = particles | ranges::cpp20::views::transform(to_index) | ranges::to<std::vector>;
    } else {
        updateChangedIndices(change);
        sasa->needs_syncing = true;
    }

    const auto neighbours_data = sasa->calcNeighbourData(spc, changed_indices);

    // update sasa areas in sasa object and update
    sasa->updateSASA(neighbours_data, changed_indices);
    const auto& new_areas = sasa->getAreas();

    for (const auto& [neighbour, index] : ranges::views::zip(neighbours_data, changed_indices)) {
        current_neighbours.at(index) = neighbour.indices;
        areas.at(index) = new_areas.at(index);
    }

    auto accumulate_energy = [this, &energy](const auto& particle) {
        energy += areas.at(indexOf(particle)) * (particle.traits().tension + cosolute_molarity * particle.traits().tfe);
    };
    ranges::cpp20::for_each(particles, accumulate_energy);
    return energy;
}

void SASAEnergy::sync(EnergyTerm* energybase_ptr, const Change& change) {
    if (auto* other = dynamic_cast<SASAEnergy*>(energybase_ptr)) {
        const auto sync_data = [this, other](const size_t changed_index) {
            this->current_neighbours.at(changed_index) = other->current_neighbours.at(changed_index);
            this->areas.at(changed_index) = other->areas.at(changed_index);
        };

        if (state == MonteCarloState::TRIAL) { //! changed_indices get updated only in TRIAL state for speedup
            ranges::for_each(this->changed_indices, sync_data);
        } else {
            ranges::for_each(other->changed_indices, sync_data);
        }

        if (sasa->needs_syncing) {
            sasa->update(other->spc, change);
        }
        other->sasa->needs_syncing = false;
        sasa->needs_syncing = false;
    }
}

TEST_CASE_TEMPLATE("[Faunus] SASAEnergy_updates", EnergyTemplate, SASAEnergyReference, SASAEnergy) {
    using doctest::Approx;
    pc::temperature = 300.0_K;
    atoms = R"([
        { "A": { "sigma": 2.0, "tfe": 1.0 } },
        { "B": { "sigma": 2.0, "tfe": 1.0 } }
    ])"_json.get<decltype(atoms)>();
    molecules = R"([
        { "M": { "structure": [
                { "A": [0.0, 0.0, 0.0] },
                { "B": [3.0, 0.0, 0.0] }
                ] }}
    ])"_json.get<decltype(molecules)>();

    Space spc;
    spc.geometry = R"( {"type": "cuboid", "length": 200} )"_json;
    json j = json::array();
    j.push_back({{"M", {{"N", 2}}}});
    InsertMoleculesInSpace::insertMolecules(j, spc);

    spc.particles.at(0).pos = {52.0, 0.0, 0.0};
    spc.particles.at(1).pos = {8.0, 0.0, 0.0};
    spc.particles.at(2).pos = {11.0, 0.0, 0.0};
    spc.particles.at(3).pos = {50.0, 0.0, 0.0};

    SUBCASE("Update everything") {
        Change change;
        change.everything = true;

        EnergyTemplate sasa_energy(spc, 1.5_molar, 1.0_angstrom, 20);
        FreeSASAEnergy ref_energy(spc, 1.5_molar, 1.0_angstrom);
        CHECK_EQ(sasa_energy.energy(change), Approx(98.0789260855));
        CHECK_EQ(ref_energy.energy(change), Approx(98.0789260855));

        for (const auto& [area, ref_area] : ranges::views::zip(sasa_energy.getAreas(), ref_energy.getAreas())) {
            CHECK_EQ(area, Approx(ref_area));
        }
    }

    SUBCASE("Partial update") {
        Change change;
        change.everything = true;
        EnergyTemplate sasa_energy(spc, 1.5_molar, 1.0_angstrom, 20);
        FreeSASAEnergy ref_energy(spc, 1.5_molar, 1.0_angstrom);

        // we must caclulate SASAs of all 4 particles before updating position of one of them
        // since all SASAs are needed in the final energy calculation.
        // these are usually obtained from previous MC step in simulations so here we need to calculate them explicitly
        sasa_energy.energy(change);
        ref_energy.energy(change);

        spc.particles.at(3).pos = {14.0, 0.0, 0.0}; // update last particle in last group
        change.everything = false;
        auto& changed_data = change.groups.emplace_back();
        changed_data.group_index = 1;
        changed_data.relative_atom_indices = {1};

        CHECK_EQ(sasa_energy.energy(change), Approx(105.7104501023));
        CHECK_EQ(ref_energy.energy(change), Approx(105.7104501023));

        for (const auto& [area, ref_area] : ranges::views::zip(sasa_energy.getAreas(), ref_energy.getAreas())) {
            CHECK_EQ(area, Approx(ref_area));
        }
    }
}

//==================== GroupCutoff ====================

GroupCutoff::GroupCutoff(Space::GeometryType& geometry) : geometry(geometry) {}

void GroupCutoff::setSingleCutoff(const double cutoff) {
    if (cutoff < std::sqrt(pc::max_value)) {
        default_cutoff_squared = cutoff * cutoff;
    } else {
        default_cutoff_squared = pc::max_value;
    }
    for (const auto& molecule1 : Faunus::molecules) {
        for (const auto& molecule2 : Faunus::molecules) {
            cutoff_squared.set(molecule1.id(), molecule2.id(), default_cutoff_squared);
        }
    }
}

double GroupCutoff::getCutoff(size_t id1, size_t id2) const {
    if (cutoff_squared.size() != Faunus::molecules.size() || id1 >= cutoff_squared.size() ||
        id2 >= cutoff_squared.size()) {
        throw std::out_of_range("cutoff matrix doesn't fit molecules");
    }
    return std::sqrt(cutoff_squared(id1, id2));
}

void from_json(const json& j, GroupCutoff& cutoff) {
    // default: no group-to-group cutoff
    cutoff.setSingleCutoff(std::sqrt(pc::max_value));

    if (const auto it = j.find("cutoff_g2g"); it != j.end()) {
        if (it->is_number()) {
            cutoff.setSingleCutoff(it->get<double>());
        } else if (it->is_object()) {
            cutoff.setSingleCutoff(it->value("default", pc::max_value));
            for (const auto& [named_pair, pair_cutoff] : it->items()) {
                if (named_pair == "default") {
                    continue;
                }
                try {
                    if (const auto molecules_names = splitConvert<std::string>(named_pair);
                        molecules_names.size() == 2) {
                        const auto& molecule1 = findMoleculeByName(molecules_names[0]);
                        const auto& molecule2 = findMoleculeByName(molecules_names[1]);
                        cutoff.cutoff_squared.set(molecule1.id(), molecule2.id(),
                                                  std::pow(pair_cutoff.get<double>(), 2));
                        faunus_logger->debug("custom cutoff for {}-{} = {} Å", molecule1.name, molecule2.name,
                                             pair_cutoff.get<double>());
                    } else {
                        throw std::runtime_error("invalid molecules names");
                    }
                } catch (const std::exception& e) {
                    throw ConfigurationError("Unable to set a custom cutoff for {}: {}", named_pair, e.what());
                }
            }
        }
    }
}

void to_json(json &j, const GroupCutoff &cutoff) {
    auto _j = json::object();
    for (auto &a : Faunus::molecules) {
        for (auto &b : Faunus::molecules) {
            if (a.id() >= b.id()) {
                if (not a.atomic && not b.atomic) {
                    auto cutoff_squared = cutoff.cutoff_squared(a.id(), b.id());
                    if (cutoff_squared < pc::max_value) {
                        _j[a.name + " " + b.name] = std::sqrt(cutoff_squared);
                    }
                }
            }
        }
    }
    if (not _j.empty()) {
        j["cutoff_g2g"] = _j;
    }
}

TEST_CASE("[Faunus] GroupCutoff") {
    Geometry::Chameleon geo;
    Faunus::atoms = R"([
        { "A": { "sigma": 4.0 } },
        { "B": { "sigma": 2.4 } }
    ])"_json.get<decltype(atoms)>();
    Faunus::molecules = R"([
        { "M": { "atoms": ["A", "B"] } },
        { "Q": { "atoms": ["A", "B"] } }
    ])"_json.get<decltype(molecules)>();

    SUBCASE("old style default") {
        GroupCutoff cutoff(geo);
        from_json(R"({"cutoff_g2g": 12.0})"_json, cutoff);
        CHECK_EQ(cutoff.getCutoff(0, 0), doctest::Approx(12.0));
        CHECK_EQ(cutoff.getCutoff(0, 1), doctest::Approx(12.0));
        CHECK_EQ(cutoff.getCutoff(1, 1), doctest::Approx(12.0));
    }

    SUBCASE("new style default") {
        GroupCutoff cutoff(geo);
        from_json(R"({"cutoff_g2g": {"default": 13.0}})"_json, cutoff);
        CHECK_EQ(cutoff.getCutoff(0, 0), doctest::Approx(13.0));
        CHECK_EQ(cutoff.getCutoff(0, 1), doctest::Approx(13.0));
        CHECK_EQ(cutoff.getCutoff(1, 1), doctest::Approx(13.0));
    }

    SUBCASE("custom") {
        GroupCutoff cutoff(geo);
        from_json(R"({"cutoff_g2g": {"default": 13.0, "M Q": 14.0}})"_json, cutoff);
        CHECK_EQ(cutoff.getCutoff(0, 0), doctest::Approx(13.0));
        CHECK_EQ(cutoff.getCutoff(0, 1), doctest::Approx(14.0));
        CHECK_EQ(cutoff.getCutoff(1, 1), doctest::Approx(13.0));
    }

    SUBCASE("custom - no default value") {
        GroupCutoff cutoff(geo);
        from_json(R"({"cutoff_g2g": {"M Q": 11.0}})"_json, cutoff);
        CHECK_EQ(cutoff.getCutoff(0, 0), doctest::Approx(std::sqrt(pc::max_value)));
        CHECK_EQ(cutoff.getCutoff(0, 1), doctest::Approx(11.0));
        CHECK_EQ(cutoff.getCutoff(1, 1), doctest::Approx(std::sqrt(pc::max_value)));
    }
}

EnergyAccumulatorBase::EnergyAccumulatorBase(double value) : value(value) {}

void EnergyAccumulatorBase::reserve([[maybe_unused]] size_t number_of_particles) {}

void EnergyAccumulatorBase::clear() { value = 0.0; }

EnergyAccumulatorBase::operator double() { return value; }

void EnergyAccumulatorBase::from_json(const json& j) {
    scheme = j.value("summation_policy", Scheme::SERIAL);
#ifndef HAS_PARALLEL_TRANSFORM_REDUCE
    if (scheme == Scheme::PARALLEL) {
        faunus_logger->warn("'parallel' summation unavailable; falling back to 'serial'");
        scheme = Scheme::SERIAL;
    }
#endif
#ifndef _OPENMP
    if (scheme == Scheme::OPENMP) {
        faunus_logger->warn("'openmp' summation unavailable; falling back to 'serial'");
        scheme = Scheme::SERIAL;
    }
#endif
    faunus_logger->debug("setting summation policy to {}", json(scheme).dump(1));
}

void EnergyAccumulatorBase::to_json(json& j) const { j["summation_policy"] = scheme; }

/*
* The only difference from ordinary Ewald is that the reciprocal energy is divided by 2
*/
double PolicyIonIonMetalSlit::reciprocalEnergy(const EwaldData& ewald_data) {
    return 0.5 * PolicyIonIon::reciprocalEnergy(ewald_data);
}

void PolicyIonIonMetalSlit::updateBox(EwaldData& ewald_data, const Point& box) const {
    // Double up z-direction to hold mirror charges
    //PolicyIonIon::updateBox(ewald_data, {box.x(), box.y(), 2.0 * box.z()});
    PolicyIonIon::updateBox(ewald_data, {box.x(), box.y(), box.z()});
}

void PolicyIonIonMetalSlit::updateComplex(EwaldData& ewald_data, const Space::GroupVector& groups) const {
    namespace rv = ranges::cpp20::views;
    auto mirror_z = getMirrorLambda(ewald_data);

    for (int k = 0; k < ewald_data.k_vectors.cols(); k++) {
        const Point& q = ewald_data.k_vectors.col(k);
        auto positions = groups | rv::join | rv::transform(&Particle::pos);
        auto mirror_positions = positions | rv::transform(mirror_z);
        auto charges = groups | rv::join | rv::transform(&Particle::charge);
        //Mirror charges should be inverted
        auto mirror_charges = charges | rv::transform([](auto charge){ return -charge; });
        // Operate on both original and mirrored positions (concatenated) and zip with charges
        ewald_data.Q_ion[k] = sumWavevector(q, ranges::views::zip(ranges::views::concat(positions, mirror_positions),
                                                                  ranges::views::concat(charges, mirror_charges)));
    }
}

/// @todo reimplement w. mirrored charges
double PolicyIonIonMetalSlit::surfaceEnergy(const EwaldData& ewald_data, const Change& change,
                                            const Space::GroupVector& groups) {
    return PolicyIonIon::surfaceEnergy(ewald_data, change, groups);
}

PolicyIonIonMetalSlit::PolicyIonIonMetalSlit()
    : PolicyIonIon() {}

/* This override adds all mirror charges */
EwaldData::Tcomplex PolicyIonIonMetalSlit::sumWavevectorGroups(const Point& wavevector, const Group& group,
                                                               const std::vector<Change::index_type>& indices,
                                                               EwaldData& ewald_data) const {
    namespace rv = ranges::cpp20::views;
    auto new_indices = indices | rv::filter([&](auto i) { return i < group.size(); }) | ranges::to_vector;
    auto mirror_z = getMirrorLambda(ewald_data);
    auto positions = group[new_indices] | rv::transform(&Particle::pos);
    auto charges = group[new_indices] | rv::transform(&Particle::charge);
    auto mirror_charges = charges | rv::transform([](auto charge){ return -charge; });
    auto mirror_positions = positions | rv::transform(mirror_z);

    return sumWavevector(wavevector, ranges::views::zip(ranges::views::concat(positions, mirror_positions),
                                                        ranges::views::concat(charges, mirror_charges)));
}

/*
* charge_total is always zero due to image charges
* Only the "real" self terms should be included, hence maybe no override is needed?
*/
double PolicyIonIonMetalSlit::selfEnergyFromChargeSums(const EwaldData& ewald_data, double charges_squared,
                                                       [[maybe_unused]] double charge_total) const {
    return PolicyIonIon::selfEnergyFromChargeSums(ewald_data, charges_squared, 0.0);
}

TEST_CASE("[Faunus] Ewald - PolicyIonIonMetalSlit") {
    using doctest::Approx;
    Space spc;
    spc.particles.resize(4);
    spc.geometry = R"( {"type": "slit", "length": [20, 20, 10]} )"_json;
    spc.particles.at(0) = R"( {"id": 0, "pos": [1,  0, 1], "q":  1.0} )"_json;
    spc.particles.at(1) = R"( {"id": 1, "pos": [0,  0, 3], "q": -1.0} )"_json;
    spc.particles.at(2) = R"( {"id": 2, "pos": [0,  4, 4], "q":  1.0} )"_json;
    spc.particles.at(3) = R"( {"id": 3, "pos": [0, -4, 2], "q":  1.0} )"_json;

    if (Faunus::molecules.empty()) {
        Faunus::molecules.resize(1);
    }
    Group g(0, spc.particles.begin(), spc.particles.end());
    spc.groups.push_back(g);

    /*auto data = static_cast<EwaldData>(R"({
                "epsr": 1.0, "alpha": 0.314129, "epss": 0.0,
                "ncutoff": 11.0, "spherical_sum": false, "cutoff": 10.0})"_json);*/
    EwaldData data(R"({
            "epsr": 1.0, "alpha": 0.314129, "epss": 0.0,
            "ncutoff": 11.0, "spherical_sum": false, "cutoff": 10.0})"_json);
    Change c;
    c.everything = true;
    data.policy = EwaldData::METALSLIT;
    auto box_length = spc.geometry.getLength();
    box_length.z() *= 2.0;
    PolicyIonIonMetalSlit policy;
    policy.updateBox(data, box_length);
    policy.updateComplex(data, spc.groups);

    auto jj = R"({"epsr": 1.0, "alpha": 0.314129, "epss": 0.0, "ncutoff": 11.0, "spherical_sum": false, "cutoff": 10.0})"_json;
    auto imageEwald = MetalSlitEwald(jj, spc);

    auto surfaceEnergy = policy.surfaceEnergy(data, c, spc.groups) / data.bjerrum_length;
    auto selfEnergy = policy.selfEnergy(data, c, spc.groups) / data.bjerrum_length;
    auto reciprocalEnergy = policy.reciprocalEnergy(data) / data.bjerrum_length;

    CHECK(surfaceEnergy == Approx(0.0));
    CHECK(selfEnergy == Approx(-0.7089132387610925));
    CHECK(reciprocalEnergy == Approx(0.17604271991696552));
    CHECK(imageEwald.completeMirrorEnergy() == Approx(-0.09724725961362456));
}

// -----------------------------------------

MetalSlitEwald::MetalSlitEwald(const json& j, const Space& spc)
    : Ewald(j, spc) {
    name = "metal slit " + name;

    auto jcopy = j;
    jcopy["type"] = "ewald";
    auto coulomb_galore = pairpotential::NewCoulombGalore();
    from_json(jcopy, coulomb_galore);
    pair_potential = coulomb_galore.getCoulombGalore();

    if (data.policy != EwaldData::METALSLIT) {
        faunus_logger->warn("{}: Invalid Ewald policy -- setting this to required METALSLIT", name);
        policy = std::make_unique<PolicyIonIonMetalSlit>();
    }
    auto box = spc.geometry.getLength();
    updateEnlargedGeometry(spc);
    policy->updateBox(data, box); // Ewald data must have enlarged box information
    //policy->updateBox(data, enlarged_geometry.getLength()); // Ewald data must have enlarged box information
    policy->updateComplex(data, spc.groups);
}

/**
 * The `enlarged_geometry` is set to 2 x the z-direction of the given Space
 */
void MetalSlitEwald::updateEnlargedGeometry(const Space& spc) {
    auto slit_dimensions = getSlitDimensions(spc);
    slit_dimensions.z() *= 2.0;
    enlarged_geometry.setLength(slit_dimensions);
}

/**
 * @return Dimensions of slit of the geometry of `spc`
 * @throw is the Space geometry is not Geometry::Slit
 */
Point MetalSlitEwald::getSlitDimensions(const Space& spc) {
    //if (auto slit_ptr = std::dynamic_pointer_cast<Geometry::Slit>(spc.geometry.asSimpleGeometry())) {
    if (auto slit_ptr = std::dynamic_pointer_cast<Geometry::Cuboid>(spc.geometry.asSimpleGeometry())) {
        return slit_ptr->getLength();
    }
    throw std::runtime_error("slit geometry required by metallic ewald");
}

/**
 * Adds the reciprocal energy as well as real <-> mirror charges
 */
double MetalSlitEwald::energy(const Change& change) {
    // optimize if only a single particle has been updated
    if (auto index_pair = change.singleParticleChange()) {
        const auto& particle = spc.groups.at(index_pair->first).at(index_pair->second);
        return Ewald::energy(change) + singleParticleMirrorEnergy(particle);
    }

    if (!change.everything) {
        faunus_logger->debug("{}: partial particle update not specialized - this may be slow", name);
    }
    return Ewald::energy(change) + completeMirrorEnergy();
}

/**
 * @todo optimize by using change object (currently everything is recalculated -> expensive)
 */
double MetalSlitEwald::completeMirrorEnergy() const {
    namespace rv = ranges::cpp20::views;
    double energy = 0.0;

    auto mirror_z = std::dynamic_pointer_cast<PolicyIonIonMetalSlit>(policy)->getMirrorLambda(data);
    auto all_particles = spc.groups | rv::join;
    auto mirror_positions = all_particles | rv::transform(&Particle::pos) | rv::transform(mirror_z);
    auto mirror_charges = all_particles | rv::transform(&Particle::charge) | rv::transform([](auto charge){ return -charge; });

    for (const auto [mirror_pos, mirror_charge] : ranges::views::zip(mirror_positions, mirror_charges)) {
        for (const auto& particle : all_particles) {
            double distance = enlarged_geometry.vdist(particle.pos, mirror_pos).norm();
            energy += pair_potential.ion_ion_energy(mirror_charge, particle.charge, distance);
        }
    }
    return 0.5 * energy;
}

/**
 * Calculates the energy of a single particle:
 *
 * - Particle <-> all mirror charges
 * - all other particles <-> mirror charge of particle
 *
 * @param particle Changed particle
 * @return Energy in kT
 */
double MetalSlitEwald::singleParticleMirrorEnergy(const Particle& particle) const {
    namespace rv = ranges::cpp20::views;
    double energy = 0.0;
    auto mirror_z = std::dynamic_pointer_cast<PolicyIonIonMetalSlit>(policy)->getMirrorLambda(data);
    auto all_particles = spc.groups | rv::join;
    //auto mirror_positions = all_particles | rv::transform(&Particle::pos) | rv::transform(mirror_z);
    //auto charges = all_particles | rv::transform(&Particle::charge);
    //auto mirror_charges = charges | rv::transform([](auto charge){ return -charge; });

    // 
    /*for (const auto [mirror_position, mirror_charge] : ranges::views::zip(mirror_positions, mirror_charges)) {
        const auto distance = enlarged_geometry.vdist(particle.pos, mirror_position).norm();
        energy += pair_potential.ion_ion_energy(particle.charge, mirror_charge, distance);
    }*/

    // Only need to calculate interaction between all other particles with new mirror charge from updated particle, 
    // since this is equal to particle with all mirror charges
    const auto mirror_position = mirror_z(particle.pos);
    const auto mirror_charge = -particle.charge;
    for (const auto& other_particle : all_particles) {
        /*if (&other_particle == &particle) {
            continue; // interaction w. self reflection already captured in above loop
        }*/
        const auto distance = enlarged_geometry.vdist(other_particle.pos, mirror_position).norm();
        energy += pair_potential.ion_ion_energy(other_particle.charge, mirror_charge, distance);
    }

    //Interaction with the self image should be divided by 2
    const auto dist = enlarged_geometry.vdist(particle.pos, mirror_position).norm();
    energy -= 0.5 * pair_potential.ion_ion_energy(particle.charge, mirror_charge, dist);

    return energy;
}

CustomGroupGroup::CustomGroupGroup(const json& j, const Space& spc)
    : spc(spc)
    , json_input_backup(j) {
    name = "custom-groupgroup";

    const auto& molecule1 = findMoleculeByName(j.at("name1").get<std::string>());
    const auto& molecule2 = findMoleculeByName(j.at("name2").get<std::string>());
    if (molecule1.isAtomic() || molecule2.isAtomic()) {
        throw ConfigurationError("molecular groups required");
    }
    molid1 = molecule1.id();
    molid2 = molecule2.id();

    auto& constants = json_input_backup["constants"];
    if (constants == nullptr) {
        constants = json::object();
    }
    constants["e0"] = pc::vacuum_permittivity;
    constants["kB"] = pc::boltzmann_constant;
    constants["kT"] = pc::kT();
    constants["Nav"] = pc::avogadro;
    constants["T"] = pc::temperature;
    expr = std::make_unique<ExprFunction<double>>();
    expr->set(json_input_backup, {{"R", &properties.mass_center_separation},
                                  {"Z1", &properties.mean_charge1},
                                  {"Z2", &properties.mean_charge2}});
}

double CustomGroupGroup::energy([[maybe_unused]] const Change& change) {
    // matches active groups with either molid1 or molid2
    auto match_groups = [&](const auto& group) -> bool {
        return group.isFull() & ((group.id == molid1) | (group.id == molid2));
    };

    auto group_group_energy = [&](auto index1, auto index2) -> double {
        const Group& group1 = spc.groups[index1];
        const Group& group2 = spc.groups[index2];
        if (((group1.id == molid1) & (group2.id == molid2)) | ((group1.id == molid2) & (group2.id == molid1))) {
            setParameters(group1, group2);
            return expr->operator()();
        }
        return 0.0;
    };

    // all indices matching either molid1 or molid2
    auto indices = spc.groups | ranges::cpp20::views::filter(match_groups) |
                   ranges::cpp20::views::transform([&](const auto& group) { return spc.getGroupIndex(group); }) |
                   ranges::to_vector;

    return for_each_unique_pair(indices.begin(), indices.end(), group_group_energy, std::plus<>());
}

void CustomGroupGroup::setParameters(const Group& group1, const Group& group2) {
    auto& mean_charge1 = mean_charges[group1.id];
    auto& mean_charge2 = mean_charges[group2.id];
    mean_charge1 += monopoleMoment(group1.begin(), group1.end());
    mean_charge2 += monopoleMoment(group2.begin(), group2.end());
    properties.mean_charge1 = mean_charge1.avg();
    properties.mean_charge2 = mean_charge2.avg();
    properties.mass_center_separation = sqrt(spc.geometry.sqdist(group1.mass_center, group2.mass_center));
}

void CustomGroupGroup::to_json(json& j) const {
    j = json_input_backup;
    j["mean_charges"] = {{"Z1", mean_charges.at(molid1).avg()}, {"Z2", mean_charges.at(molid2).avg()}};
};

} // end of namespace Faunus::Energy
