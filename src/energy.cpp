#include "energy.h"
#include "penalty.h"
#include "potentials.h"
#include "externalpotential.h"
#include <functional>

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

namespace Faunus::Energy {

EwaldData::EwaldData(const json &j) {
    alpha = j.at("alpha");                          // damping-parameter
    r_cutoff = j.at("cutoff");                      // real space cut-off
    n_cutoff = j.at("ncutoff");                     // reciprocal space cut-off
    use_spherical_sum = j.value("spherical_sum", true); // Using spherical summation of k-vectors in reciprocal space?
    bjerrum_length = pc::bjerrumLength(j.at("epsr"));
    surface_dielectric_constant = j.value("epss", 0.0); // dielectric constant of surrounding medium
    const_inf =
        (surface_dielectric_constant < 1) ? 0 : 1; // if unphysical (<1) use epsr infinity for surrounding medium
    kappa = j.value("kappa", 0.0);
    kappa_squared = kappa * kappa;

    if (j.count("kcutoff")) {
        faunus_logger->warn("`kcutoff` is deprecated, use `ncutoff` instead");
        n_cutoff = j.at("kcutoff");
    } else {
        n_cutoff = j.at("ncutoff");
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

void to_json(json &j, const EwaldData &d) {
    j = {{"lB", d.bjerrum_length},
         {"epss", d.surface_dielectric_constant},
         {"alpha", d.alpha},
         {"cutoff", d.r_cutoff},
         {"ncutoff", d.n_cutoff},
         {"wavefunctions", d.k_vectors.cols()},
         {"spherical_sum", d.use_spherical_sum},
         {"kappa", d.kappa},
         {"ewaldscheme", d.policy}};
}

TEST_CASE("[Faunus] Ewald - EwaldData") {
    using doctest::Approx;

    Space spc;
    EwaldData data(R"({
                "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "ncutoff": 11.0, "spherical_sum": true, "cutoff": 5.0})"_json);

    CHECK(data.policy == EwaldData::PBC);
    CHECK(data.const_inf == 1);
    CHECK(data.alpha == 0.894427190999916);

    // Check number of wave-vectors using PBC
    PolicyIonIon ionion;
    ionion.updateBox(data, Point(10, 10, 10));
    CHECK(data.k_vectors.cols() == 2975);
    CHECK(data.Q_ion.size() == data.k_vectors.cols());

    // Check number of wave-vectors using IPBC
    data.policy = EwaldData::IPBC;
    PolicyIonIonIPBC ionionIPBC;
    ionionIPBC.updateBox(data, Point(10, 10, 10));
    CHECK(data.k_vectors.cols() == 846);
    CHECK(data.Q_ion.size() == data.k_vectors.cols());
}

//----------------- Ewald Policies -------------------

std::shared_ptr<EwaldPolicyBase> EwaldPolicyBase::makePolicy(EwaldData::Policies policy) {
    switch (policy) {
    case EwaldData::PBC:
        return std::make_shared<PolicyIonIon>();
    case EwaldData::PBCEigen:
        return std::make_shared<PolicyIonIonEigen>();
    case EwaldData::IPBC:
        return std::make_shared<PolicyIonIonIPBC>();
    case EwaldData::IPBCEigen:
        return std::make_shared<PolicyIonIonIPBCEigen>();
    case EwaldData::INVALID:
        throw std::runtime_error("invalid Ewald policy");
    }
    return nullptr;
}

PolicyIonIon::PolicyIonIon() { cite = "doi:10.1063/1.481216"; }
PolicyIonIonIPBC::PolicyIonIonIPBC() { cite = "doi:10/css8"; }

/**
 * Resize k-vectors according to current variables and box length
 */
void PolicyIonIon::updateBox(EwaldData &d, const Point &box) const {
    assert(d.policy == EwaldData::PBC or d.policy == EwaldData::PBCEigen);
    d.box_length = box;
    int n_cutoff_ceil = ceil(d.n_cutoff);
    d.check_k2_zero = 0.1 * std::pow(2 * pc::pi / d.box_length.maxCoeff(), 2);
    int k_vector_size = (2 * n_cutoff_ceil + 1) * (2 * n_cutoff_ceil + 1) * (2 * n_cutoff_ceil + 1) - 1;
    if (k_vector_size == 0) {
        d.k_vectors.resize(3, 1);
        d.Aks.resize(1);
        d.k_vectors.col(0) = Point(1, 0, 0); // Just so it is not the zero-vector
        d.Aks[0] = 0;
        d.num_kvectors = 1;
        d.Q_ion.resize(1);
        d.Q_dipole.resize(1);
    } else {
        double nc2 = d.n_cutoff * d.n_cutoff;
        d.k_vectors.resize(3, k_vector_size);
        d.Aks.resize(k_vector_size);
        d.num_kvectors = 0;
        d.k_vectors.setZero();
        d.Aks.setZero();
        int start_value = 1;
        for (int nx = 0; nx <= n_cutoff_ceil; nx++) {
            double dnx2 = double(nx * nx);
            double factor = (nx > 0) ? 2.0 : 1.0; // optimization of PBC Ewald (and
                                                  // always the case for IPBC Ewald)
            for (int ny = -n_cutoff_ceil * start_value; ny <= n_cutoff_ceil; ny++) {
                double dny2 = double(ny * ny);
                for (int nz = -n_cutoff_ceil * start_value; nz <= n_cutoff_ceil; nz++) {
                    Point kv = 2 * pc::pi * Point(nx, ny, nz).cwiseQuotient(d.box_length);
                    double k2 = kv.squaredNorm() + d.kappa_squared; // last term is only for Yukawa-Ewald
                    if (k2 < d.check_k2_zero)                       // Check if k2 != 0
                        continue;
                    if (d.use_spherical_sum) {
                        double dnz2 = double(nz * nz);
                        if ((dnx2 + dny2 + dnz2) / nc2 > 1)
                            continue;
                    }
                    d.k_vectors.col(d.num_kvectors) = kv;
                    d.Aks[d.num_kvectors] = factor * exp(-k2 / (4 * d.alpha * d.alpha)) / k2;
                    d.num_kvectors++;
                }
            }
        }
        d.Q_ion.resize(d.num_kvectors);
        d.Q_dipole.resize(d.num_kvectors);
        d.Aks.conservativeResize(d.num_kvectors);
        d.k_vectors.conservativeResize(3, d.num_kvectors);
    }
}

/**
 * @todo Add OpenMP pragma to first loop
 */
void PolicyIonIon::updateComplex(EwaldData &data, Space::Tgvec &groups) const {
    for (int k = 0; k < data.k_vectors.cols(); k++) {
        const Point &q = data.k_vectors.col(k);
        EwaldData::Tcomplex Q(0, 0);
        for (auto &g : groups) {       // loop over molecules
            for (auto &particle : g) { // loop over active particles
                double qr = q.dot(particle.pos);
                Q += particle.charge * EwaldData::Tcomplex(std::cos(qr),
                                                           std::sin(qr)); // 'Q^q', see eq. 25 in ref.
            }
            data.Q_ion[k] = Q;
        }
    }
}

void PolicyIonIonEigen::updateComplex(EwaldData &data, Space::Tgvec &groups) const {
    auto [pos, charge] = mapGroupsToEigen(groups);                             // throws if inactive particles
    Eigen::MatrixXd kr = pos.matrix() * data.k_vectors;                        // ( N x 3 ) * ( 3 x K ) = N x K
    data.Q_ion.real() = (kr.array().cos().colwise() * charge).colwise().sum(); // real part of 'Q^q', see eq. 25 in ref.
    data.Q_ion.imag() = kr.array().sin().colwise().sum(); // imaginary part of 'Q^q', see eq. 25 in ref.
}

void PolicyIonIon::updateComplex(EwaldData &d, Change &change, Space::Tgvec &groups, Space::Tgvec &oldgroups) const {
    assert(groups.size() == oldgroups.size());
    for (int k = 0; k < d.k_vectors.cols(); k++) {
        auto &Q = d.Q_ion[k];
        const Point &q = d.k_vectors.col(k);

        for (auto &changed_group : change.groups) {
            auto &g_new = groups.at(changed_group.index);
            auto &g_old = oldgroups.at(changed_group.index);
            for (auto i : changed_group.atoms) {
                if (i < g_new.size()) {
                    double qr = q.dot(g_new[i].pos);
                    Q += g_new[i].charge * EwaldData::Tcomplex(std::cos(qr), std::sin(qr));
                }
                if (i < g_old.size()) {
                    double qr = q.dot(g_old[i].pos);
                    Q -= g_old[i].charge * EwaldData::Tcomplex(std::cos(qr), std::sin(qr));
                }
            }
        }
    }
}

TEST_CASE("[Faunus] Ewald - IonIonPolicy") {
    using doctest::Approx;
    Space spc;
    spc.p.resize(2);
    spc.geo = R"( {"type": "cuboid", "length": 10} )"_json;
    spc.p[0] = R"( {"pos": [0,0,0], "q": 1.0} )"_json;
    spc.p[1] = R"( {"pos": [1,0,0], "q": -1.0} )"_json;
    Group<Particle> g(spc.p.begin(), spc.p.end());
    spc.groups.push_back(g);

    EwaldData data = R"({
                "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "ncutoff": 11.0, "spherical_sum": true, "cutoff": 5.0})"_json;
    Change c;
    c.all = true;
    data.policy = EwaldData::PBC;

    SUBCASE("PBC") {
        PolicyIonIon ionion;
        ionion.updateBox(data, spc.geo.getLength());
        ionion.updateComplex(data, spc.groups);
        CHECK(ionion.selfEnergy(data, c, spc.groups) == Approx(-1.0092530088080642 * data.bjerrum_length));
        CHECK(ionion.surfaceEnergy(data, c, spc.groups) == Approx(0.0020943951023931952 * data.bjerrum_length));
        CHECK(ionion.reciprocalEnergy(data) == Approx(0.21303063979675319 * data.bjerrum_length));
    }

    SUBCASE("PBCEigen") {
        PolicyIonIonEigen ionion;
        ionion.updateBox(data, spc.geo.getLength());
        ionion.updateComplex(data, spc.groups);
        CHECK(ionion.selfEnergy(data, c, spc.groups) == Approx(-1.0092530088080642 * data.bjerrum_length));
        CHECK(ionion.surfaceEnergy(data, c, spc.groups) == Approx(0.0020943951023931952 * data.bjerrum_length));
        CHECK(ionion.reciprocalEnergy(data) == Approx(0.21303063979675319 * data.bjerrum_length));
    }

    SUBCASE("IPBC") {
        PolicyIonIonIPBC ionion;
        data.policy = EwaldData::IPBC;
        ionion.updateBox(data, spc.geo.getLength());
        ionion.updateComplex(data, spc.groups);
        CHECK(ionion.selfEnergy(data, c, spc.groups) == Approx(-1.0092530088080642 * data.bjerrum_length));
        CHECK(ionion.surfaceEnergy(data, c, spc.groups) == Approx(0.0020943951023931952 * data.bjerrum_length));
        CHECK(ionion.reciprocalEnergy(data) == Approx(0.0865107467 * data.bjerrum_length));
    }

    // IPBCEigen is under construction
    /*SUBCASE("IPBCEigen") {
        PolicyIonIonIPBCEigen ionion();
        data.type = EwaldData::IPBCEigen;
        ionion.updateBox(data, spc.geo.getLength());
        ionion.updateComplex(data, spc.groups);
        CHECK(ionion.selfEnergy(data, c, spc.groups) == Approx(-1.0092530088080642 * data.lB));
        CHECK(ionion.surfaceEnergy(data, c, spc.groups) == Approx(0.0020943951023931952 * data.lB));
        CHECK(ionion.reciprocalEnergy(data) == Approx(0.0865107467 * data.lB));
    }*/
}

TEST_CASE("[Faunus] Ewald - IonIonPolicy Benchmarks") {
    Space spc;
    spc.geo = R"( {"type": "cuboid", "length": 80} )"_json;
    spc.p.resize(200);
    for (auto &p : spc.p) {
        p.charge = 1.0;
        p.pos = (random() - 0.5) * spc.geo.getLength();
    }
    Group<Particle> g(spc.p.begin(), spc.p.end());
    spc.groups.push_back(g);

    EwaldData data(R"({
                "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "ncutoff": 11.0, "spherical_sum": true, "cutoff": 9.0})"_json);
    Change c;
    c.all = true;
    data.policy = EwaldData::PBC;

    {
        PolicyIonIon pbc;
        PolicyIonIonEigen pbc_eigen;
        pbc.updateBox(data, spc.geo.getLength());
        pbc_eigen.updateBox(data, spc.geo.getLength());

        ankerl::nanobench::Config bench;
        bench.minEpochIterations(20);
        bench.run("PBC", [&] { pbc.updateComplex(data, spc.groups); }).doNotOptimizeAway();
        bench.run("PBCEigen", [&] { pbc_eigen.updateComplex(data, spc.groups); }).doNotOptimizeAway();
    }
}

//----------------- IPBC Ewald -------------------

/**
 * Resize k-vectors according to current variables and box length
 */
void PolicyIonIonIPBC::updateBox(EwaldData &data, const Point &box) const {
    assert(data.policy == EwaldData::IPBC or data.policy == EwaldData::IPBCEigen);
    data.box_length = box;
    int ncc = std::ceil(data.n_cutoff);
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
        double nc2 = data.n_cutoff * data.n_cutoff;
        data.k_vectors.resize(3, k_vector_size);
        data.Aks.resize(k_vector_size);
        data.num_kvectors = 0;
        data.k_vectors.setZero();
        data.Aks.setZero();
        int start_value = 0;
        for (int nx = 0; nx <= ncc; nx++) {
            double dnx2 = double(nx * nx);
            double xfactor = (nx > 0) ? 2.0 : 1.0; // optimization of PBC Ewald
            for (int ny = -ncc * start_value; ny <= ncc; ny++) {
                double dny2 = double(ny * ny);
                double yfactor = (ny > 0) ? 2.0 : 1.0; // optimization of PBC Ewald
                for (int nz = -ncc * start_value; nz <= ncc; nz++) {
                    double factor = xfactor * yfactor;
                    if (nz > 0)
                        factor *= 2;
                    Point kv = 2 * pc::pi * Point(nx, ny, nz).cwiseQuotient(data.box_length);
                    double k2 = kv.squaredNorm() + data.kappa_squared; // last term is only for Yukawa-Ewald
                    if (k2 < data.check_k2_zero)                       // Check if k2 != 0
                        continue;
                    if (data.use_spherical_sum) {
                        double dnz2 = double(nz * nz);
                        if ((dnx2 + dny2 + dnz2) / nc2 > 1)
                            continue;
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

void PolicyIonIonIPBC::updateComplex(EwaldData &d, Space::Tgvec &groups) const {
    assert(d.policy == EwaldData::IPBC or d.policy == EwaldData::IPBCEigen);
    for (int k = 0; k < d.k_vectors.cols(); k++) {
        const Point &q = d.k_vectors.col(k);
        EwaldData::Tcomplex Q(0, 0);
        for (auto &g : groups) {
            for (auto &particle : g) {
                Q += q.cwiseProduct(particle.pos).array().cos().prod() * particle.charge; // see eq. 2 in doi:10/css8
            }
        }
        d.Q_ion[k] = Q;
    }
}

void PolicyIonIonIPBCEigen::updateComplex(EwaldData &d, Space::Tgvec &groups) const {
    assert(d.policy == EwaldData::IPBC or d.policy == EwaldData::IPBCEigen);
    auto [pos, charge] = mapGroupsToEigen(groups); // throws if inactive particles
    d.Q_ion.real() = (d.k_vectors.array().cwiseProduct(pos).array().cos().prod() * charge)
                         .colwise()
                         .sum(); // see eq. 2 in doi:10/css8
}

void PolicyIonIonIPBC::updateComplex(EwaldData &d, Change &change, Space::Tgvec &groups,
                                     Space::Tgvec &oldgroups) const {
    assert(d.policy == EwaldData::IPBC or d.policy == EwaldData::IPBCEigen);
    assert(groups.size() == oldgroups.size());

    for (int k = 0; k < d.k_vectors.cols(); k++) {
        auto &Q = d.Q_ion[k];
        const Point &q = d.k_vectors.col(k);
        for (auto &changed_group : change.groups) {
            auto &g_new = groups.at(changed_group.index);
            auto &g_old = oldgroups.at(changed_group.index);
            for (auto i : changed_group.atoms) {
                if (i < g_new.size())
                    Q += q.cwiseProduct(g_new[i].pos).array().cos().prod() * g_new[i].charge;
                if (i < g_old.size())
                    Q -= q.cwiseProduct(g_old[i].pos).array().cos().prod() * g_old[i].charge;
            }
        }
    }
}

double PolicyIonIon::surfaceEnergy(const EwaldData &d, Change &change, Space::Tgvec &groups) {
    if (d.const_inf < 0.5)
        return 0;
    Point qr(0, 0, 0);
    if (change.all or change.dV) {
        for (auto &g : groups) {
            for (auto &particle : g) {
                qr += particle.charge * particle.pos;
            }
        }
    } else if (change.groups.size() > 0) {
        for (auto &changed_group : change.groups) {
            auto &g = groups.at(changed_group.index);
            for (auto i : changed_group.atoms) {
                if (i < g.size()) {
                    qr += g[i].charge * g[i].pos;
                }
            }
        }
    }
    double volume = d.box_length.prod();
    return d.const_inf * 2 * pc::pi / ((2 * d.surface_dielectric_constant + 1) * volume) * qr.dot(qr) *
           d.bjerrum_length;
}

double PolicyIonIon::selfEnergy(const EwaldData &d, Change &change, Space::Tgvec &groups) {
    double charges_squared = 0;
    double charge_total = 0;
    if (change.dN) {
        for (auto &changed_group : change.groups) {
            auto &g = groups.at(changed_group.index);
            for (auto i : changed_group.atoms) {
                if (i < g.size()) {
                    charges_squared += std::pow(g[i].charge, 2);
                    charge_total += g[i].charge;
                }
            }
        }
    } else if (change.all and not change.dV) {
        for (auto &g : groups) {
            for (auto &particle : g) {
                charges_squared += particle.charge * particle.charge;
                charge_total += particle.charge;
            }
        }
    }
    double Vcc = -pc::pi / 2.0 / d.alpha / d.alpha / ( d.box_length[0] * d.box_length[1] * d.box_length[2] ) * charge_total * charge_total; // compensate with neutralizing background (if non-zero total charge in system)
    double beta = d.kappa / (2.0 * d.alpha);
    if( beta > 1e-6 )
        Vcc *= ( 1.0 - exp( -beta * beta ) ) / beta / beta; // same as above but for Yukawa-systems
    return ( -d.alpha * charges_squared / std::sqrt(pc::pi) * (std::exp(-beta*beta) + std::sqrt(pc::pi) * beta * std::erf(beta)) + Vcc) * d.bjerrum_length;
}

/**
 * Updates the reciprocal space terms 'Q^q' and 'A_k'.
 * See eqs. 24 and 25 in ref. for PBC Ewald, and eq. 2 in doi:10/css8 for IPBC Ewald.
 */
double PolicyIonIon::reciprocalEnergy(const EwaldData &d) {
    double energy = 0;
    for (int k = 0; k < d.Q_ion.size(); k++) {
        energy += d.Aks[k] * std::norm(d.Q_ion[k]);
    }
    return 2 * pc::pi * energy * d.bjerrum_length / d.box_length.prod();
}

double PolicyIonIonEigen::reciprocalEnergy(const EwaldData &d) {
    double energy = d.Aks.cwiseProduct(d.Q_ion.cwiseAbs2()).sum();
    return 2 * pc::pi * d.bjerrum_length * energy / d.box_length.prod();
}

Ewald::Ewald(const json &j, Space &spc) : data(j), spc(spc) {
    name = "ewald";
    policy = EwaldPolicyBase::makePolicy(data.policy);
    citation_information = policy->cite;
    init();
}

void Ewald::init() {
    policy->updateBox(data, spc.geo.getLength());
    policy->updateComplex(data, spc.groups); // brute force. todo: be selective
}

double Ewald::energy(Change &change) {
    double u = 0;
    if (change) {
        // If the state is NEW_MONTE_CARLO_STATE (trial state), then update all k-vectors
        if (key == TRIAL_MONTE_CARLO_STATE) {
            if (change.all or change.dV) { // everything changes
                policy->updateBox(data, spc.geo.getLength());
                policy->updateComplex(data, spc.groups); // update all (expensive!)
            } else { // much cheaper partial update
              if (change.groups.size() > 0) {
                assert(old_groups != nullptr);
                policy->updateComplex(data, change, spc.groups, *old_groups);
              }
            }
        }
        // the selfEnergy() is omitted as this is added as a separate term in `Hamiltonian`
        // (The pair-potential is responsible for this)
        u = policy->surfaceEnergy(data, change, spc.groups) + policy->reciprocalEnergy(data);
    }
    return u;
}

/**
 * @param forces Destination force vector
 *
 * Calculate forces from reciprocal space. Note that
 * the destination force vector will *not* be zeroed
 * before addition.
 */
void Ewald::force(std::vector<Point> &forces) {
    assert(forces.size() == spc.p.size());
    const double volume = spc.geo.getVolume();

    // Surface contribution
    Point total_dipole_moment = {0.0, 0.0, 0.0};
    for (auto &particle : spc.p) {
        auto mu = particle.hasExtension() ? particle.getExt().mu * particle.getExt().mulen : Point(0, 0, 0);
        total_dipole_moment += particle.pos * particle.charge + mu;
    }

    auto force = forces.begin(); // iterator to force vector on first particle

    assert(data.k_vectors.cols() == data.Q_ion.size());
    data.Q_dipole.resize(data.Q_ion.size());

    for (auto &particle : spc.p) { // loop over particles
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
 * @todo Implement a sync() function in EwaldData to selectively copy information
 */
void Ewald::sync(Energybase* energybase, const Change& change) {
    if (auto* other = dynamic_cast<Ewald*>(energybase)) {
        if (other->key == ACCEPTED_MONTE_CARLO_STATE) {
            old_groups = &(
                other->spc
                    .groups); // give NEW_MONTE_CARLO_STATE access to OLD_MONTE_CARLO_STATE space for optimized updates
        }

        // hard-coded sync; should be expanded when dipolar ewald is supported
        if (change.all or change.dV) {
            other->data.Q_dipole.resize(0); // dipoles are currently unsupported
            data = other->data;
        } else {
            data.Q_ion = other->data.Q_ion;
        }
    } else {
        throw std::runtime_error("sync error");
    }
}

void Ewald::to_json(json &j) const { j = data; }

double Example2D::energy(Change &) {
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

Example2D::Example2D(const json &j, Space &spc) : particle(spc.p.at(0).pos) {
    scale_energy = j.value("scale", 1.0);
    use_2d = j.value("2D", true);
    name = "Example2D";
}
void Example2D::to_json(json &j) const {
    j["scale"] = scale_energy;
    j["2D"] = use_2d;
}

double ContainerOverlap::energy(Change& change) {
    if (change && spc.geo.type != Geometry::CUBOID) { // no need to check in PBC systems
        // *all* groups
        if (change.dV or change.all) {
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
                               [&](const auto& position) { return spc.geo.collision(position); });
    return outside ? pc::infty : 0.0;
}

/**
 * @brief Check a single group based on change object
 * @return true is any particle in group is outside; false otherwise
 */
bool ContainerOverlap::groupIsOutsideContainer(const Change::data& group_change) const {
    const auto& group = spc.groups.at(group_change.index);
    // *all* atoms
    if (group_change.all) {
        return std::any_of(group.begin(), group.end(),
                           [&](auto const& particle) { return spc.geo.collision(particle.pos); });
    }
    // *subset* of atoms
    for (const auto particle_index : group_change.atoms) {
        if (particle_index < group.size()) { // condition due to speciation move?
            if (spc.geo.collision((group.begin() + particle_index)->pos)) {
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
    citation_information = "Frenkel & Smith 2nd Ed (Eq. 5.4.13)";

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
    _roundjson(j, 5);
}

/**
 * @brief Calculates the energy contribution p × V / kT - (N + 1) × ln(V)
 * @return energy in kT
 */
double Isobaric::energy(Change& change) {
    if (change.dV || change.all || change.dN) {
        auto group_is_active = [](const auto& group) { return !group.empty(); };
        auto count_particles = [](const auto& group) { return group.isAtomic() ? group.size() : 1; };
        auto particles_per_group = spc.groups | ranges::cpp20::views::filter(group_is_active) |
                                   ranges::cpp20::views::transform(count_particles);
        auto number_of_particles = std::accumulate(particles_per_group.begin(), particles_per_group.end(), 0);
        const auto volume = spc.geo.getVolume();
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
        CHECK(j.at("P/atm").get<double>() == doctest::Approx(0.5));
        Isobaric(json({{"P/bar", 0.4}}), spc).to_json(j);
        CHECK(j.at("P/bar").get<double>() == doctest::Approx(0.4));
        Isobaric(json({{"P/kT", 0.3}}), spc).to_json(j);
        CHECK(j.at("P/kT").get<double>() == doctest::Approx(0.3));
        Isobaric(json({{"P/mM", 0.2}}), spc).to_json(j);
        CHECK(j.at("P/mM").get<double>() == doctest::Approx(0.2));
        Isobaric(json({{"P/Pa", 0.1}}), spc).to_json(j);
        CHECK(j.at("P/Pa").get<double>() == doctest::Approx(0.1));
    }
}

Constrain::Constrain(const json &j, Space &spc) {
    using namespace Faunus::ReactionCoordinate;
    name = "constrain";
    type = j.at("type").get<std::string>();
    rc = ReactionCoordinate::createReactionCoordinate({{type, j}}, spc);
}

double Constrain::energy(Change &change) {
    if (change) {
        double val = (*rc)();     // calculate reaction coordinate
        if (not rc->inRange(val)) // is it within allowed range?
            return pc::infty;     // if not, return infinite energy
    }
    return 0;
}
void Constrain::to_json(json &j) const {
    j = json(*rc).at(type);
    j.erase("resolution");
    j["type"] = type;
}

void Bonded::updateGroupBonds(const Space::Tgroup& group) {
    const auto first_particle_index = spc.getFirstParticleIndex(group);
    const auto group_index = spc.getGroupIndex(group);
    auto& bonds = internal_bonds[group_index];              // access or insert
    for (const auto& generic_bond : group.traits().bonds) { // generic bonds defined in topology
        const auto& bond = bonds.template push_back<Potential::BondData>(generic_bond->clone());
        bond->shiftIndices(first_particle_index); // shift to absolute particle index
        bond->setEnergyFunction(spc.p);
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
    auto bond_energy = [&](const auto& bond) { return bond->energyFunc(spc.geo.getDistanceFunc()); };
    return std::transform_reduce(bonds.begin(), bonds.end(), 0.0, std::plus<>(), bond_energy);
#else
    double energy = 0.0;
    for (const auto& bond : bonds) {
        energy += bond->energyFunc(spc.geo.getDistanceFunc());
    }
    return energy;
#endif
}

Bonded::Bonded(const Space& spc, const BondVector& external_bonds = BondVector())
    : spc(spc), external_bonds(external_bonds) {
    name = "bonded";
    updateInternalBonds();
    for (auto& bond : this->external_bonds) {
        bond->setEnergyFunction(spc.p);
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

double Bonded::energy(Change &change) {
    double energy = 0.0;
    if (change) {
        energy += sumBondEnergy(external_bonds);
        if (change.all || change.dV) { // calc. for everything!
            for (const auto& [group_index, bonds] : internal_bonds) {
                if (!spc.groups[group_index].empty()) {
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

double Bonded::internalGroupEnergy(const Change::data& changed) {
    using namespace ranges::cpp20::views; // @todo cpp20 --> std::ranges
    double energy = 0.0;
    const auto& group = spc.groups.at(changed.index);
    if (changed.internal && !group.empty()) {
        const auto& bonds = internal_bonds.at(changed.index);
        if (changed.all) { // all internal positions updated
            energy += sumBondEnergy(bonds);
        } else { // only partial update of affected atoms
            const auto first_particle_index = spc.getFirstParticleIndex(group);
            auto particle_indices =
                changed.atoms | transform([first_particle_index](auto i) { return i + first_particle_index; });
            energy += sum_energy(bonds, particle_indices);
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
    auto distance_function = spc.geo.getDistanceFunc();

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
            emplace_back<Energy::Ewald>(_j, spc);
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
    if (spc.geo.type != Geometry::CUBOID) {
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

double Hamiltonian::energy(Change& change) {
    latest_energies.clear();
    for (auto& energy_ptr : energy_terms) {
        energy_ptr->key = key; // is this needed?
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

void Hamiltonian::sync(Energybase* other_hamiltonian, const Change& change) {
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
std::shared_ptr<Energybase> Hamiltonian::createEnergy(Space& spc, const std::string& name, const json& j) {
    using namespace Potential;
    using CoulombLJ = CombinedPairPotential<NewCoulombGalore, LennardJones>;
    using CoulombWCA = CombinedPairPotential<NewCoulombGalore, WeeksChandlerAndersen>;
    using PrimitiveModelWCA = CombinedPairPotential<Coulomb, WeeksChandlerAndersen>;
    using PrimitiveModel = CombinedPairPotential<Coulomb, HardSphere>;

    // only a single pairing policy and cutoff scheme so far
    using PairingPolicy = GroupPairing<GroupPairingPolicy<GroupCutoff>>;

    try {
        if (name == "nonbonded_coulomblj" || name == "nonbonded_newcoulomblj") {
            return std::make_shared<Nonbonded<PairEnergy<CoulombLJ, false>, PairingPolicy>>(j, spc, *this);
        }
        if (name == "nonbonded_coulomblj_EM") {
            return std::make_shared<NonbondedCached<PairEnergy<CoulombLJ, false>, PairingPolicy>>(j, spc, *this);
        }
        if (name == "nonbonded_splined") {
            return std::make_shared<Nonbonded<PairEnergy<Potential::SplinedPotential, false>, PairingPolicy>>(j, spc,
                                                                                                              *this);
        }
        if (name == "nonbonded" || name == "nonbonded_exact") {
            return std::make_shared<Nonbonded<PairEnergy<Potential::FunctorPotential, true>, PairingPolicy>>(j, spc,
                                                                                                             *this);
        }
        if (name == "nonbonded_cached") {
            return std::make_shared<NonbondedCached<PairEnergy<Potential::SplinedPotential>, PairingPolicy>>(j, spc,
                                                                                                             *this);
        }
        if (name == "nonbonded_coulombwca") {
            return std::make_shared<Nonbonded<PairEnergy<CoulombWCA, false>, PairingPolicy>>(j, spc, *this);
        }
        if (name == "nonbonded_pm" || name == "nonbonded_coulombhs") {
            return std::make_shared<Nonbonded<PairEnergy<PrimitiveModel, false>, PairingPolicy>>(j, spc, *this);
        }
        if (name == "nonbonded_pmwca") {
            return std::make_shared<Nonbonded<PairEnergy<PrimitiveModelWCA, false>, PairingPolicy>>(j, spc, *this);
        }
        if (name == "bonded") {
            return std::make_shared<Bonded>(j, spc);
        }
        if (name == "customexternal") {
            return std::make_shared<CustomExternal>(j, spc);
        }
        if (name == "akesson") {
            return std::make_shared<ExternalAkesson>(j, spc);
        }
        if (name == "confine") {
            return std::make_shared<Confine>(j, spc);
        }
        if (name == "constrain") {
            return std::make_shared<Constrain>(j, spc);
        }
        if (name == "example2d") {
            return std::make_shared<Example2D>(j, spc);
        }
        if (name == "isobaric") {
            return std::make_shared<Isobaric>(j, spc);
        }
        if (name == "penalty") {
#ifdef ENABLE_MPI
            return std::make_shared<PenaltyMPI>(j, spc);
#else
            return std::make_shared<Penalty>(j, spc);
#endif
        }
        if (name == "sasa") {
#if defined ENABLE_FREESASA
            return std::make_shared<SASAEnergy>(j, spc);
#else
            throw ConfigurationError("faunus not compiled with sasa support");
#endif
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

SASAEnergy::SASAEnergy(Space &spc, double cosolute_concentration, double probe_radius)
    : spc(spc), cosolute_concentration(cosolute_concentration)
{
    name = "sasa"; // todo predecessor constructor
    citation_information = "doi:10.12688/f1000research.7931.1"; // todo predecessor constructor
    parameters = freesasa_default_parameters;
    parameters.probe_radius = probe_radius;
    init();
}

SASAEnergy::SASAEnergy(const json &j, Space &spc)
    : SASAEnergy(spc, j.value("molarity", 0.0) * 1.0_molar, j.value("radius", 1.4) * 1.0_angstrom) {}

void SASAEnergy::updatePositions([[gnu::unused]] const ParticleVector &p) {
    assert(p.size() == spc.positions().size());
    positions.clear();
    for(auto pos: spc.positions()) {
        auto xyz = pos.data();
        positions.insert(positions.end(), xyz, xyz+3);
    }
}

void SASAEnergy::updateRadii(const ParticleVector &p) {
    radii.resize(p.size());
    std::transform(p.begin(), p.end(), radii.begin(),
                   [](auto &a) { return atoms[a.id].sigma * 0.5; });
}

void SASAEnergy::updateSASA(const ParticleVector &p, const Change &) {
    updateRadii(p);
    updatePositions(p);
    auto result = freesasa_calc_coord(positions.data(), radii.data(), p.size(), &parameters);
    if(result) {
        assert(result->n_atoms == p.size());
        sasa.clear();
        sasa.insert(sasa.begin(), result->sasa, result->sasa + result->n_atoms); // copy
        assert(sasa.size() == p.size());
        freesasa_result_free(result);
    } else {
        throw std::runtime_error("FreeSASA failed");
    }
}

void SASAEnergy::init() {
    auto box = spc.geo.getLength();
    auto box_pbc = box;
    spc.geo.boundary(box_pbc);
    if(box_pbc != box) {
        faunus_logger->error("PBC applied, but PBC not implemented for FreeSASA. Expect unphysical results.");
    }
    Change change;
    change.all = true;
    updateSASA(spc.p, change);
}

double SASAEnergy::energy(Change &change) {
    double u = 0, A = 0;
    updateSASA(spc.p, change); // ideally we want
    for (size_t i = 0; i < spc.p.size(); ++i) {
        auto &a = atoms[spc.p[i].id];
        u += sasa[i] * (a.tension + cosolute_concentration * a.tfe);
        A += sasa[i];
    }
    avgArea += A; // sample average area for accepted confs.
    return u;
}

void SASAEnergy::sync(Energybase* basePtr, const Change& c) {
    auto other = dynamic_cast<decltype(this)>(basePtr);
    if (other) {
        if (c.all || c.dV) {
            radii = other->radii;
            positions = other->positions;
            sasa = other->sasa;
        } else {
            for (auto &d : c.groups) {
                int offset = std::distance(spc.p.begin(), spc.groups.at(d.index).begin());
                for (int j : d.atoms) {
                    int i = j + offset;
                    radii[i] = other->radii[i];
                    sasa[i] = other->sasa[i];
                    for(size_t k = 0; k < 3; ++k) {
                        positions[3*i + k] = spc.positions()[i][k];
                    }
                }
            }
        }
    }
}

void SASAEnergy::to_json(json &j) const {
    using namespace u8;
    j["molarity"] = cosolute_concentration / 1.0_molar;
    j["radius"] = parameters.probe_radius / 1.0_angstrom;
    j[bracket("SASA") + "/" + angstrom + squared] = avgArea.avg() / 1.0_angstrom;
    _roundjson(j, 5); // set json output precision
}

TEST_CASE("[Faunus] FreeSASA") {
    using doctest::Approx;
    Change change; // change object telling that a full energy calculation
    change.all = true;
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
    spc.p[0].pos = {0.0, 0.0, 0.0};
    spc.p[1].pos = {0.0, 0.0, 20.0};

    SUBCASE("Separated atoms") {
        SASAEnergy sasa(spc, 1.5_molar, 1.4_angstrom);
        CHECK(sasa.energy(change) == Approx(4 * pc::pi * (3.4 * 3.4 + 2.6 * 2.6) * 1.5 * 1.0_kJmol));
    }

    SUBCASE("Intersecting atoms") {
        SASAEnergy sasa(spc, 1.5_molar, 1.4_angstrom);
        std::vector<double> distance = {0.0, 2.5, 5.0, 7.5, 10.0};
        std::vector<double> sasa_energy = {87.3576, 100.4612, 127.3487, 138.4422, 138.4422};
        for (size_t i = 0; i < distance.size(); ++i) {
            spc.p[1].pos = {0.0, 0.0, distance[i]};
            CHECK(sasa.energy(change) == Approx(sasa_energy[i]).epsilon(0.02));
        }
    }

    SUBCASE("PBC") {}
}
#endif

//==================== GroupCutoff ====================

GroupCutoff::GroupCutoff(Space::Tgeometry &geometry) : geometry(geometry) {}

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
                    if (const auto molecules_names = words2vec<std::string>(named_pair); molecules_names.size() == 2) {
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
        CHECK(cutoff.getCutoff(0, 0) == doctest::Approx(12.0));
        CHECK(cutoff.getCutoff(0, 1) == doctest::Approx(12.0));
        CHECK(cutoff.getCutoff(1, 1) == doctest::Approx(12.0));
    }

    SUBCASE("new style default") {
        GroupCutoff cutoff(geo);
        from_json(R"({"cutoff_g2g": {"default": 13.0}})"_json, cutoff);
        CHECK(cutoff.getCutoff(0, 0) == doctest::Approx(13.0));
        CHECK(cutoff.getCutoff(0, 1) == doctest::Approx(13.0));
        CHECK(cutoff.getCutoff(1, 1) == doctest::Approx(13.0));
    }

    SUBCASE("custom") {
        GroupCutoff cutoff(geo);
        from_json(R"({"cutoff_g2g": {"default": 13.0, "M Q": 14.0}})"_json, cutoff);
        CHECK(cutoff.getCutoff(0, 0) == doctest::Approx(13.0));
        CHECK(cutoff.getCutoff(0, 1) == doctest::Approx(14.0));
        CHECK(cutoff.getCutoff(1, 1) == doctest::Approx(13.0));
    }

    SUBCASE("custom - no default value") {
        GroupCutoff cutoff(geo);
        from_json(R"({"cutoff_g2g": {"M Q": 11.0}})"_json, cutoff);
        CHECK(cutoff.getCutoff(0, 0) == doctest::Approx(std::sqrt(pc::max_value)));
        CHECK(cutoff.getCutoff(0, 1) == doctest::Approx(11.0));
        CHECK(cutoff.getCutoff(1, 1) == doctest::Approx(std::sqrt(pc::max_value)));
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
        scheme = SERIAL;
    }
#endif
#ifndef _OPENMP
    if (scheme == Scheme::OPENMP) {
        faunus_logger->warn("'openmp' summation unavailable; falling back to 'serial'");
        scheme = SERIAL;
    }
#endif
    faunus_logger->debug("setting summation policy to {}", json(scheme).dump(1));
}

void EnergyAccumulatorBase::to_json(json& j) const { j["summation_policy"] = scheme; }

} // end of namespace Faunus::Energy
