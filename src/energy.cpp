#include "energy.h"
#include "penalty.h"
#include "potentials.h"
#include "externalpotential.h"

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
void Ewald::sync(Energybase *energybase_pointer, Change &change) {
    auto other = dynamic_cast<decltype(this)>(energybase_pointer);
    assert(other);
    if (other->key == ACCEPTED_MONTE_CARLO_STATE) {
        old_groups =
            &(other->spc
                  .groups); // give NEW_MONTE_CARLO_STATE access to OLD_MONTE_CARLO_STATE space for optimized updates
    }

    // hard-coded sync; should be expanded when dipolar ewald is supported
    if (change.all or change.dV) {
        other->data.Q_dipole.resize(0); // dipoles are currently unsupported
        data = other->data;
    } else {
        data.Q_ion = other->data.Q_ion;
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

double ContainerOverlap::energy(Change &change) {
    // if (spc.geo.type not_eq Geometry::CUBOID) // cuboid have PBC in all directions
    if (change) {
        // all groups have been updated
        if (change.dV or change.all) {
            for (auto &g : spc.groups) // loop over *all* groups in system
                for (auto &p : g)      // loop over *all* active particles in group
                    if (spc.geo.collision(p.pos))
                        return pc::infty;
            return 0;
        }

        // only a subset of groups have been updated
        for (auto &d : change.groups) {
            auto &g = spc.groups[d.index];
            // all atoms were updated
            if (d.all) {
                for (auto &p : g) // loop over *all* active particles in group
                    if (spc.geo.collision(p.pos))
                        return pc::infty;
            } else
                // only a subset of atoms were updated
                for (int i : d.atoms) // loop over specific atoms
                    if (i < g.size())
                        if (spc.geo.collision((g.begin() + i)->pos))
                            return pc::infty;
        }
    }
    return 0;
}

// ------------- Isobaric ---------------

Isobaric::Isobaric(const json &j, Space &spc) : spc(spc) {
    name = "isobaric";
    citation_information = "Frenkel & Smith 2nd Ed (Eq. 5.4.13)";
    P = j.value("P/mM", 0.0) * 1.0_millimolar;
    if (P < 1e-10) {
        P = j.value("P/Pa", 0.0) * 1.0_Pa;
        if (P < 1e-10)
            P = j.at("P/atm").get<double>() * 1.0_atm;
    }
}

double Isobaric::energy(Change &change) {
    if (change.dV || change.all || change.dN) {
        size_t N = 0;
        for (auto &g : spc.groups) {
            if (!g.empty()) {
                if (g.atomic)
                    N += g.size();
                else
                    N++;
            }
        }
        double V = spc.geo.getVolume();
        return P * V - (N+1)*std::log(V);
    } else
        return 0;
}
void Isobaric::to_json(json &j) const {
    j["P/atm"] = P / 1.0_atm;
    j["P/mM"] = P / 1.0_millimolar;
    j["P/Pa"] = P / 1.0_Pa;
    _roundjson(j, 5);
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
void Bonded::update_intra() {
    using namespace Potential;
    intra.clear();
    for (size_t i = 0; i < spc.groups.size(); i++) {
        const auto &group = spc.groups[i];
        for (const auto &bond : Faunus::molecules.at(group.id).bonds) {
            intra[i].push_back<BondData>(bond->clone()); // deep copy BondData from MoleculeData
            intra[i].back()->shift(std::distance(spc.p.begin(), group.begin()));
            Potential::setBondEnergyFunction(intra[i].back(), spc.p);
        }
    }
}
double Bonded::sum_energy(const Bonded::BondVector &bonds) const {
    double energy = 0;
    for (const auto &bond : bonds) {
        assert(bond->hasEnergyFunction());
        energy += bond->energyFunc(spc.geo.getDistanceFunc());
    }
    return energy;
}

Bonded::Bonded(const json &j, Space &spc) : spc(spc) {
    name = "bonded";
    update_intra();
    if (j.is_object())
        if (j.count("bondlist") == 1)
            inter = j["bondlist"].get<BondVector>();
    for (auto &i : inter) // set all energy functions
        Potential::setBondEnergyFunction(i, spc.p);
}
void Bonded::to_json(json &j) const {
    if (!inter.empty())
        j["bondlist"] = inter;
    if (!intra.empty()) {
        json &_j = j["bondlist-intramolecular"];
        _j = json::array();
        for (auto &i : intra)
            for (auto &b : i.second)
                _j.push_back(b);
    }
}
double Bonded::energy(Change &change) {
    double energy = 0;
    if (change) {
        energy += sum_energy(inter); // energy of inter-molecular bonds

        if (change.all || change.dV) {              // compute all active groups
            for (const auto &i : intra) {           // energies of intra-molecular bonds
                if (!spc.groups[i.first].empty()) { // add only if group is active
                    energy += sum_energy(i.second);
                }
            }
        } else { // compute only the affected groups
            for (const auto &changed : change.groups) {
                const auto &intra_group = intra[changed.index];
                if (changed.internal) {
                    if (changed.all) { // all internal positions updated
                        if (not spc.groups[changed.index].empty())
                            energy += sum_energy(intra_group);
                    } else { // only partial update of affected atoms
                        // an offset is the index of the first particle in the group
                        const int offset = std::distance(spc.p.begin(), spc.groups[changed.index].begin());
                        // add an offset to the group atom indices to get the absolute indices
                        auto particle_indices =
                            changed.atoms | ranges::cpp20::views::transform([offset](auto i) { return i + offset; });
                        energy += sum_energy(intra_group, particle_indices);
                    }
                }
            }
        } // for-loop over groups
    }
    return energy;
}

/**
 * @param forces Target force vector for *all* particles in the system
 *
 * Each element in `force` represent the force on a particle and this
 * updates (add) the bonded force.
 *
 * - loop over groups and their internal bonds and _add_ force
 * - loop over inter-molecular bonds and _add_ force
 *
 * Force unit: kT/Å
 *
 * @warning Untested
 */
void Bonded::force(std::vector<Point> &forces) {
    auto distance_function = spc.geo.getDistanceFunc();
    for (const auto &[group_index, bonds] : intra) {                     // loop over all intra-molecular bonds
        const auto &group = spc.groups[group_index];                     // this is the group we're currently working on
        for (const auto &bond : bonds) {                                 // loop over all bonds in group
            assert(bond->forceFunc != nullptr);                          // the force function must be implemented
            const auto bond_forces = bond->forceFunc(distance_function); // get forces on each atom in bond
            assert(bond->index.size() == bond_forces.size());
            for (size_t i = 0; i < bond->index.size(); i++) { // loop over atom index in bond (relative to group begin)
                const auto absolute_index = std::distance(spc.p.begin(), group.begin()) + bond->index[i];
                assert(absolute_index < forces.size());
                forces[absolute_index] += bond_forces[i]; // add to overall force
            }
        }
    }

    for (const auto &bond : inter) {                                 // loop over inter-molecular bonds
        assert(bond->forceFunc != nullptr);                          // the force function must be implemented
        const auto bond_forces = bond->forceFunc(distance_function); // get forces on each atom in bond
        assert(bond_forces.size() == bond->index.size());
        for (size_t i = 0; i < bond->index.size(); i++) { // loop over atom index in bond (absolute index)
            forces[bond->index[i]] += bond_forces[i];     // add to overall force
        }
    }
}

//---------- Hamiltonian ------------

void Hamiltonian::to_json(json &j) const {
    for (auto i : this->vec)
        j.push_back(*i);
}

void Hamiltonian::addEwald(const json &j, Space &spc) {
    // note this will currently not detect multipolar energies ore
    // deeply nested "coulomb" pair-potentials
    json _j;
    if (j.count("default") == 1) { // try to detect FunctorPotential
        for (auto &i : j["default"]) {
            if (i.count("coulomb") == 1) {
                _j = i["coulomb"];
                break;
            }
        }
    } else if (j.count("coulomb") == 1)
        _j = j["coulomb"];
    else
        return;

    if (_j.count("type")) {
        if (_j.at("type") == "ewald") {
            faunus_logger->debug("adding Ewald reciprocal and surface energy terms");
            emplace_back<Energy::Ewald>(_j, spc);
        }
    }
}

void Hamiltonian::force(PointVector &forces) {
    for (auto energy_ptr : this->vec) { // loop over terms in Hamiltonian
        energy_ptr->force(forces);      // and update forces
    }
}

Hamiltonian::Hamiltonian(Space& spc, const json& j) {
    using namespace Potential;

    typedef CombinedPairPotential<NewCoulombGalore, LennardJones> CoulombLJ; // temporary name
    typedef CombinedPairPotential<NewCoulombGalore, WeeksChandlerAndersen> CoulombWCA;
    typedef CombinedPairPotential<Coulomb, WeeksChandlerAndersen> PrimitiveModelWCA;
    typedef CombinedPairPotential<Coulomb, HardSphere> PrimitiveModel;

    name = "hamiltonian";

    // add container overlap energy for non-cuboidal geometries
    if (spc.geo.type != Geometry::CUBOID) {
        emplace_back<Energy::ContainerOverlap>(spc);
    }

    // only a single pairing policy and cutoff scheme so far
    typedef GroupPairing<GroupPairingPolicy<GroupCutoff>> PairingPolicy;

    if (!j.is_array()) {
        throw ConfigurationError("energy: json array expected");
    }
    for (auto& j_energy : j) { // loop over energy list
        try {
            const auto& [key, j_params] = jsonSingleItem(j_energy);
            try {
                if (key == "nonbonded_coulomblj" || key == "nonbonded_newcoulomblj") {
                    emplace_back<Nonbonded<PairEnergy<CoulombLJ, false>, PairingPolicy>>(j_params, spc, *this);
                } else if (key == "nonbonded_coulomblj_EM") {
                    emplace_back<NonbondedCached<PairEnergy<CoulombLJ, false>, PairingPolicy>>(j_params, spc, *this);
                } else if (key == "nonbonded_splined") {
                    emplace_back<Nonbonded<PairEnergy<SplinedPotential, false>, PairingPolicy>>(j_params, spc, *this);
                } else if (key == "nonbonded" || key == "nonbonded_exact") {
                    emplace_back<Nonbonded<PairEnergy<FunctorPotential, true>, PairingPolicy>>(j_params, spc, *this);
                } else if (key == "nonbonded_cached") {
                    emplace_back<NonbondedCached<PairEnergy<SplinedPotential>, PairingPolicy>>(j_params, spc, *this);
                } else if (key == "nonbonded_coulombwca") {
                    emplace_back<Nonbonded<PairEnergy<CoulombWCA, false>, PairingPolicy>>(j_params, spc, *this);
                } else if (key == "nonbonded_pm" or key == "nonbonded_coulombhs") {
                    emplace_back<Nonbonded<PairEnergy<PrimitiveModel, false>, PairingPolicy>>(j_params, spc, *this);
                } else if (key == "nonbonded_pmwca") {
                    emplace_back<Nonbonded<PairEnergy<PrimitiveModelWCA, false>, PairingPolicy>>(j_params, spc, *this);
                } else if (key == "bonded") {
                    emplace_back<Bonded>(j_params, spc);
                } else if (key == "customexternal") {
                    emplace_back<CustomExternal>(j_params, spc);
                } else if (key == "akesson") {
                    emplace_back<ExternalAkesson>(j_params, spc);
                } else if (key == "confine") {
                    emplace_back<Confine>(j_params, spc);
                } else if (key == "constrain") {
                    emplace_back<Constrain>(j_params, spc);
                } else if (key == "example2d") {
                    emplace_back<Example2D>(j_params, spc);
                } else if (key == "isobaric") {
                    emplace_back<Isobaric>(j_params, spc);
                } else if (key == "penalty") {
#ifdef ENABLE_MPI
                    emplace_back<PenaltyMPI>(j_params, spc);
#else
                    emplace_back<Penalty>(j_params, spc);
#endif
                } else if (key == "sasa") {
#if defined ENABLE_FREESASA
                    emplace_back<SASAEnergy>(j_params, spc);
#else
                    throw ConfigurationError("not included - recompile Faunus with FreeSASA support");
#endif
                } else if (key == "maxenergy") {
                    // looks like an unfortumate json scheme decision that requires special handling here
                    maxenergy = j_params.get<double>();
                } else {
                    throw ConfigurationError("unknown energy");
                }
                // append additional energies to the if-chain

                // this should be moved into `Nonbonded` and added when appropriate
                // Nonbonded now has access to Hamiltonian (*this) and can therefore
                // add energy terms
                addEwald(j_params, spc); // add reciprocal Ewald terms if appropriate
            } catch (std::exception& e) {
                usageTip.pick(key);
                throw ConfigurationError("'{}': {}", key, e.what());
            }
        } catch (std::exception& e) {
            throw ConfigurationError("energy: {}", e.what()).attachJson(j_energy);
        }
    }

    // Check if there are molecules with bonds and warn if "bonded" has not been added
    for (auto& molecule : Faunus::molecules) {
        if (!molecule.bonds.empty() && find<Energy::Bonded>().empty()) {
            faunus_logger->warn("{} bonds specified in topology but missing in energy", molecule.name);
        }
    }
}

double Hamiltonian::energy(Change &change) {
    double du = 0;
    for (auto i : this->vec) { // loop over terms in Hamiltonian
        i->key = key;
        i->timer.start(); // time each term
        du += i->energy(change);
        i->timer.stop();
        if (du >= maxenergy)
            break; // stop summing energies
    }
    return du;
}
void Hamiltonian::init() {
    for (auto i : this->vec)
        i->init();
}
void Hamiltonian::sync(Energybase *basePtr, Change &change) {
    auto other = dynamic_cast<decltype(this)>(basePtr);
    if (other)
        if (other->size() == size()) {
            for (size_t i = 0; i < size(); i++)
                this->vec[i]->sync(other->vec[i].get(), change);
            return;
        }
    throw std::runtime_error("hamiltonian mismatch");
}

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

void SASAEnergy::sync(Energybase *basePtr, Change &c) {
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

void from_json(const json &j, GroupCutoff &cutoff) {
    // disable all group-to-group cutoffs by setting infinity
    for (auto &i : Faunus::molecules) {
        for (auto &j : Faunus::molecules) {
            cutoff.cutoff_squared.set(i.id(), j.id(), pc::max_value);
        }
    }

    auto it = j.find("cutoff_g2g");
    if (it != j.end()) {
        if (it->is_number()) {
            // old style input w. only a single cutoff
            cutoff.default_cutoff_squared = std::pow(it->get<double>(), 2);
            for (auto &i : Faunus::molecules) {
                for (auto &j : Faunus::molecules) {
                    cutoff.cutoff_squared.set(i.id(), j.id(), cutoff.default_cutoff_squared);
                }
            }
        }
        else if (it->is_object()) {
            // new style input w. multiple cutoffs between molecules
            // ensure that there is a default, fallback cutoff
            cutoff.default_cutoff_squared = std::pow(it->at("default").get<double>(), 2);
            for (auto &i : Faunus::molecules) {
                for (auto &j : Faunus::molecules) {
                    cutoff.cutoff_squared.set(i.id(), j.id(), cutoff.default_cutoff_squared);
                }
            }
            // loop for space separated molecule pairs in keys
            for (auto &i : it->items()) {
                try {
                    if(const auto molecules_names = words2vec<std::string>(i.key()); molecules_names.size() == 2) {
                        const auto molecule1 = findMoleculeByName(molecules_names[0]);
                        const auto molecule2 = findMoleculeByName(molecules_names[1]);
                        cutoff.cutoff_squared.set(molecule1.id(), molecule2.id(), std::pow(i.value().get<double>(), 2));
                    } else {
                        throw std::runtime_error("invalid molecules names");
                    }
                } catch(const std::exception &e) {
                    faunus_logger->warn("Unable to set a custom cutoff “{}” for “{}”.", i.value().dump(), i.key());
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
} // end of namespace Faunus::Energy
