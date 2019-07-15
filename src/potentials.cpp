#include "potentials.h"
#include "multipole.h"

namespace Faunus {
namespace Potential {

// =============== PairMixer ===============

TPairMatrixPtr PairMixer::createPairMatrix(const std::vector<AtomData> &atoms) {
    size_t n = atoms.size(); // number of atom types
    TPairMatrixPtr matrix = std::make_shared<TPairMatrix>(n, n);
    for (auto &i : atoms) {
        for (auto &j : atoms) {
            (*matrix)(i.id(), j.id()) = combinator(extractor(i), extractor(j));
        }
    }
    return matrix;
}

TPairMatrixPtr PairMixer::createPairMatrix(const std::vector<AtomData> &atoms,
        const std::vector<InteractionData> &interactions) {
    TPairMatrixPtr matrix = PairMixer::createPairMatrix(atoms);
    auto dimension = std::min(matrix->rows(), matrix->cols());
    for (auto &i : interactions) {
        if (i.atom_id[0] >= 0 && i.atom_id[1] >= 0 && i.atom_id[0] < dimension && i.atom_id[1] < dimension) {
            // interaction is always symmetric
            // treat it as a homogeneous interaction if some extra transformation is supposed,
            // e.g., to square or multiply by a constant
            (*matrix)(i.atom_id[0], i.atom_id[1]) = (*matrix)(i.atom_id[1], i.atom_id[0]) =
                    combinator(extractor(i.interaction), extractor(i.interaction));
        } else {
            throw std::runtime_error("atomtype index out of range");
        }
    }
    return matrix;
}

void from_json(const json &j, std::vector<InteractionData> &interactions) {
    auto &custom_list = j.at("custom");
    if(! custom_list.is_object() && ! custom_list.is_array()) {
        throw CustomInteractionException("custom parameters syntax error");
    }

    for (auto custom_pair = custom_list.begin(); custom_pair != custom_list.end(); ++custom_pair) {
        // atomdata is an array with items ecapsulated as objects hence we emulate here
        AtomData a = custom_list.is_object() ? json::object({{custom_pair.key(), *custom_pair}}) : (*custom_pair);
        auto atoms_name = words2vec<std::string>(a.name);
        if (atoms_name.size() == 2) {
            auto atom0 = findName(atoms, atoms_name[0]);
            auto atom1 = findName(atoms, atoms_name[1]);
            if (atom0 == atoms.end() or atom1 == atoms.end()) {
                throw CustomInteractionException(
                    ("unknown atom(s): ["s + atoms_name[0] + " " + atoms_name[1] + "]"));
            }
            interactions.push_back({{atom0->id(), atom1->id()}, a});
        } else {
            throw CustomInteractionException("custom parameters require exactly two space-separated atoms");
        }
    }
}

void to_json(json &j, const std::vector<InteractionData> &interactions) {
    if(! interactions.empty()) {
        auto &j_custom = j["custom"];
        for (auto &i : interactions) {
            j_custom = i.interaction;
        }
    }
}


// =============== PairPotentialBase ===============

Point PairPotentialBase::force(const Particle &, const Particle &, double, const Point &) {
    assert(false && "We should never reach this point!");
    return {0, 0, 0};
}

void RepulsionR3::from_json(const json &j) {
    f = j.value("prefactor", 1.0);
    e = j.value("lj-prefactor", 1.0);
    s = j.value("sigma", 1.0);
}

void RepulsionR3::to_json(json &j) const { j = {{"prefactor", f}, {"lj-prefactor", e}, {"sigma", s}}; }

void CosAttract::to_json(json &j) const {
    j = {{"eps", eps / 1.0_kJmol}, {"rc", rc / 1.0_angstrom}, {"wc", wc / 1.0_angstrom}};
}

void CosAttract::from_json(const json &j) {
    eps = j.at("eps").get<double>() * 1.0_kJmol;
    rc = j.at("rc").get<double>() * 1.0_angstrom;
    wc = j.at("wc").get<double>() * 1.0_angstrom;
    rc2 = rc * rc;
    c = pc::pi / 2 / wc;
    rcwc2 = pow((rc + wc), 2);
}

// -------------- CoulombGalore ---------------

void CoulombGalore::sfYukawa(const json &j) {
    kappa = 1.0 / j.at("debyelength").get<double>();
    I = kappa * kappa / (8.0 * lB * pc::pi * pc::Nav / 1e27);
    table = sf.generate([&](double q) { return std::exp(-q * rc * kappa) - std::exp(-kappa * rc); }, 0, 1); // q=r/Rc
    // we could also fill in some info std::string or JSON output...
}

void CoulombGalore::sfYukawaPoisson(const json &j) {
    kappa = 1.0 / j.at("debyelength").get<double>();
    I = kappa * kappa / (8.0 * lB * pc::pi * pc::Nav / 1e27);
    C = j.value("C", 3);
    D = j.value("D", 3);
    if ((C < 1) || (D < 1))
        throw std::runtime_error("`C` and `D` must be larger than zero");
    table = sf.generate(
        [&](double q) {
            double qt = (1.0 - exp(2.0 * kappa * rc * q)) / (1.0 - exp(2.0 * kappa * rc));
            double tmp = 0.0;
            for (int c = 0; c < C; c++)
                tmp += double(factorial(D - 1 + c)) / double(factorial(D - 1)) / double(factorial(c)) * double(C - c) /
                       double(C) * pow(qt, double(c));
            return pow(1.0 - qt, double(D) + 1.0) * tmp;
        },
        0, 1);
    calcDielectric = [&](double M2V) { return 1 + 3 * M2V; };
    selfenergy_prefactor = -double(C + D) / double(C);
}

void CoulombGalore::sfReactionField(const json &j) {
    epsrf = j.at("eps_rf");
    table = sf.generate(
        [&](double q) {
            return 1 + ((epsrf - epsr) / (2 * epsrf + epsr)) * q * q * q - 3 * (epsrf / (2 * epsrf + epsr)) * q;
        },
        0, 1);
    calcDielectric = [&](double M2V) {
        if (epsrf > 1e10)
            return 1 + 3 * M2V;
        if (fabs(epsrf - epsr) < 1e-6)
            return 2.25 * M2V + 0.25 + 0.75 * sqrt(9 * M2V * M2V + 2 * M2V + 1);
        if (fabs(epsrf - 1.0) < 1e-6)
            return (2 * M2V + 1) / (1 - M2V);
        return 0.5 * (2 * epsrf - 1 + sqrt(-72 * M2V * M2V * epsrf + 4 * epsrf * epsrf + 4 * epsrf + 1)) /
               (3 * M2V - 1); // Needs to be checked!
        // return (6*M2V*epsrf + 2*epsrf + 1.0)/(1.0 + 2*epsrf - 3*M2V); // Is OK when epsr=1.0
    };
    selfenergy_prefactor = 1.5 * epsrf / (2.0 * epsrf + epsr); // Correct?!, see Eq.14 in DOI: 10.1021/jp510612w
    // we could also fill in some info std::string or JSON output...
}

void CoulombGalore::sfQpotential(const json &j) {
    order = j.value("order", 300);
    table = sf.generate([&](double q) { return qPochhammerSymbol(q, 1, order); }, 0, 1);
    calcDielectric = [&](double M2V) { return 1 + 3 * M2V; };
    selfenergy_prefactor = -1.0;
}

void CoulombGalore::sfYonezawa(const json &j) {
    alpha = j.at("alpha");
    table = sf.generate([&](double q) { return 1 - std::erfc(alpha * rc) * q + q * q; }, 0, 1);
    calcDielectric = [&](double M2V) { return 1 + 3 * M2V; };
    selfenergy_prefactor = -0.5*(erfc(alpha*rc) + 2.0*alpha*rc / sqrt(pc::pi));
}

void CoulombGalore::sfFanourgakis(const json &) {
    table =
        sf.generate([&](double q) { return 1 - 1.75 * q + 5.25 * pow(q, 5) - 7 * pow(q, 6) + 2.5 * pow(q, 7); }, 0, 1);
    calcDielectric = [&](double M2V) { return 1 + 3 * M2V; };
    selfenergy_prefactor = -0.875;
}

void CoulombGalore::sfPoisson(const json &j) {
    C = j.value("C", 3);
    D = j.value("D", 3);
    if ((C < 1) || (D < 1))
        throw std::runtime_error("`C` and `D` must be larger than zero");
    table = sf.generate(
        [&](double q) {
            double tmp = 0.0;
            for (int c = 0; c < C; c++)
                tmp += double(factorial(D - 1 + c)) / double(factorial(D - 1)) / double(factorial(c)) * double(C - c) /
                       double(C) * pow(q, double(c));
            return pow(1.0 - q, double(D) + 1.0) * tmp;
        },
        0, 1);
    calcDielectric = [&](double M2V) { return 1 + 3 * M2V; };
    selfenergy_prefactor = -double(C + D) / double(C);
}

void CoulombGalore::sfFennel(const json &j) {
    alpha = j.at("alpha");
    table = sf.generate(
        [&](double q) {
            return (
                erfc(alpha * rc * q) - std::erfc(alpha * rc) * q +
                (q - 1.0) * q *
                    (std::erfc(alpha * rc) + 2 * alpha * rc / std::sqrt(pc::pi) * std::exp(-alpha * alpha * rc * rc)));
        },
        0, 1);
    calcDielectric = [&](double M2V) {
        double T = erf(alpha * rc) -
                   (2 / (3 * sqrt(pc::pi))) * exp(-alpha * alpha * rc * rc) *
                       (alpha * alpha * rc * rc * alpha * alpha * rc * rc + 2.0 * alpha * alpha * rc * rc + 3.0);
        return (((T + 2.0) * M2V + 1.0) / ((T - 1.0) * M2V + 1.0));
    };
    selfenergy_prefactor = -(erfc(alpha*rc) + alpha*rc / sqrt(pc::pi) * (1.0 + exp(-alpha*alpha*rc*rc)));
}

void CoulombGalore::sfEwald(const json &j) { // is all this true for kappa \ne 0 ?
    alpha = j.at("alpha");
    kappa = j.value("kappa",0.0);
    table = sf.generate( [&](double q) { return (std::erfc(alpha*rc*q + kappa/2.0/alpha)*std::exp(kappa*rc*q) + std::erfc(alpha*rc*q - kappa/2.0/alpha)*std::exp(-kappa*rc*q)  )/2.0; }, 0, 1 ); // Yukawa potential
    //table = sf.generate( [&](double q) { return std::erfc(alpha*rc*q); }, 0, 1 ); // pure Coulomb potential
    calcDielectric = [&](double M2V) {
        double T = std::erf(alpha * rc) -
                   (2 / (3 * sqrt(pc::pi))) * std::exp(-alpha * alpha * rc * rc) * (2 * alpha * alpha * rc * rc + 3);
        return ((T + 2.0) * M2V + 1) / ((T - 1) * M2V + 1);
    };
    selfenergy_prefactor = alpha * rc / sqrt(pc::pi);
}

void CoulombGalore::sfWolf(const json &j) {
    alpha = j.at("alpha");
    table = sf.generate([&](double q) { return (erfc(alpha * rc * q) - erfc(alpha * rc) * q); }, 0, 1);
    calcDielectric = [&](double M2V) {
        double T = erf(alpha * rc) -
                   (2 / (3 * sqrt(pc::pi))) * exp(-alpha * alpha * rc * rc) * (2.0 * alpha * alpha * rc * rc + 3.0);
        return (((T + 2.0) * M2V + 1.0) / ((T - 1.0) * M2V + 1.0));
    };
    selfenergy_prefactor = -0.5*(erfc(alpha*rc) + 2.0*alpha*rc / sqrt(pc::pi));
}

void CoulombGalore::sfPlain(const json &, double val) {
    table = sf.generate([&](double) { return val; }, 0, 1);
    calcDielectric = [&](double M2V) { return (2.0 * M2V + 1.0) / (1.0 - M2V); };
    selfenergy_prefactor = 0.0;
}

void CoulombGalore::from_json(const json &j) {
    try {
        kappa = 0.0;
        type = j.at("type");
        rc = j.at("cutoff").get<double>();
        rc2 = rc * rc;
        rc1i = 1 / rc;
        epsr = j.at("epsr").get<double>();
        lB = pc::lB(epsr); // Bjerrum length

        depsdt = j.value("depsdt", -0.368 * pc::temperature / epsr);
        sf.setTolerance(j.value("utol", 1e-5), j.value("ftol", 1e-2));

        if (type == "yukawapoisson")
            sfYukawaPoisson(j);
        else if (type == "reactionfield")
            sfReactionField(j);
        else if (type == "fanourgakis")
            sfFanourgakis(j);
        else if (type == "qpotential")
            sfQpotential(j);
        else if (type == "yonezawa")
            sfYonezawa(j);
        else if (type == "poisson")
            sfPoisson(j);
        else if (type == "yukawa")
            sfYukawa(j);
        else if (type == "fennel")
            sfFennel(j);
        else if (type == "plain")
            sfPlain(j, 1);
        else if (type == "ewald")
            sfEwald(j);
        else if (type == "none")
            sfPlain(j, 0);
        else if (type == "wolf")
            sfWolf(j);

        ecs = std::make_shared<PairMatrix<double>>();
        for (auto &i : atoms)
            for (auto &j : atoms) {
                double tmpi = kappa * i.sigma / 2.0;
                double tmpj = kappa * j.sigma / 2.0;
                double ecsi = 1.0;
                double ecsj = 1.0;
                if (tmpi > 1e-6)
                    ecsi = std::sinh(tmpi) / tmpi;
                if (tmpj > 1e-6)
                    ecsj = std::sinh(tmpj) / tmpj;
                ecs->set(i.id(), j.id(), ecsi * ecsj);
            }

        if (table.empty())
            throw std::runtime_error(name + ": unknown coulomb type '" + type + "'");

        // Set particle self-energy function. For reasons yet to be understood,
        // rc etc. cannot be lambda captured by reference, but must be hard-copied,
        // here into `factor`
        selfEnergy = [factor = selfenergy_prefactor * lB / rc](const Particle &a) {
            return a.charge * a.charge * factor;
        };
    }

    catch (std::exception &e) {
        std::cerr << "CoulombGalore error: " << e.what();
        throw;
    }
}

double CoulombGalore::dielectric_constant(double M2V) { return calcDielectric(M2V); }

void CoulombGalore::to_json(json &j) const {
    using namespace u8;
    j["epsr"] = epsr;
    j["T" + partial + epsilon_m + "/" + partial + "T"] = depsdt;
    j["lB"] = lB;
    j["cutoff"] = rc;
    j["type"] = type;
    if (type == "yukawa" || type == "yukawapoisson") {
        j["debyelength"] = 1.0 / kappa;
        j["ionic strength"] = I;
    }
    if (type == "yukawapoisson" || type == "poisson") {
        j["C"] = C;
        j["D"] = D;
    }
    if (type == "qpotential")
        j["order"] = order;
    if (type == "yonezawa" || type == "fennel" || type == "wolf" || type == "ewald")
        j["alpha"] = alpha;
    if (type == "reactionfield") {
        if (epsrf > 1e10)
            j[epsilon_m + "_rf"] = 2e10;
        else
            j[epsilon_m + "_rf"] = epsrf;
    }
    _roundjson(j, 5);
}

//-------------- DipoleDipoleGalore ----------------

void DipoleDipoleGalore::sfReactionField(const json &j) { // Preliminary, needs to be checked!
    epsrf = j.at("epsrf");
    tableA = sfA.generate([&](double) { return 1.0; }, 0, 1);
    tableB = sfB.generate( [&](double q) { return -(2*(epsrf-epsr)/(2*epsrf+epsr))/epsr*q*q*q; },0,1 );
    calcDielectric = [&](double M2V) {
        if(epsrf > 1e10)
            return 1 + 3*M2V;
        if(fabs(epsrf-epsr) < 1e-6)
            return 2.25*M2V + 0.25 + 0.75*sqrt(9*M2V*M2V + 2*M2V + 1);
        if(fabs(epsrf-1.0) < 1e-6)
            return ( 2*M2V + 1 ) / ( 1 - M2V );
        return 0.5 * ( 2*epsrf - 1 + sqrt( -72*M2V*M2V*epsrf + 4*epsrf*epsrf + 4*epsrf + 1) ) / ( 3*M2V-1 ); // Needs to be checked!
        //return (6*M2V*epsrf + 2*epsrf + 1.0)/(1.0 + 2*epsrf - 3*M2V); // Is OK when epsr=1.0
        };
    selfenergy_prefactor = 2.0*(epsr - epsrf)/(2.0*epsrf + epsr); // Preliminary, needs to be checked!
}

void DipoleDipoleGalore::sfQ2potential(const json &j) { // Preliminary, needs to be checked!
    order = j.at("order");
    tableA = sfA.generate( [&](double q) { return qPochhammerSymbol(q,3,order);  },0,1 );
    tableB = sfB.generate([&](double) { return 0.0; }, 0, 1);
    calcDielectric = [&](double M2V) { return (2*M2V + 1.0)/(1.0 - M2V); };
    selfenergy_prefactor = -1.0;
}

void DipoleDipoleGalore::sfQ0potential(const json &j) { // Preliminary, needs to be checked!
    order = j.at("order");
    tableA = sfA.generate( [&](double q) { return dipoleDipoleQ2Help(q,0,order); },0,1 );
    tableB = sfB.generate( [&](double q) { return dipoleDipoleQ2Help(q,0,order,false); },0,1 );
    calcDielectric = [&](double M2V) { return 1 + 3*M2V; };
    selfenergy_prefactor = -1.0;
}

void DipoleDipoleGalore::sfFanourgakis(const json &) { // Preliminary, needs to be checked!
    tableA = sfA.generate( [&](double q) { return ( 1.0 + 14.0*pow(q,5) - 35.0*pow(q,6) + 20.0*pow(q,7) ); },0,1 );
    tableB = sfB.generate( [&](double q) { return 35.0*pow(q,5)*pow( 1.0 - q,2.0 ); },0,1 );
    calcDielectric = [&](double M2V) { return 1 + 3*M2V; };
    selfenergy_prefactor = 0.0; // Seems so but is it really correct? Check!
}

void DipoleDipoleGalore::sfFennell(const json &j) {
    alpha = j.at("alpha");
    double ar = alpha*rc;

    tableA = sfA.generate( [&](double q) {
      double kq = ar*q;
      double a = erfc(kq) + 4.0*kq*exp(-kq*kq)*(kq*kq + 1.5)/(3*std::sqrt(pc::pi));

      double kqc = ar; // using q=1 (i.e. r=Rc)
      double kqc2 = kqc*kqc;
      double ac = erfc(kqc) + 4.0*kqc*exp(-kqc2)*(kqc2 + 1.5)/(3*std::sqrt(pc::pi));

      double dac = -( 8.0 * ( kqc * ( kqc2*kqc2 + 1.5*kqc2 + 2.25 ) * exp(-kqc2) + 9.0 * 0.125 * std::sqrt(pc::pi) * erfc(kqc) ) ) / ( 3.0 * std::sqrt(pc::pi) );

      double q3 = q*q*q;
      return ( a - ac*q3 - (q-1.0)*dac*q3 ); },0,1 );
    tableB = sfB.generate( [&](double q) {
      double kq = ar*q;
      double b = 4.0*kq*kq*kq*exp(-kq*kq)/(3.0*std::sqrt(pc::pi));

      double kqc = ar; // using q=1 (i.e. r=Rc)
      double kqc2 = kqc*kqc;
      double bc = 4.0*kqc2*kqc*exp(-kqc2)/(3.0*std::sqrt(pc::pi));

      double dbc = -8.0*kqc2*kqc2*kqc*exp(-kqc2)/(3.0*std::sqrt(pc::pi));

      double q3 = q*q*q; // compensate for later multiplication with r^-3
      return ( b - bc*q3 - (q-1.0)*dbc*q3 ); },0,1 );
    calcDielectric = [&](double M2V) { double T = erf(ar) - (2 / (3 * sqrt(pc::pi))) * exp(-ar*ar) * ar * (ar*ar*ar*ar + 2.0 * ar*ar + 3.0);
        return (((T + 2.0) * M2V + 1.0)/ ((T - 1.0) * M2V + 1.0)); };
    selfenergy_prefactor = -0.5*( erfc(alpha*rc) + 2.0*alpha*rc/sqrt(pc::pi)*exp(-alpha*alpha*rc*rc) + (4.0/3.0)*pow(alpha*rc,3.0)/sqrt(pc::pi) );
}

void DipoleDipoleGalore::sfWolf(const json &j) {
    alpha = j.at("alpha");
    double ar = alpha*rc;

    tableA = sfA.generate( [&](double q) {
      double kq = ar*q;
      double a = erfc(kq) + 4.0*kq*exp(-kq*kq)*(kq*kq + 1.5)/(3*std::sqrt(pc::pi));

      double kqc = ar; // using q=1 (i.e. r=Rc)
      double ac = erfc(kqc) + 4.0*kqc*exp(-kqc*kqc)*(kqc*kqc + 1.5)/(3*std::sqrt(pc::pi));

      double q3 = q*q*q;
      return ( a - ac*q3 ); },0,1 );
    tableB = sfB.generate( [&](double q) {
      double kq = ar*q;
      double b = 4.0*kq*kq*kq*exp(-kq*kq)/(3.0*std::sqrt(pc::pi));

      double kqc = ar; // using q=1 (i.e. r=Rc)
      double bc = 4.0*kqc*kqc*kqc*exp(-kqc*kqc)/(3.0*std::sqrt(pc::pi));

      double q3 = q*q*q; // compensate for later multiplication with r^-3
      return ( b - bc*q3 ); },0,1 );
    calcDielectric = [&](double M2V) {
        double T = erf(ar) - (2 / (3 * sqrt(pc::pi))) * exp(-ar*ar) * (2.0 * ar*ar + 3.0); // check!
        return (((T + 2.0) * M2V + 1.0) / ((T - 1.0) * M2V + 1.0));
    };
    selfenergy_prefactor = -0.5*( erfc(ar) + 2.0*ar/sqrt(pc::pi)*exp(-ar*ar) + (4.0/3.0)*pow(ar,3.0)/sqrt(pc::pi) );
}

void DipoleDipoleGalore::sfEwald(const json &j) {
    alpha = j.at("alpha");
    double ar = alpha*rc;

    tableA = sfA.generate( [&](double q) { return ( erfc(ar*q) + 4.0*ar*q*exp(-ar*q*ar*q)*(ar*q*ar*q + 1.5)/(3*std::sqrt(pc::pi)) ); },0,1 );
    tableB = sfB.generate( [&](double q) { return ( 4.0*ar*q*ar*q*ar*q*exp(-ar*q*ar*q)/(3.0*std::sqrt(pc::pi)) ); },0,1 );
    calcDielectric = [&](double M2V) { return 1 + 3*M2V; };
    selfenergy_prefactor = -2.0/3.0*pow(alpha,3.0)/std::sqrt(pc::pi);
}

void DipoleDipoleGalore::sfPlain(const json &, double val) {
    tableA = sfA.generate([&](double) { return val; }, 0, 1);
    tableB = sfB.generate([&](double) { return 0; }, 0, 1);
    calcDielectric = [&](double M2V) { return (2.0 * M2V + 1.0) / (1.0 - M2V); };
    selfenergy_prefactor = 0.0;
}

void DipoleDipoleGalore::from_json(const json &j) {
    try {
        kappa = 0.0;
        type = j.at("type");
        rc = j.at("cutoff").get<double>();
        rc2 = rc * rc;
        rc1i = 1 / rc;
        epsr = j.at("epsr").get<double>();
        lB = pc::lB(epsr);

        depsdt = j.value("depsdt", -0.368 * pc::temperature / epsr);
        sfA.setTolerance(j.value("utol", 1e-5), j.value("ftol", 1e-2));
        sfB.setTolerance(j.value("utol", 1e-5), j.value("ftol", 1e-2));

        if (type == "plain")
            sfPlain(j, 1);
        else if (type == "none")
            sfPlain(j, 0);
        else if (type == "ewald")
            sfEwald(j);
        else if (type == "wolf")
            sfWolf(j);
        else if (type == "fennell")
            sfFennell(j);
        else if (type == "fanourgakis")
            sfFanourgakis(j);
        else if (type == "qpotential")
            sfQ0potential(j);
        else if (type == "q2potential")
            sfQ2potential(j);
        else if (type == "reactionfield")
            sfReactionField(j);

        if (tableA.empty() || tableB.empty())
            throw std::runtime_error(name + ": unknown dipole-dipole type '" + type + "'");

        // Set particle self-energy function. For reasons yet to be understood,
        // rc, rc2 etc. cannot be captured to reference, but must be hard-copied,
        // here into `factor`
        selfEnergy = [factor = selfenergy_prefactor * lB / (rc * rc2)](const Particle &a) -> double {
            if (a.hasExtension())
                return a.getExt().mulen * a.getExt().mulen * factor;
            return 0.0;
        };
    }

    catch (std::exception &e) {
        throw std::runtime_error(name + ": " + e.what());
    }
}

double DipoleDipoleGalore::dielectric_constant(double M2V) { return calcDielectric(M2V); }

void DipoleDipoleGalore::to_json(json &j) const {
    using namespace u8;
    j["epsr"] = epsr;
    j["T" + partial + epsilon_m + "/" + partial + "T"] = depsdt;
    j["lB"] = lB;
    j["cutoff"] = rc;
    j["type"] = type;
    if (type == "wolf" || type == "fennell" || type == "ewald")
        j["alpha"] = alpha;
    else if (type == "qpotential" || type == "q2potential")
        j["order"] = order;
    else if (type == "reactionfield")
        j["epsrf"] = epsrf;

    _roundjson(j, 5);
}

void Coulomb::to_json(json &j) const {
    j["epsr"] = pc::lB2epsr(lB);
    j["lB"] = lB;
}

void Coulomb::from_json(const json &j) { lB = pc::lB(j.at("epsr")); }

void DipoleDipole::to_json(json &j) const {
    j["epsr"] = pc::lB2epsr(lB);
    j["lB"] = lB;
}

void DipoleDipole::from_json(const json &j) { lB = pc::lB(j.at("epsr")); }

void FENE::from_json(const json &j) {
    k = j.at("stiffness");
    r02 = std::pow(double(j.at("maxsep")), 2);
    r02inv = 1 / r02;
}

void FENE::to_json(json &j) const { j = {{"stiffness", k}, {"maxsep", std::sqrt(r02)}}; }

void to_json(json &j, const PairPotentialBase &base) {
    base.name.empty() ? base.to_json(j) : base.to_json(j[base.name]);
}

void from_json(const json &j, PairPotentialBase &base) {
    try {
        if (not base.name.empty()) {
            if (j.count(base.name) == 1) {
                base.from_json(j.at(base.name));
                return;
            }
        }
        base.from_json(j);
    } catch (std::exception &e) {
        throw std::runtime_error("pairpotential error for " + base.name + ": " + e.what() + usageTip[base.name]);
    }
}

void from_json(const json &j, ParametersTable &m) {
    std::function<std::pair<double, double>(double, double, double, double)> mixerFunc;

    auto mixer = j.at("mixing").get<std::string>();
    if (mixer == "LB")
        m.mixer = ParametersTable::LB;
    if (mixer == "LBSW")
        m.mixer = ParametersTable::LBSW;
    if (mixer == "HE")
        m.mixer = ParametersTable::HE;

    size_t n = atoms.size(); // number of atom types
    switch (m.mixer) {
    case ParametersTable::LB:
        mixerFunc = [](double s1, double s2, double e1, double e2) {
            return std::pair<double, double>({(s1 + s2) / 2, std::sqrt(e1 * e2)});
        };
        m.s2.resize(n);  // not required...
        m.eps.resize(n); // ...but possible reduced mem. fragmentation
        break;
    case ParametersTable::LBSW:
        mixerFunc = [](double s1, double s2, double e1, double e2) {
            return std::pair<double, double>({s1 + s2, std::sqrt(e1 * e2)});
        };
        m.th.resize(n);  // not required...
        m.esw.resize(n); // ...but possible reduced mem. fragmentation
        break;
    case ParametersTable::HE:
        mixerFunc = [](double s1, double s2, double e1, double e2) {
            return std::pair<double, double>({s1 + s2, std::sqrt(e1 * e2)});
        };
        m.hd.resize(n);  // not required...
        m.ehe.resize(n); // ...but possible reduced mem. fragmentation
        break;
    default:
        throw std::runtime_error("unknown mixing rule");
    }

    for (auto &i : atoms)
        for (auto &j : atoms) {
            switch (m.mixer) {
            case ParametersTable::LB:
                double sigma, epsilon; // mixed values
                std::tie(sigma, epsilon) = mixerFunc(i.sigma, j.sigma, i.eps, j.eps);
                m.s2.set(i.id(), j.id(), sigma * sigma);
                m.eps.set(i.id(), j.id(), 4 * epsilon); // should already be in kT
                break;
            case ParametersTable::LBSW:
                double threshold, depth; // mixed values
                std::tie(threshold, depth) =
                    mixerFunc(0.5 * i.sigma + i.squarewell_threshold, 0.5 * j.sigma + j.squarewell_threshold,
                              i.squarewell_depth, j.squarewell_depth);
                m.th.set(i.id(), j.id(), threshold * threshold);
                m.esw.set(i.id(), j.id(), depth); // should already be in kT
                break;
            case ParametersTable::HE:
                double hdd, eh; // mixed values
                std::tie(hdd, eh) = mixerFunc(i.hdr, j.hdr, i.eps_hertz, j.eps_hertz);
                m.hd.set(i.id(), j.id(), hdd);
                m.ehe.set(i.id(), j.id(), eh); // should already be in kT
                break;
            default:
                throw std::runtime_error("unknown mixing rule");
            }
        }

    // custom eps/sigma for specific pairs
    if (j.count("custom") == 1) {
        auto &_j = j.at("custom");
        if (_j.is_object()) {
            for (auto it = _j.begin(); it != _j.end(); ++it) {
                auto v = words2vec<std::string>(it.key());
                if (v.size() == 2) {
                    auto it1 = findName(atoms, v[0]);
                    auto it2 = findName(atoms, v[1]);
                    if (it1 == atoms.end() or it2 == atoms.end())
                        throw std::runtime_error("unknown atom(s): ["s + v[0] + " " + v[1] + "]");
                    int id1 = it1->id();
                    int id2 = it2->id();

                    switch (m.mixer) {
                    case ParametersTable::LB:
                        m.s2.set(id1, id2, std::pow(it.value().at("sigma").get<double>(), 2));
                        m.eps.set(id1, id2, 4 * it.value().at("eps").get<double>() * 1.0_kJmol);
                        break;
                    case ParametersTable::LBSW:
                        m.th.set(id1, id2, it.value().at("sigma_sw").get<double>());
                        m.esw.set(id1, id2, it.value().at("eps_sw").get<double>() * 1.0_kJmol);
                        break;
                    case ParametersTable::HE:
                        m.hd.set(id1, id2, it.value().at("hdd").get<double>());
                        m.ehe.set(id1, id2, it.value().at("eps_hertz").get<double>() * 1.0_kJmol);
                        break;
                    default:
                        throw std::runtime_error("unknown mixing rule");
                    }
                } else {
                    switch (m.mixer) {
                    case ParametersTable::LB:
                        throw std::runtime_error(
                            "custom epsilon/sigma parameters require exactly two space-separated atoms");
                    case ParametersTable::LBSW:
                        throw std::runtime_error(
                            "custom eps_sw/sigma_sw parameters require exactly two space-separated atoms");
                    case ParametersTable::HE:
                        throw std::runtime_error(
                            "custom eps_hertz/hdd parameters require exactly two space-separated atoms");
                    default:
                        throw std::runtime_error("unknown mixing rule");
                    }
                }
            }
        } else {
            switch (m.mixer) {
            case ParametersTable::LB:
                throw std::runtime_error("custom sigma/epsilon syntax error");
            case ParametersTable::LBSW:
                throw std::runtime_error("custom eps_sw/sigma_sw syntax error");
            case ParametersTable::HE:
                throw std::runtime_error("custom eps_hertz/hdd syntax error");
            default:
                throw std::runtime_error("unknown mixing rule");
            }
        }
    }
}
void to_json(json &j, const ParametersTable &m) {
    auto &_j = j["custom"];
    switch (m.mixer) {
    case ParametersTable::LB:
        j["mixing"] = "LB";
        j["epsilon unit"] = "kJ/mol";
        for (size_t i = 0; i < m.eps.size(); i++)
            for (size_t j = 0; j < m.eps.size(); j++)
                if (i >= j) {
                    auto str = atoms[i].name + " " + atoms[j].name;
                    _j[str] = {{"eps", m.eps(i, j) / 4.0_kJmol}, {"sigma", std::sqrt(m.s2(i, j))}};
                    _roundjson(_j[str], 5);
                }
        break;
    case ParametersTable::LBSW:
        j["mixing"] = "LBSW";
        j["depth unit"] = "kJ/mol";
        for (size_t i = 0; i < m.esw.size(); i++)
            for (size_t j = 0; j < m.esw.size(); j++)
                if (i >= j) {
                    auto str = atoms[i].name + " " + atoms[j].name;
                    _j[str] = {{"eps_sw", m.esw(i, j) / 1.0_kJmol}, {"sigma_sw", m.th(i, j)}};
                    _roundjson(_j[str], 5);
                }
        break;
    case ParametersTable::HE:
        j["mixing"] = "HE";
        j["eps_hertz unit"] = "kJ/mol";
        for (size_t i = 0; i < m.ehe.size(); i++)
            for (size_t j = 0; j < m.ehe.size(); j++)
                if (i >= j) {
                    auto str = atoms[i].name + " " + atoms[j].name;
                    _j[str] = {{"eps_hertz", m.ehe(i, j) / 1.0_kJmol}, {"hdd", m.hd(i, j)}};
                    _roundjson(_j[str], 5);
                }
        break;
    default:
        throw std::runtime_error("unknown mixing rule");
    }
}

void SASApotential::from_json(const json &j) {
    assertKeys(j, {"shift", "molarity", "radius"});
    shift = j.value("shift", true);
    conc = j.at("molarity").get<double>() * 1.0_molar;
    proberadius = j.value("radius", 1.4) * 1.0_angstrom;
}

void SASApotential::to_json(json &j) const {
    j["molarity"] = conc / 1.0_molar;
    j["radius"] = proberadius / 1.0_angstrom;
    j["shift"] = shift;
}

double SASApotential::area(double R, double r, double d_squared) const {
    R += proberadius;
    r += proberadius;
    double area = 4 * pc::pi * (R * R + r * r); // full volume of both spheres
    double offset = (shift ? area : 0);
    if (d_squared > (R + r) * (R + r))
        return area - offset;
    if (r > R)
        std::swap(r, R);
    double d = sqrt(d_squared);
    if (d + r <= R)
        return 4 * pc::pi * R * R - offset;          // full volume of biggest sphere
    double h1 = (r - R + d) * (r + R - d) / (2 * d); // height of spherical caps
    double h2 = (R - r + d) * (R + r - d) / (2 * d); // comprising intersecting lens
    return area - 2 * pc::pi * (R * h1 + r * h2) - offset;
}

void CustomPairPotential::from_json(const json &j) {
    Rc2 = j.value("cutoff", pc::infty);
    Rc2 = Rc2 * Rc2;
    jin = j;
    auto &_j = jin["constants"];
    if (_j == nullptr)
        _j = json::object();
    _j["e0"] = pc::e0;
    _j["kB"] = pc::kB;
    _j["kT"] = pc::kT();
    _j["Nav"] = pc::Nav;
    _j["Rc"] = std::sqrt(Rc2);
    _j["T"] = pc::temperature;
    expr.set(jin, {{"r", &d->r}, {"q1", &d->q1}, {"q2", &d->q2}, {"s1", &d->s1}, {"s2", &d->s2}});
}

void CustomPairPotential::to_json(json &j) const {
    j = jin;
    if (std::isfinite(Rc2))
        j["cutoff"] = std::sqrt(Rc2);
}


// =============== Dummy ===============

Dummy::Dummy() { name = "dummy"; }
void Dummy::from_json(const json &) {}
void Dummy::to_json(json &) const {}

// =============== LennardJones ===============

void LennardJones::initPairMatrices() {
    TCombinatorFunc func_sigma_squared, func_epsilon_quadrupled;
    auto faunus_logger = spdlog::get("faunus");
    if (combination_rule == LorentzBerthelot) {
        faunus_logger->debug("Lorentz-Berthelot combination rules in effect for the {} potential.", name);
        func_sigma_squared = &PairMixer::combArithmeticSquared;
        func_epsilon_quadrupled = [](double e1, double e2) -> double { return 4 * std::sqrt(e1 * e2); };
    } else {
        throw std::runtime_error("Undefined combination rule for the Lennard-Jones potential.");
        // todo error
    }

    sigma_squared = PairMixer([](const AtomData &a) -> double { return a.sigma; }, func_sigma_squared).
            createPairMatrix(atoms, custom_pairs);
    epsilon_quadrupled = PairMixer([](const AtomData &a) -> double { return a.eps; }, func_epsilon_quadrupled).
            createPairMatrix(atoms, custom_pairs);
    faunus_logger->debug("Pair matrices for {} sigma ({}×{}) and epsilon ({}×{}) created using {} custom pairs.", name,
                         sigma_squared->rows(), sigma_squared->cols(),
                         epsilon_quadrupled->rows(), epsilon_quadrupled->cols(), custom_pairs.size());
}

void LennardJones::from_json(const json &j) {
    auto mixing = j.at("mixing").get<std::string>();
    if (mixing == "LB") {
        combination_rule = LorentzBerthelot;
    } else {
        throw std::runtime_error("unknown mixing rule for Lennard-Jones potential");
        // todo error
    }
    if (j.count("custom") == 1) {
        custom_pairs = j;
    }
    initPairMatrices();
}

void LennardJones::to_json(json &j) const {
    if (combination_rule == LorentzBerthelot) {
        j["mixing"] = "LB";
    } else {
        //todo error
    }
    if (!custom_pairs.empty()) {
        j["custom"] = custom_pairs;
    }
}


// =============== HardSphere ===============

void HardSphere::initPairMatrices() {
    auto faunus_logger = spdlog::get("faunus");
    faunus_logger->debug("Arithmetic mean combination rule in effect for the {} potential.", name);
    sigma_squared = PairMixer([](const AtomData &a) -> double { return a.sigma; }, &PairMixer::combArithmeticSquared).
            createPairMatrix(atoms, custom_pairs);
    faunus_logger->debug("Pair matrix for {} sigma ({}×{}) created using {} custom pairs.", name,
                         sigma_squared->rows(), sigma_squared->cols(), custom_pairs.size());
}

void HardSphere::from_json(const json &j) {
    if (j.count("custom") == 1) {
        custom_pairs = j;
    }
    initPairMatrices();
}

void HardSphere::to_json(json &j) const {
    if (!custom_pairs.empty()) {
        j["custom"] = custom_pairs;
    }
}


// =============== Hertz ===============

void Hertz::initPairMatrices() {
    auto faunus_logger = spdlog::get("faunus");
    faunus_logger->debug("Hertz combination rules in effect for the {} potential.", name);
    hydrodynamic_diameter = PairMixer([](const AtomData &a) -> double { return a.hdr; }, &PairMixer::combSum).
            createPairMatrix(atoms, custom_pairs);
    epsilon_hertz = PairMixer([](const AtomData &a) -> double { return a.eps_hertz; }, &PairMixer::combGeometric).
            createPairMatrix(atoms, custom_pairs);
    faunus_logger->debug(
            "Pair matrix for {} hydrodynamic radius ({}×{}) and epsilon ({}×{}) created using {} custom pairs.", name,
            hydrodynamic_diameter->rows(), hydrodynamic_diameter->cols(), epsilon_hertz->rows(), epsilon_hertz->cols(),
            custom_pairs.size());
}

void Hertz::from_json(const json &j) {
    if (j.count("custom") == 1) {
        custom_pairs = j;
    }
    initPairMatrices();
}

void Hertz::to_json(json &j) const {
    if (!custom_pairs.empty()) {
        j["custom"] = custom_pairs;
    }
}

void SquareWell::to_json(json &j) const { j = *m; }
void SquareWell::from_json(const json &j) {
    *m = j;
    if (m->th.size() == 0)
        throw std::runtime_error("unknown mixing rule for Square-well potential");
}
void Polarizability::from_json(const json &j) {
    epsr = j.at("epsr").get<double>();
    double lB = pc::lB(epsr);
    for (auto &i : atoms) {
        for (auto &j : atoms) {
            m_neutral->set(i.id(), j.id(), -3 * i.alphax * pow(0.5 * i.sigma, 3) * j.alphax * pow(0.5 * j.sigma, 3));
            m_charged->set(i.id(), j.id(),
                           -lB / 2 *
                               (pow(i.charge, 2) * j.alphax * pow(0.5 * j.sigma, 3) +
                                pow(j.charge, 2) * i.alphax * pow(0.5 * i.sigma, 3)));
        }
    }
}

//----------------- FunctorPotential ---------------------

void FunctorPotential::registerSelfEnergy(PairPotentialBase *pot) {
    if (pot->selfEnergy) {
        self_energy_vector.push_back(pot->selfEnergy);
        assert(self_energy_vector.back());
    }
}

FunctorPotential::uFunc FunctorPotential::combineFunc(const json &j) {
    uFunc u = [](const Particle &, const Particle &, const Point &) { return 0.0; };
    if (j.is_array()) {
        for (auto &i : j) { // loop over all defined potentials in array
            if (i.is_object() and (i.size() == 1)) {
                for (auto it : i.items()) {
                    uFunc _u = nullptr;
                    try {
                        if (it.key() == "custom")
                            _u = CustomPairPotential() = it.value();

                        // add Coulomb potential and self-energy
                        // terms if not already added
                        else if (it.key() == "coulomb") {
                            _u = std::get<0>(potlist) = i;
                            if (not have_monopole_self_energy) {
                                registerSelfEnergy(&std::get<0>(potlist));
                                have_monopole_self_energy = true;
                            }
                        } else if (it.key() == "cos2")
                            _u = std::get<1>(potlist) = i;
                        else if (it.key() == "polar")
                            _u = std::get<2>(potlist) = i;
                        else if (it.key() == "hardsphere")
                            _u = std::get<3>(potlist) = i;
                        else if (it.key() == "lennardjones")
                            _u = std::get<4>(potlist) = i;
                        else if (it.key() == "repulsionr3")
                            _u = std::get<5>(potlist) = i;
                        else if (it.key() == "sasa")
                            _u = std::get<6>(potlist) = i;
                        else if (it.key() == "wca")
                            _u = std::get<7>(potlist) = i;
                        else if (it.key() == "pm")
                            _u = std::get<8>(potlist) = it.value();
                        else if (it.key() == "pmwca")
                            _u = std::get<9>(potlist) = it.value();
                        else if (it.key() == "hertz")
                            _u = std::get<10>(potlist) = i;
                        else if (it.key() == "squarewell")
                            _u = std::get<11>(potlist) = i;
                        else if (it.key() == "dipoledipole") {
                            isotropic = false; // potential is now angular dependent
                            _u = std::get<12>(potlist) = i;
                            if (not have_dipole_self_energy) {
                                registerSelfEnergy(&std::get<12>(potlist));
                                have_dipole_self_energy = true;
                            }
                        } else if (it.key() == "stockmayer") {
                            _u = std::get<13>(potlist) = it.value();
                            isotropic = false; // potential is now angular dependent
                            if (not have_dipole_self_energy) {
                                registerSelfEnergy(&std::get<13>(potlist));
                                have_dipole_self_energy = true;
                            }
                        }
                        // place additional potentials here...
                    } catch (std::exception &e) {
                        throw std::runtime_error("Error adding energy '" + it.key() + "': " + e.what() +
                                                 usageTip[it.key()]);
                    }

                    if (_u != nullptr) // if found, sum them into new function object
                        u = [u, _u](const Particle &a, const Particle &b, const Point &r) {
                            return u(a, b, r) + _u(a, b, r);
                        };
                    else
                        throw std::runtime_error("unknown pair-potential: " + it.key());
                }
            }
        }
    } else
        throw std::runtime_error("dictionary of potentials required");

    // set self energy function
    if (self_energy_vector.empty())
        selfEnergy = nullptr;
    else
        // why must self_energy_vector be copied? Using [=] or [&]
        // lambda capturing causes bad access on copied objects
        selfEnergy = [vec = self_energy_vector](const Particle &p) {
            double sum = 0;
            for (auto &func : vec) {
                assert(func);
                sum += func(p);
            }
            return sum;
        };
    return u;
}

void FunctorPotential::to_json(json &j) const {
    j = _j;
    j["selfenergy"] = {{"monopole", have_monopole_self_energy}, {"dipole", have_dipole_self_energy}};
}

void FunctorPotential::from_json(const json &j) {
    _j = j;
    umatrix = decltype(umatrix)(atoms.size(), combineFunc(j.at("default")));
    for (auto it = j.begin(); it != j.end(); ++it) {
        auto atompair = words2vec<std::string>(it.key()); // is this for a pair of atoms?
        if (atompair.size() == 2) {
            auto ids = names2ids(atoms, atompair);
            umatrix.set(ids[0], ids[1], combineFunc(it.value()));
        }
    }
}

//---------------- TabulatedPotential ---------------------

void TabulatedPotential::from_json(const json &j) {
    FunctorPotential::from_json(j);

    // if user specifies an anisotropic potential, make sure to bail out
    if (not isotropic)
        throw std::runtime_error("only isotropic pair-potentials can be splined");

    tblt.setTolerance(j.value("utol", 1e-5), j.value("ftol", 1e-2));
    double u_at_rmin = j.value("u_at_rmin", 20);
    double u_at_rmax = j.value("u_at_rmax", 1e-6);
    hardsphere = j.value("hardsphere", false);

    // build matrix of spline data, each element corresponding
    // to a pair of atom types
    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t k = 0; k <= i; ++k) {
            if (atoms[i].implicit == false and atoms[k].implicit == false) {
                Particle a = atoms.at(i);
                Particle b = atoms.at(k);
                double rmin2 = .5 * (atoms[i].sigma + atoms[k].sigma);
                rmin2 = rmin2 * rmin2;
                double rmax2 = rmin2 * 100;
                auto it = j.find("cutoff_g2g");
                if (j.count("rmax") == 1) {
                    rmax2 = std::pow(j.at("rmax").get<double>(), 2);
                } else if (it != j.end()) {
                    if (it->is_number())
                        rmax2 = std::pow(it->get<double>(), 2);
                    else if (it->is_object())
                        rmax2 = std::pow(it->at("default").get<double>(), 2);
                }

                // adjust lower splining distance to match
                // the given energy threshold (u_at_min2)
                double dr = 1e-2;
                while (rmin2 >= dr) {
                    double u = std::fabs(this->umatrix(i, k)(a, b, {0, 0, sqrt(rmin2)}));
                    if (u > u_at_rmin * 1.1)
                        rmin2 = rmin2 + dr;
                    else if (u < u_at_rmin / 1.1)
                        rmin2 = rmin2 - dr;
                    else
                        break;
                }

                assert(rmin2 >= 0);

                while (rmax2 >= dr) {
                    double u = std::fabs(this->umatrix(i, k)(a, b, {0, 0, sqrt(rmax2)}));
                    if (u > u_at_rmax)
                        rmax2 = rmax2 + dr;
                    else
                        break;
                }

                assert(rmin2 < rmax2);

                Ttable knotdata = tblt.generate(
                    [&](double r2) {
                        return this->umatrix(i, k)(a, b, {0, 0, sqrt(r2)});
                    },
                    rmin2, rmax2);

                // assert if potential is negative for r<rmin
                if (tblt.eval(knotdata, knotdata.rmin2 + dr) < 0)
                    knotdata.isNegativeBelowRmin = true;

                tmatrix.set(i, k, knotdata);
                if (j.value("to_disk", false)) {
                    std::ofstream f(atoms[i].name + "-" + atoms[k].name + "_tabulated.dat"); // output file
                    f << "# r splined exact\n";
                    Point r = {dr, 0, 0}; // variable distance vector between particle a and b
                    for (; r.x() < sqrt(rmax2); r.x() += dr)
                        f << r.x() << " " << operator()(a, b, r) << " " << this->umatrix(i, k)(a, b, r) << "\n";
                }
            }
        }
    }
}
} // namespace Potential
} // namespace Faunus
