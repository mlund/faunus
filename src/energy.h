#pragma once

#include "core.h"
#include "geometry.h"
#include "space.h"
#include "potentials.h"
#include "multipole.h"
#include "penalty.h"
#include "mpi.h"
#include <Eigen/Dense>
#include <set>

#ifdef FAU_POWERSASA
#include <power_sasa.h>
#endif

namespace Faunus {
    namespace Energy {

        class Energybase {
            public:
                enum keys {OLD, NEW, NONE};
                keys key=NONE;
                std::string name;
                std::string cite;
                virtual double energy(Change&)=0; //!< energy due to change
                inline virtual void to_json(json &j) const {}; //!< json output
                inline virtual void sync(Energybase*, Change&) {}
        };

        void to_json(json &j, const Energybase &base) {
            assert(!base.name.empty());
            if (!base.cite.empty())
                j[base.name]["reference"] = base.cite;
            base.to_json( j[base.name] );
        } //!< Converts any energy class to json object

        /**
         * This holds Ewald setup and must *not* depend on particle type, nor depend on Space
         */
        struct EwaldData {
            typedef std::complex<double> Tcomplex;
            Eigen::Matrix3Xd kVectors; // k-vectors, 3xK
            Eigen::VectorXd Aks;       // 1xK, to minimize computational effort (Eq.24,DOI:10.1063/1.481216)
            Eigen::VectorXcd Qion, Qdip; // 1xK
            double alpha, rc, kc, check_k2_zero, lB;
            double const_inf, eps_surf;
            bool spherical_sum=true;
            bool ipbc=false;
            int kVectorsInUse=0;
            Point L; //!< Box dimensions

            void update(const Point &box) {
                L = box;
                int kcc = std::ceil(kc);
                check_k2_zero = 0.1*std::pow(2*pc::pi/L.maxCoeff(), 2);
                int kVectorsLength = (2*kcc+1) * (2*kcc+1) * (2*kcc+1) - 1;
                if (kVectorsLength == 0) {
                    kVectors.resize(3,1); 
                    Aks.resize(1);
                    kVectors.col(0) = Point(1,0,0); // Just so it is not the zero-vector
                    Aks[0] = 0;
                    kVectorsInUse = 1;
                    Qion.resize(1);
                    Qdip.resize(1);
                } else {
                    double kc2 = kc*kc;
                    kVectors.resize(3, kVectorsLength); 
                    Aks.resize(kVectorsLength);
                    kVectorsInUse = 0;
                    kVectors.setZero();
                    Aks.setZero();
                    int startValue = 1 - int(ipbc);
                    double factor = 1;
                    for (int kx = 0; kx <= kcc; kx++) {
                        double dkx2 = double(kx*kx);
                        for (int ky = -kcc*startValue; ky <= kcc; ky++) {
                            double dky2 = double(ky*ky);
                            for (int kz = -kcc*startValue; kz <= kcc; kz++) {
                                double factor = 1.0;
                                if(kx > 0)
                                    factor *= 2;
                                if(ky > 0 && ipbc)
                                    factor *= 2;
                                if(kz > 0 && ipbc)
                                    factor *= 2;
                                double dkz2 = double(kz*kz);
                                Point kv = 2*pc::pi*Point(kx/L.x(),ky/L.y(),kz/L.z());
                                double k2 = kv.dot(kv);
                                if (k2 < check_k2_zero) // Check if k2 != 0
                                    continue;
                                if (spherical_sum)
                                    if( (dkx2/kc2) + (dky2/kc2) + (dkz2/kc2) > 1)
                                        continue;
                                kVectors.col(kVectorsInUse) = kv; 
                                Aks[kVectorsInUse] = factor*std::exp(-k2/(4*alpha*alpha))/k2;
                                kVectorsInUse++;
                            }
                        }
                    }
                    Qion.resize(kVectorsInUse);
                    Qdip.resize(kVectorsInUse);
                    Aks.conservativeResize(kVectorsInUse);
                    kVectors.conservativeResize(3,kVectorsInUse);
                }
            }
        };

        void from_json(const json &j, EwaldData &d) {
            d.alpha = j.at("alpha");
            d.rc = j.at("cutoff");
            d.kc = j.at("kcutoff");
            d.ipbc = j.value("ipbc", false);
            d.spherical_sum = j.value("spherical_sum", true);
            d.lB = pc::lB( j.at("epsr") );
	    d.eps_surf = j.value("epss", 0.0);
            d.const_inf = (d.eps_surf < 1) ? 0 : 1; // if unphysical (<1) use epsr infinity for surrounding medium
        }

        void to_json(json &j, const EwaldData &d) {
            j = {{"lB", d.lB}, {"ipbc", d.ipbc}, {"epss", d.eps_surf},
                {"alpha", d.alpha}, {"cutoff", d.rc}, {"kcutoff", d.kc},
                {"wavefunctions", d.kVectors.cols()}, {"spherical_sum", d.spherical_sum}};
        }

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] Ewald - EwaldData") 
        {
            using doctest::Approx;

            EwaldData data = R"({
                "ipbc": false, "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "kcutoff": 11.0, "spherical_sum": true, "cutoff": 5.0})"_json;

            data.update( Point(10,10,10) );

            CHECK(data.ipbc == false);
            CHECK(data.const_inf == 1);
            CHECK(data.alpha == 0.894427190999916);
            CHECK(data.kVectors.cols() == 2975); 
            CHECK(data.Qion.size() == data.kVectors.cols());

            data.ipbc=true;
            data.update( Point(10,10,10) );
            CHECK(data.kVectors.cols() == 846); 
            CHECK(data.Qion.size() == data.kVectors.cols());
        }
#endif
 
        /** @brief recipe or policies for ion-ion ewald */
        template<class Tspace, bool eigenopt=false /** use Eigen matrix ops where possible */>
            struct PolicyIonIon  {
                typedef typename Tspace::Tpvec::iterator iter;
                Tspace *spc;
                Tspace *old=nullptr; // set only if key==NEW at first call to `sync()`

                PolicyIonIon(Tspace &spc) : spc(&spc) {}

                void updateComplex(EwaldData &data) const {
                    if (eigenopt)
                        if (data.ipbc==false) {
                            auto pos = asEigenMatrix(spc->p.begin(), spc->p.end(), &Tspace::Tparticle::pos); //  Nx3
                            auto charge = asEigenVector(spc->p.begin(), spc->p.end(), &Tspace::Tparticle::charge); // Nx1
                            Eigen::MatrixXd kr = pos.matrix() * data.kVectors; // Nx3 * 3xK = NxK
                            data.Qion.real() = (kr.array().cos().colwise()*charge).colwise().sum();
                            data.Qion.imag() = kr.array().sin().colwise().sum();
                            return;
                        }
                    for (int k=0; k<data.kVectors.cols(); k++) {
                        const Point& kv = data.kVectors.col(k);
                        EwaldData::Tcomplex Q(0,0);
                        if (data.ipbc)
                            for (auto &i : spc->p)
                                Q += kv.cwiseProduct(i.pos).array().cos().prod() * i.charge;
                        else
                            for (auto &i : spc->p) {
                                double dot = kv.dot(i.pos);
                                Q += i.charge * EwaldData::Tcomplex( std::cos(dot), std::sin(dot) );
                            }
                        data.Qion[k] = Q;
                    }
                } //!< Update all k vectors

                void updateComplex(EwaldData &data, iter begin, iter end) const {
                    assert(old!=nullptr);
                    assert(spc->p.size() == old->p.size());
                    size_t ibeg = std::distance(spc->p.begin(), begin); // it->index
                    size_t iend = std::distance(spc->p.begin(), end);   // it->index
                    for (int k=0; k<data.kVectors.cols(); k++) {
                        auto& Q = data.Qion[k];
                        Point q = data.kVectors.col(k);
                        if (data.ipbc)
                            for (size_t i=ibeg; i<=iend; i++) {
                                Q +=  q.cwiseProduct( spc->p[i].pos ).array().cos().prod() * spc->p[i].charge;
                                Q -=  q.cwiseProduct( old->p[i].pos ).array().cos().prod() * old->p[i].charge;
                            }
                        else
                            for (size_t i=ibeg; i<=iend; i++) {
                                double _new = q.dot(spc->p[i].pos);
                                double _old = q.dot(old->p[i].pos);
                                Q += spc->p[i].charge * EwaldData::Tcomplex( std::cos(_new), std::sin(_new) );
                                Q -= old->p[i].charge * EwaldData::Tcomplex( std::cos(_old), std::sin(_old) );
                            }
                    }
                } //!< Optimized update of k subset. Require access to old positions through `old` pointer

                double selfEnergy(const EwaldData &d) {
                    double E = 0;
                    for (auto& i : spc->p)
                        E += i.charge * i.charge;
                    return -d.alpha*E / std::sqrt(pc::pi) * d.lB;
                }

                double surfaceEnergy(const EwaldData &d) {
                    if (d.const_inf < 0.5)
                        return 0;
                    Point qr(0,0,0);
                    for (auto &i : spc->p)
                        qr += i.charge*i.pos;
                    return d.const_inf * 2 * pc::pi / ( (2*d.eps_surf+1) * spc->geo.getVolume() ) * qr.dot(qr) * d.lB;
                }

                double reciprocalEnergy(const EwaldData &d) {
                    double E = 0;
                    if (eigenopt) // known at compile time
                        E = d.Aks.cwiseProduct( d.Qion.cwiseAbs2() ).sum();
                    else
                        for (size_t k=0; k<d.Qion.size(); k++)
                            E += d.Aks[k] * std::norm( d.Qion[k] );
                    return 2 * pc::pi / spc->geo.getVolume() * E * d.lB;
                }
            };

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] Ewald - IonIonPolicy") 
        {
            using doctest::Approx;
            typedef Space<Geometry::Cuboid, Particle<Charge,Dipole>> Tspace;

            Tspace spc;
            spc.p.resize(2);
            spc.geo  = R"( {"length": 10} )"_json;
            spc.p[0] = R"( {"pos": [0,0,0], "q": 1.0} )"_json;
            spc.p[1] = R"( {"pos": [1,0,0], "q": -1.0} )"_json;
  
            PolicyIonIon<Tspace> ionion(spc);
            EwaldData data = R"({
                "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "kcutoff": 11.0, "spherical_sum": true, "cutoff": 5.0})"_json;

            data.ipbc = false; // PBC Ewald (http://dx.doi.org/10.1063/1.481216)
            data.update( spc.geo.getLength() );
            ionion.updateComplex( data );
            CHECK( ionion.selfEnergy(data) == Approx(-1.0092530088080642*data.lB) );
            CHECK( ionion.surfaceEnergy(data) == Approx(0.0020943951023931952*data.lB) );
            CHECK( ionion.reciprocalEnergy(data) == Approx(0.21303063979675319*data.lB) );

            data.ipbc = true; // IPBC Ewald
            data.update( spc.geo.getLength() );
            ionion.updateComplex( data );
            CHECK( ionion.selfEnergy(data) == Approx(-1.0092530088080642*data.lB) );
            CHECK( ionion.surfaceEnergy(data) == Approx(0.0020943951023931952*data.lB) );
            CHECK( ionion.reciprocalEnergy(data) == Approx(0.0865107467*data.lB) );
        }
#endif

        /** @brief Ewald summation reciprocal energy */
        template<class Tspace, class Policy=PolicyIonIon<Tspace>>
            class Ewald : public Energybase {
                private:
                    EwaldData data;
                    Policy policy;
                public:
                    Tspace& spc;

                    Ewald(const json &j, Tspace &spc) : spc(spc), policy(spc) {
                        name = "ewald";
                        data = j;
                        data.update( spc.geo.getLength() );
                        policy.updateComplex(data); // brute force. todo: be selective
                    }

                    double energy(Change &change) override {
                        double u=0;
                        if (!change.empty()) {
                            // If the state is NEW (trial state), then update all k-vectors
                            if (key==NEW) {
                                if (change.all || change.dV) {       // everything changes
                                    data.update( spc.geo.getLength() );
                                    policy.updateComplex(data);    // update all (expensive!)
                                }
                                else {
                                    if (change.groups.size()==1) { // exactly one group is moved
                                        auto& d = change.groups[0];
                                        auto& g = spc.groups[d.index];
                                        if (d.atoms.size()==1)     // exactly one atom is moved
                                            policy.updateComplex(data, g.begin()+d.atoms[0], g.begin()+d.atoms[0]);
                                        else
                                            policy.updateComplex(data, g.begin(), g.end());
                                    } else
                                        policy.updateComplex(data);
                                }
                            }
                            u = policy.selfEnergy(data) + policy.surfaceEnergy(data) + policy.reciprocalEnergy(data); 
                        }
                        return u;
                    }

                    void sync(Energybase *basePtr, Change &change) override {
                        auto other = dynamic_cast<decltype(this)>(basePtr);
                        assert(other);
                        if (other->key==OLD)
                            policy.old = &(other->spc); // give NEW access to OLD space for optimized updates
                        data = other->data; // copy everything!

                    } //!< Called after a move is rejected/accepted as well as before simulation

                    void to_json(json &j) const override {
                        j = data;
                    }
            };

        template<typename Tspace>
            class Isobaric : public Energybase {
                private:
                    Tspace& spc;
                    double P; // P/kT
                public:
                    Isobaric(const json &j, Tspace &spc) : spc(spc) {
                        name = "isobaric";
                        cite = "Frenkel & Smith 2nd Ed (Eq. 5.4.13)";
                        P = j.value("P/mM", 0.0) * 1.0_mM;
                        if (P<1e-10) {
                            P = j.value("P/Pa", 0.0) * 1.0_Pa;
                            if (P<1e-10)
                                P = j.at("P/atm").get<double>() * 1.0_atm;
                        }
                    }
                    double energy(Change &change) override {
                        if (change.dV || change.all) {
                            double V = spc.geo.getVolume();
                            size_t N=0;
                            for (auto &g : spc.groups)
                                if (!g.empty()) {
                                    if (g.atomic)
                                        N += g.size();
                                    else
                                        N++;
                                }
                            return P*V-(N+1)*std::log(V);
                        } else return 0; 
                    }
                    void to_json(json &j) const override {
                        j["P/atm"] = P / 1.0_atm;
                        j["P/mM"] = P / 1.0_mM;
                        j["P/Pa"] = P / 1.0_Pa;
                        _roundjson(j,5);
                    }
            };

        template<typename Tspace>
            class ExternalPotential : public Energybase {
                protected:
                    typedef typename Tspace::Tpvec Tpvec;
                    typedef typename Tspace::Tparticle Tparticle;
                    bool COM=false; // apply on center-of-mass
                    Tspace& spc;
                    std::set<int> molids; // molecules to act upon
                    std::function<double(const Tparticle&)> func=nullptr; // energy of single particle
                    std::vector<std::string> _names;

                    template<class Tparticle>
                        double _energy(const Group<Tparticle> &g) const {
                            double u=0;
                            if (molids.find(g.id)!=molids.end()) {
                                if (COM) { // apply only to center of mass
                                    Tparticle dummy;
                                    dummy.pos = g.cm;
                                    u = func(dummy);
                                } else {
                                    for (auto &p : g) {
                                        u += func(p);
                                        if (std::isnan(u))
                                            break;
                                    }
                                }
                            }
                            return u;
                        } //!< External potential on a single particle
                public:
                    ExternalPotential(const json &j, Tspace &spc) : spc(spc) {
                        name="external";
                        COM = j.value("com", false);
                        _names = j.at("molecules").get<decltype(_names)>(); // molecule names
                        auto _ids = names2ids(molecules<Tpvec>, _names);     // names --> molids
                        molids = std::set<int>(_ids.begin(), _ids.end());    // vector --> set
                        if (molids.empty() || molids.size()!=_names.size() )
                            throw std::runtime_error(name + ": molecule list is empty");
                    }

                    double energy(Change &change) override {
                        assert(func!=nullptr);
                        double u=0;
                        if (change.dV || change.all) {
                            for (auto &g : spc.groups) { // check all groups
                                u += _energy(g);
                                if (std::isnan(u))
                                    break;
                            }
                        } else
                            for (auto &d : change.groups) {
                                auto &g = spc.groups.at(d.index); // check specified groups
                                if (d.all || COM)  // check all atoms in group
                                    u += _energy(g);
                                else {       // check only specified atoms in group
                                    if (molids.find(g.id)!=molids.end())
                                        for (auto i : d.atoms)
                                            u += func( *(g.begin()+i) );
                                }
                                if (std::isnan(u))
                                    break;
                            }
                        return u;
                    }

                    void to_json(json &j) const override {
                        j["molecules"] = _names;
                        j["com"] = COM;
                    }
            }; //!< Base class for external potentials, acting on particles

        template<typename Tspace, typename base=ExternalPotential<Tspace>>
            class Confine : public base {
                public:
                    enum Variant {sphere, cylinder, cuboid, none};
                    Variant type=none;

                private:
                    Point origo={0,0,0}, dir={1,1,1};
                    Point low, high;
                    double radius, k;
                    bool scale=false;
                    std::map<std::string, Variant> m = {
                        {"sphere", sphere}, {"cylinder", cylinder}, {"cuboid", cuboid}
                    };

                public:
                    Confine(const json &j, Tspace &spc) : base(j,spc) {
                        base::name = "confine";
                        k = value_inf(j, "k") * 1.0_kJmol; // get floating point; allow inf/-inf
                        type = m.at( j.at("type") );

                        if (type==sphere || type==cylinder) {
                            radius = j.at("radius");
                            origo = j.value("origo", origo);
                            scale = j.value("scale", scale);
                            if (type==cylinder)
                                dir = {1,1,0};
                            base::func = [&radius=radius, origo=origo, k=k, dir=dir](const typename base::Tparticle &p) {
                                double d2 = (origo-p.pos).cwiseProduct(dir).squaredNorm() - radius*radius;
                                if (d2>0)
                                    return 0.5*k*d2;
                                return 0.0;
                            };

                            // If volume is scaled, also scale the confining radius by adding a trigger
                            // to `Space::scaleVolume()`
                            if (scale)
                                spc.scaleVolumeTriggers.push_back( [&radius=radius](Tspace &spc, double Vold, double Vnew) {
                                        radius *= std::cbrt(Vnew/Vold); } );
                        }

                        if (type==cuboid) {
                            low = j.at("low").get<Point>();
                            high = j.at("high").get<Point>();
                            base::func = [low=low, high=high, k=k](const typename base::Tparticle &p) {
                                double u=0;
                                Point d = low-p.pos;
                                for (int i=0; i<3; ++i)
                                    if (d[i]>0) u+=d[i]*d[i];
                                d = p.pos-high;
                                for (int i=0; i<3; ++i)
                                    if (d[i]>0) u+=d[i]*d[i];
                                return 0.5*k*u;
                            };
                        }
                    }

                    void to_json(json &j) const override {
                        if (type==cuboid)
                            j = {{"low", low}, {"high", high}};
                        if (type==sphere || type==cylinder)
                            j = {{"radius", radius}};
                        if (type==sphere) {
                            j["origo"] = origo;
                            j["scale"] = scale;
                        }
                        for (auto &i : m)
                            if (i.second==type)
                                j["type"] = i.first;
                        j["k"] = k/1.0_kJmol;
                        base::to_json(j);
                        _roundjson(j,5);
                    }
            }; //!< Confine particles to a sub-region of the simulation container

        /*
         * The keys of the `intra` map are group index and the values
         * is a vector of `BondData`. For bonds between groups, fill
         * in `inter` which is evaluated for every update of call to
         * `energy`.
         *
         * @todo Optimize.
         */
        template<typename Tspace>
            class Bonded : public Energybase {
                private:
                    Tspace& spc;
                    typedef typename Tspace::Tpvec Tpvec;
                    typedef std::vector<Potential::BondData> BondVector;
                    BondVector inter;  // inter-molecular bonds
                    std::map<int,BondVector> intra; // intra-molecular bonds

                    void update() {
                        intra.clear();
                        for (size_t i=0; i<spc.groups.size(); i++) {
                            if (!spc.groups.empty()) {
                                auto &g = spc.groups[i];
                                intra[i] = molecules<Tpvec>.at(g.id).bonds;
                                for (auto &b : intra[i])
                                    b.shift( std::distance(spc.p.begin(), g.begin()) );
                            }
                        }
                    } // finds and adds all intra-molecular bonds of active molecules 

                    double sum( const BondVector &v ) const {
                        double u=0;
                        for (auto &b : v)
                            u += b.energy(spc.p, spc.geo.distanceFunc);
                        return u;
                    } // sum energy in vector of BondData

                public:
                    Bonded(const json &j, Tspace &spc) : spc(spc) {
                        name = "bonded";
                        update();
                        if (j.is_object())
                            if (j.count("bondlist")==1)
                                inter = j["bondlist"].get<BondVector>();
                    }

                    void to_json(json &j) const override {
                        if (!inter.empty())
                            j["bondlist"] = inter;
                        if (!intra.empty()) {
                            json& _j = j["bondlist-intramolecular"];
                            _j = json::array();
                            for (auto &i : intra)
                                for (auto &b : i.second)
                                    _j.push_back(b);
                        }
                    }

                    double energy(Change &c) override {
                        double u=0;
                        if ( !c.empty() ) {
                            u = sum(inter); // energy of inter-molecular bonds
                            if ( c.all || c.dV ) {
                                for (auto& i : intra) // energy of intra-molecular bonds
                                    if (!spc.groups[i.first].empty()) // add only if group is active
                                        u += sum(i.second);
                            } else
                                for (auto &d : c.groups)
                                    if (d.internal) 
                                        u += sum( intra[d.index] );
                        }
                        return u;
                    }; // brute force -- refine this!
            };

        /**
         * @brief Nonbonded energy using a pair-potential
         */
        template<typename Tspace, typename Tpairpot>
            class Nonbonded : public Energybase {
                private:
                    double g2gcnt=0, g2gskip=0;
                protected:
                    typedef typename Tspace::Tgroup Tgroup;
                    double Rc2_g2g=pc::infty;

                    void to_json(json &j) const override {
                        j["pairpot"] = pairpot;
                        j["cutoff_g2g"] = std::sqrt(Rc2_g2g);
                    }

                    template<typename T>
                        inline bool cut(const T &g1, const T &g2) {
                            g2gcnt++;
                            if (g1.atomic || g2.atomic)
                                return false;
                            if ( spc.geo.sqdist(g1.cm, g2.cm)<Rc2_g2g )
                                return false;
                            g2gskip++;
                            return true;
                        } //!< true if group<->group interaction can be skipped

                    template<typename T>
                        inline double i2i(const T &a, const T &b) {
                            assert(&a!=&b && "a and b cannot be the same particle");
                            return pairpot(a, b, spc.geo.vdist(a.pos, b.pos));
                        }

                    /*
                     * Internal energy in group, calculating all with all or, if `index`
                     * is given, only a subset. Index specifies the internal index (starting
                     * at zero) of changed particles within the group.
                     */
                    double g_internal(const Tgroup &g, const std::vector<int> &index=std::vector<int>()) {
                        using namespace ranges;
                        double u=0;
                        if (index.empty()) // assume that all atoms have changed
                            for ( auto i = g.begin(); i != g.end(); ++i )
                                for ( auto j=i; ++j != g.end(); )
                                    u += i2i(*i, *j);
                        else { // only a subset have changed
                            auto fixed = view::ints( 0, int(g.size()) )
                                | view::remove_if(
                                        [&index](int i){return std::binary_search(index.begin(), index.end(), i);});
                            for (int i : index) // moved<->static
                                for (int j : fixed)
                                    u += i2i( *(g.begin()+i), *(g.begin()+j));
                            for (int i : index) // moved<->moved
                                for (int j : index)
                                    if (j>i)
                                        u += i2i( *(g.begin()+i), *(g.begin()+j));
                        }
                        return u;
                    }


                    /*
                     * Calculates the interaction energy of a particle, `i`,
                     * and checks (1) if it is already part of Space, or (2)
                     * external to space.
                     */
                    double i2all(const typename Tspace::Tparticle &i) {
                        double u=0;
                        auto it = spc.findGroupContaining(i); // iterator to group
                        if (it!=spc.groups.end()) {    // check if i belongs to group in space
                            for (auto &g : spc.groups) // i with all other particles
                                if (&g!=&(*it))        // avoid self-interaction
                                    if (!cut(g, *it))  // check g2g cut-off
                                        for (auto &j : g) // loop over particles in other group
                                            u += i2i(i,j);
                            for (auto &j : *it)        // i with all particles in own group
                                if (&j!=&i) 
                                    u += i2i(i,j);
                        } else // particle does not belong to any group
                            for (auto &g : spc.groups) // i with all other *active* particles
                                for (auto &j : g)      // (this will include only active particles)
                                    u += i2i(i,j);
                        return u;
                    }

                    /*
                     * Group-to-group energy. A subset of `g1` can be given with `index` which refers
                     * to the internal index (starting at zero) of the first group, `g1`.
                     */
                    virtual double g2g(const Tgroup &g1, const Tgroup &g2, const std::vector<int> &index=std::vector<int>()) {
                        double u = 0;
                        if (!cut(g1,g2)) {
                            if (index.empty()) // if index is empty, assume all in g1 have changed
                                for (auto &i : g1)
                                    for (auto &j : g2)
                                        u += i2i(i,j);
                            else // only a subset of g1
                                for (auto i : index)
                                    for (auto &j : g2)
                                        u += i2i( *(g1.begin()+i), j);
                        }
                        return u;
                    }

                public:
                    Tspace& spc;   //!< Space to operate on
                    Tpairpot pairpot; //!< Pair potential

                    Nonbonded(const json &j, Tspace &spc) : spc(spc) {
                        name="nonbonded";
                        pairpot = j;
                        Rc2_g2g = std::pow( j.value("cutoff_g2g", pc::infty), 2);
                    }

                    double energy(Change &change) override {
                        using namespace ranges;
                        double u=0;

                        if (!change.empty()) {

                            if (change.dV) {
#pragma omp parallel for reduction (+:u) schedule (dynamic) 
                                for ( auto i = spc.groups.begin(); i < spc.groups.end(); ++i ) {
                                    for ( auto j=i; ++j != spc.groups.end(); )
                                        u += g2g( *i, *j );
                                    if (i->atomic)
                                        u += g_internal(*i);
                                }
                                return u;
                            }

                            // did everything change?
                            if (change.all) {
#pragma omp parallel for reduction (+:u) schedule (dynamic) 
                                for ( auto i = spc.groups.begin(); i < spc.groups.end(); ++i ) {
                                    for ( auto j=i; ++j != spc.groups.end(); )
                                        u += g2g( *i, *j );
                                    u += g_internal(*i);
                                }
                                // more todo here...
                                return u;
                            }

                            // if exactly ONE molecule is changed
                            if (change.groups.size()==1) { 
                                auto& d = change.groups[0];
                                auto gindex = spc.groups.at(d.index).to_index(spc.p.begin()).first;
                                if (d.atoms.size()==1) // exactly one atom has moved
                                    return i2all(spc.p.at(gindex+d.atoms[0]));
                                auto& g1 = spc.groups.at(d.index);
                                for (auto &g2 : spc.groups)
                                    if (&g1 != &g2)
                                        u += g2g(g1, g2, d.atoms);
                                if (d.internal)
                                    u += g_internal(g1, d.atoms);
                                return u;
                            }

                            auto moved = change.touchedGroupIndex(); // index of moved groups
                            auto fixed = view::ints( 0, int(spc.groups.size()) )
                                | view::remove_if(
                                        [&moved](int i){return std::binary_search(moved.begin(), moved.end(), i);}
                                        ); // index of static groups

                            // moved<->moved
                            for ( auto i = moved.begin(); i != moved.end(); ++i )
                                for ( auto j=i; ++j != moved.end(); )
                                    u += g2g( spc.groups[*i], spc.groups[*j] );

                            // moved<->static
                            for ( auto i : moved)
                                for ( auto j : fixed)
                                    u += g2g(spc.groups[i], spc.groups[j]);

                            // more todo!
                        }
                        return u;
                    }

            }; //!< Nonbonded, pair-wise additive energy term

        template<typename Tspace, typename Tpairpot>
            class NonbondedCached : public Nonbonded<Tspace,Tpairpot> {
                private:
                    typedef Nonbonded<Tspace,Tpairpot> base;
                    typedef typename Tspace::Tgroup Tgroup;
                    Eigen::MatrixXf cache;

                    double g2g(const Tgroup &g1, const Tgroup &g2, const std::vector<int> &index=std::vector<int>()) override {
                        int i = &g1 - &base::spc.groups.front();
                        int j = &g2 - &base::spc.groups.front();
                        if (j<i)
                            std::swap(i,j);
                        if (base::key==Energybase::NEW) {        // if this is from the trial system,
                            double u = 0;
                            if (!base::cut(g1,g2)) {
                                for (auto &i : g1)
                                    for (auto &j : g2)
                                        u += base::i2i(i,j);
                            }
                            cache(i,j) = u;
                        }
                        return cache(i,j);                     // return (cached) value
                    }

                public:
                    NonbondedCached(const json &j, Tspace &spc) : base(j,spc) {
                        base::name += "EM";
                        cache.resize( spc.groups.size(), spc.groups.size() );
                        cache.setZero();
                        for ( auto i = base::spc.groups.begin(); i < base::spc.groups.end(); ++i ) {
                            for ( auto j=i; ++j != base::spc.groups.end(); ) {
                                int k = &(*i) - &base::spc.groups.front();
                                int l = &(*j) - &base::spc.groups.front();
                                if (l<k)
                                    std::swap(k,l);
                                double u = 0;
                                if (!base::cut(*i,*j)) {
                                    for (auto &k : *i)
                                        for (auto &l : *j)
                                            u += base::i2i(k,l);
                                }
                                cache(k,l) = u;
                            }
                        }
                    }

                    double energy(Change &change) override {
                        using namespace ranges;
                        double u=0;

                        if (!change.empty()) {

                            if (change.all || change.dV) {
#pragma omp parallel for reduction (+:u) schedule (dynamic) 
                                for ( auto i = base::spc.groups.begin(); i < base::spc.groups.end(); ++i ) {
                                    for ( auto j=i; ++j != base::spc.groups.end(); )
                                        u += g2g( *i, *j );
                                }
                                return u;
                            }

                            // if exactly ONE molecule is changed
                            if (change.groups.size()==1) { 
                                auto& d = change.groups[0];
                                auto& g1 = base::spc.groups.at(d.index);
                                for (auto &g2 : base::spc.groups) {
                                    if (&g1 != &g2)
                                        u += g2g(g1, g2, d.atoms);
                                }
                                return u;
                            }

                            auto moved = change.touchedGroupIndex(); // index of moved groups
                            auto fixed = view::ints( 0, int(base::spc.groups.size()) )
                                | view::remove_if(
                                        [&moved](int i){return std::binary_search(moved.begin(), moved.end(), i);}
                                        ); // index of static groups

                            // moved<->moved
                            for ( auto i = moved.begin(); i != moved.end(); ++i )
                                for ( auto j=i; ++j != moved.end(); ) {
                                    u += g2g( base::spc.groups[*i], base::spc.groups[*j] );
                            }
                            // moved<->static
                            for ( auto i : moved)
                                for ( auto j : fixed)
                                    u += g2g(base::spc.groups[i], base::spc.groups[j]);

                            // more todo!
                        }
                        return u;
                    }

                    void sync(Energybase *basePtr, Change &change) override {
                        auto other = dynamic_cast<decltype(this)>(basePtr);
                        assert(other);
                        if (change.all || change.dV)
                            cache.triangularView<Eigen::StrictlyUpper>() = (other->cache).template triangularView<Eigen::StrictlyUpper>();
                        else
                            for (auto &d : change.groups) {
                                for (size_t i=0; i<d.index; i++)
                                    cache(i,d.index) = other->cache(i,d.index);
                                for (size_t i=d.index+1; i<base::spc.groups.size(); i++)
                                    cache(d.index,i) = other->cache(d.index,i);
                            }
                    } //!< Copy energy matrix from other
            }; //!< Nonbonded with cached energies (Energy Matrix)

        /**
         * `udelta` is the total change of updating the energy function. If
         * not handled this will appear as an energy drift (which it is!). To
         * avoid this, this term is added to the energy but since it's the
         * same in both the trial and old state energies it will not affect
         * MC move acceptance.
         */
        template<typename Tspace>
            class Penalty : public Energybase {
                protected:
                    typedef typename Tspace::Tparticle Tparticle;
                    typedef typename Tspace::Tgroup Tgroup;
                    typedef typename Tspace::Tpvec Tpvec;
                    typedef typename std::shared_ptr<ReactionCoordinate::ReactionCoordinateBase> Tcoord;

                    Tspace &spc;
                    bool nodrift;
                    bool quiet;
                    size_t dim=0;
                    size_t cnt=0;       // number of calls to `sync()`
                    size_t nupdate;     // update frequency [steps]
                    size_t samplings;
                    size_t nconv=0;
                    double udelta=0;    // total energy change of updating penalty function
                    double scale;       // scaling factor for f0
                    double f0;          // penalty increment
                    std::string file, hisfile;
                    std::vector<Tcoord> rcvec; // vector of reaction coordinate functions
                    std::vector<double> coord; // latest reaction coordinate

                    Table<int> histo;
                    Table<double> penalty;
                public:
                    Penalty(const json &j, Tspace &spc) : spc(spc) {
                        using namespace ReactionCoordinate;
                        name = "penalty";
                        f0 = j.value("f0", 0.5);
                        scale = j.value("scale", 0.8);
                        quiet = j.value("quiet", true);
                        nupdate = j.value("update", 0);
                        samplings = j.value("samplings", 1);
                        nodrift = j.value("nodrift", true);
                        file = j.at("file").get<std::string>();
                        hisfile = j.value("histogram", "penalty-histogram.dat");
                        std::vector<double> binwidth, min, max;

                        if (scale<0 || scale>1)
                            throw std::runtime_error("`scale` must be in the interval [0:1]");

                        for (auto &i : j.at("coords"))
                            if (i.is_object())
                                if (i.size()==1) {
                                    std::shared_ptr<ReactionCoordinate::ReactionCoordinateBase> rc=nullptr;
                                    for (auto it=i.begin(); it!=i.end(); ++it) {
                                        if (it.key()=="atom")
                                            rc = std::make_shared<AtomProperty>(it.value(), spc);
                                        if (it.key()=="system")
                                            rc = std::make_shared<SystemProperty>(it.value(), spc);
                                        if (it.key()=="cmcm")
                                            rc = std::make_shared<MassCenterSeparation>(it.value(), spc);
                                        if (it.key()=="angle")
                                            rc = std::make_shared<PrincipalAxisAngle>(it.value(), spc);
                                        if (rc!=nullptr) {
                                            if (rc->min>=rc->max || rc->binwidth<=0)
                                                throw std::runtime_error("min<max and binwidth>0 required for '" + it.key() + "'");
                                            rcvec.push_back(rc);
                                            binwidth.push_back( rc->binwidth );
                                            min.push_back( rc->min );
                                            max.push_back( rc->max );
                                        } else
                                            throw std::runtime_error("unknown coordinate type '" + it.key() + "'");
                                    }
                                }
                        dim = binwidth.size();
                        if (dim<1 || dim>2)
                            throw std::runtime_error("minimum one maximum two coordinates required");

                        coord.resize(2,0);
                        histo.reInitializer(binwidth, min, max);
                        penalty.reInitializer(binwidth, min, max);

                        std::ifstream f(MPI::prefix+file);
                        if (f) {
                            cout << "Loading penalty function '" << MPI::prefix+file << "'" << endl;
                            std::string hash;
                            f >> hash >> f0 >> samplings;
                            cout << "f0 " << f0 << " samplings " << samplings << endl;
                            for (int row=0; row<penalty.rows(); row++)
                                for (int col=0; col<penalty.cols(); col++)
                                    if (!f.eof()) 
                                        f >> penalty(row,col);
                                    else
                                        throw std::runtime_error("penalty file dimension mismatch");
                            cout << "maxCoeff " << penalty.maxCoeff() << endl;
                        }
                     }

                    virtual ~Penalty() {
                        std::ofstream f1(MPI::prefix + file), f2(MPI::prefix + hisfile);
                        if (f1) f1 << "# " << f0 << " " << samplings << "\n" << penalty.array() - penalty.minCoeff() << endl;
                        if (f2) f2 << histo << endl;
                        cout << "nconv is " << nconv << endl;
                        // add function to save to numpy-friendly file...
                    }

                    void to_json(json &j) const override {
                        j["file"] = file;
                        j["scale"] = scale;
                        j["update"] = nupdate;
                        j["nodrift"] = nodrift;
                        j["histogram"] = hisfile;
                        j["f0_final"] = f0;
                        j["nconv"] = nconv;
                        auto& _j = j["coords"] = json::array();
                        for (auto rc : rcvec) {
                            json t;
                            t[rc->name] = *rc;
                            _j.push_back(t);
                        }
                     }

                    double energy(Change &change) override {
                        assert(rcvec.size()<=coord.size());
                        double u=0;
                        coord.resize( rcvec.size() );
                        if (!change.empty()) {
                            for (size_t i=0; i<rcvec.size(); i++) {
                                coord.at(i) = rcvec[i]->operator()();
                                if (!rcvec[i]->inRange(coord[i]))
                                    return pc::infty;
                            }
                            penalty.to_index(coord);
                            u = penalty[coord];
                        }
                        return (nodrift) ? u - udelta : u;
                    }

                    virtual void update(const std::vector<double> &c) {
                        if (++cnt % nupdate == 0 && f0>0) {
                            bool b = histo.minCoeff() >= samplings;
                            if (b && f0>0) {
                                double min = penalty.minCoeff();
                                penalty = penalty.array() - min;
                                if (!quiet)
                                    cout << "Barriers/kT. Penalty=" << penalty.maxCoeff()
                                        << " Histogram=" << std::log(double(histo.maxCoeff())/histo.minCoeff())
                                        << endl;
                                f0 = f0 * scale; // reduce penalty energy
                                samplings = std::ceil( samplings / scale );
                                histo.setZero();
                                udelta += -min;
                                nconv += 1;
                            }
                        }
                        coord = c;
                        histo[coord]++;
                        penalty[coord] += f0;
                        udelta += f0; 
                    }

                    void sync(Energybase *basePtr, Change &change) override {
                        auto other = dynamic_cast<decltype(this)>(basePtr);
                        assert(other);
                        update(other->coord);
                        other->update(other->coord);
                    } // @todo: this double the MPI communication
            };

#ifdef ENABLE_MPI
        template<typename Tspace, typename Base=Penalty<Tspace>>
            struct PenaltyMPI : public Base {
                using Base::samplings;
                using Base::penalty;
                using Base::udelta;
                using Base::scale;
                using Base::histo;
                using Base::coord;
                using Base::cnt;
                using Base::f0;
                using Base::file;
                using Base::hisfile;
                using Base::nconv;

                Eigen::VectorXi weights; // array w. mininum histogram counts
                Eigen::VectorXd buffer; // receive buffer for penalty functions

                PenaltyMPI(const json &j, Tspace &spc) : Base(j,spc) {
                    weights.resize( MPI::mpi.nproc() );
                    buffer.resize( penalty.size()*MPI::mpi.nproc() );
                }

                void update(const std::vector<double> &c) override {
                    using namespace Faunus::MPI;
                    double uold = penalty[c];
                    if (++cnt % this->nupdate == 0 && f0>0) {
                        int min = histo.minCoeff();
                        MPI_Barrier(mpi.comm);
                        MPI_Allgather(&min, 1, MPI_INT, weights.data(), 1, MPI_INT, mpi.comm);

                        if ( weights.maxCoeff() > samplings ) {
                            MPI_Gather(penalty.data(), penalty.size(), MPI_DOUBLE,
                                    buffer.data(), penalty.size(), MPI_DOUBLE, 0, mpi.comm);

                            if (mpi.isMaster()) {
                                penalty.setZero();
                                for (int i=0; i<mpi.nproc(); i++)
                                    penalty += Eigen::Map<Eigen::MatrixXd>(
                                            buffer.data()+i*penalty.size(), penalty.rows(), penalty.cols() ) 
                                            / double(mpi.nproc());
                                penalty = penalty.array() - penalty.minCoeff();
                            }

                            MPI_Bcast(penalty.data(), penalty.size(), MPI_DOUBLE, 0, mpi.comm);
                            nconv += 1;
                            std::ofstream f3(MPI::prefix + std::to_string(nconv) + file);
                            if (f3) f3 << "# " << f0 << " " << samplings << "\n" << penalty.array() << endl;
                            std::ofstream f4(MPI::prefix + std::to_string(nconv) + hisfile);
                            if (f4) f4 << histo << endl;
                            if (min>0 && !this->quiet)
                                cout << "Barriers/kT. Penalty=" << penalty.maxCoeff()
                                    << " Histogram=" << std::log(double(histo.maxCoeff())/histo.minCoeff()) << endl;

                            histo.setZero();
                            f0 = f0 * scale; // reduce penalty energy
                            samplings = std::ceil( samplings / scale );
                        }
                    }
                    coord = c;
                    histo[coord]++;
                    penalty[coord] += f0;
                    udelta += penalty[coord] - uold;
                } //!< Average penalty function across all nodes
    }; //!< Penalty function with MPI exchange
#endif

#ifdef FAU_POWERSASA
        template<class Tspace>
            class SASAEnergy : public Energybase {
                typedef typename Tspace::Tparticle Tparticle;
                typedef typename Tspace::Tpvec Tpvec;
                Tspace& spc;
                std::vector<float> sasa, radii; 
                std::vector<Point> coords;
                double probe; // sasa probe radius (angstrom)
                double conc=0;// co-solute concentration (mol/l)
                Average<double> avgArea; // average surface area
                std::shared_ptr<POWERSASA::PowerSasa<float,Point>> ps;

                void updateSASA(const Tpvec &p) {
                    radii.resize(p.size());
                    coords.resize(p.size());
                    std::transform(p.begin(), p.end(), coords.begin(), [](auto &a){ return a.pos;});
                    std::transform(p.begin(), p.end(), radii.begin(),
                            [this](auto &a){ return atoms<Tparticle>[a.id].sigma*0.5 + this->probe;});

                    ps->update_coords(coords, radii); // slowest step!

                    for (size_t i=0; i<p.size(); i++) {
                        auto &a = atoms<Tparticle>[p[i].id];
                        if (std::fabs(a.tfe)>1e-9 || std::fabs(a.tension)>1e-9)
                            ps->calc_sasa_single(i);
                    }
                    sasa = ps->getSasa();
                    assert(sasa.size()==p.size());
                }

                void to_json(json &j) const override {
                    using namespace u8;
                    j["molarity"] = conc / 1.0_molar;
                    j["radius"] = probe / 1.0_angstrom;
                    j[bracket("SASA")+"/"+angstrom+squared] = avgArea.avg() / 1.0_angstrom;
                    _roundjson(j,5);
                }

                public:
                SASAEnergy(const json &j, Tspace &spc) : spc(spc) {
                    name = "sasa";
                    cite = "doi:10.1002/jcc.21844";
                    probe = j.value("radius", 1.4) * 1.0_angstrom;
                    conc = j.value("molarity", conc) * 1.0_molar;

                    radii.resize(spc.p.size());
                    coords.resize(spc.p.size());
                    std::transform(spc.p.begin(), spc.p.end(), coords.begin(), [](auto &a){ return a.pos;});
                    std::transform(spc.p.begin(), spc.p.end(), radii.begin(),
                            [this](auto &a){ return atoms<Tparticle>[a.id].sigma*0.5 + this->probe;});

                    ps = std::make_shared<POWERSASA::PowerSasa<float,Point>>(coords,radii);
                }

                double energy(Change &change) override {
                    double u=0, A=0;
                    updateSASA(spc.p);
                    for (size_t i=0; i<sasa.size(); ++i) {
                        auto &a = atoms<Tparticle>[ spc.p[i].id ];
                        u += sasa[i] * (a.tension + conc * a.tfe);
                        A += sasa[i];
                    }
                    avgArea+=A; // sample average area for accepted confs. only
                    return u;
                }
            }; //!< SASA energy from transfer free energies
#endif

        struct Example2D : public Energybase {
            Point& i; // reference to 1st particle in the system
            template<typename Tspace>
                Example2D(const json &j, Tspace &spc): i(spc.p.at(0).pos) { name = "Example2D"; }
            double energy(Change &change) override {
                double s=1+std::sin(2*pc::pi*i.x())+std::cos(2*pc::pi*i.y());
                if (i.x()>=-2.00 && i.x()<=-1.25) return 1*s;
                if (i.x()>=-1.25 && i.x()<=-0.25) return 2*s;
                if (i.x()>=-0.25 && i.x()<= 0.75) return 3*s;
                if (i.x()>= 0.75 && i.x()<= 1.75) return 4*s;
                if (i.x()>= 1.75 && i.x()<= 2.00) return 5*s;
                return 1e10;
            }
        };

        template<typename Tspace>
            class Hamiltonian : public Energybase, public BasePointerVector<Energybase> {
                protected:
                    typedef typename Tspace::Tparticle Tparticle;
                    void to_json(json &j) const override {
                        for (auto i : this->vec)
                            j.push_back(*i);
                    }

                    void addEwald(const json &j, Tspace &spc) {
                        if (j.count("coulomb")==1)
                            if (j["coulomb"].at("type")=="ewald")
                                push_back<Energy::Ewald<Tspace>>(j["coulomb"], spc);
                    } //!< Adds an instance of reciprocal space Ewald energies (if appropriate)

                public:
                    Hamiltonian(Tspace &spc, const json &j) {
                        using namespace Potential;

                        typedef CombinedPairPotential<CoulombGalore,LennardJones<Tparticle>> CoulombLJ;
                        typedef CombinedPairPotential<CoulombGalore,HardSphere<Tparticle>> CoulombHS;
                        typedef CombinedPairPotential<CoulombGalore,WeeksChandlerAndersen<Tparticle>> CoulombWCA;
                        typedef CombinedPairPotential<Coulomb,WeeksChandlerAndersen<Tparticle>> PrimitiveModelWCA;

                        Energybase::name="hamiltonian";
                        for (auto &m : j.at("energy")) {// loop over move list
                            size_t oldsize = vec.size();
                            for (auto it=m.begin(); it!=m.end(); ++it) {
                                try {
                                    if (it.key()=="nonbonded_coulomblj")
                                        push_back<Energy::Nonbonded<Tspace,CoulombLJ>>(it.value(), spc);

                                    if (it.key()=="nonbonded")
                                        push_back<Energy::Nonbonded<Tspace,FunctorPotential<typename Tspace::Tparticle>>>(it.value(), spc);

                                    if (it.key()=="nonbonded_coulombhs")
                                        push_back<Energy::Nonbonded<Tspace,CoulombHS>>(it.value(), spc);

                                    if (it.key()=="nonbonded_coulombwca")
                                        push_back<Energy::Nonbonded<Tspace,CoulombWCA>>(it.value(), spc);

                                    if (it.key()=="nonbonded_pmwca")
                                        push_back<Energy::Nonbonded<Tspace,PrimitiveModelWCA>>(it.value(), spc);

                                    if (it.key()=="nonbonded_deserno")
                                        push_back<Energy::NonbondedCached<Tspace,DesernoMembrane<typename Tspace::Tparticle>>>(it.value(), spc);

                                    if (it.key()=="nonbonded_desernoAA")
                                        push_back<Energy::NonbondedCached<Tspace,DesernoMembraneAA<typename Tspace::Tparticle>>>(it.value(), spc);

                                    if (it.key()=="bonded")
                                        push_back<Energy::Bonded<Tspace>>(it.value(), spc);

                                    if (it.key()=="confine")
                                        push_back<Energy::Confine<Tspace>>(it.value(), spc);

                                    if (it.key()=="example2d")
                                        push_back<Energy::Example2D>(it.value(), spc);

                                    if (it.key()=="isobaric")
                                        push_back<Energy::Isobaric<Tspace>>(it.value(), spc);

                                    if (it.key()=="penalty")
#ifdef ENABLE_MPI
                                        push_back<Energy::PenaltyMPI<Tspace>>(it.value(), spc);
#else
                                        push_back<Energy::Penalty<Tspace>>(it.value(), spc);
#endif
#ifdef ENABLE_POWERSASA
                                    if (it.key()=="sasa")
                                        push_back<Energy::SASAEnergy<Tspace>>(it.value(), spc);
#endif
                                    // additional energies go here...

                                    addEwald(it.value(), spc); // add reciprocal Ewald terms if appropriate

                                    if (vec.size()==oldsize)
                                        std::cerr << "warning: ignoring unknown energy '" << it.key() << "'" << endl;

                                } catch (std::exception &e) {
                                    throw std::runtime_error("Error adding energy '" + it.key() + "': " + e.what());
                                }
                            }
                        }
                    }

                    double energy(Change &change) override {
                        double du=0;
                        for (auto i : this->vec) {
                            i->key=key;
                            du += i->energy(change);
                        }
                        return du;
                    } //!< Energy due to changes

                    void sync(Hamiltonian &other, Change &change) {
                        assert(other.size()==size());
                        for (size_t i=0; i<size(); i++)
                            this->vec[i]->sync( other.vec[i].get(), change);
                    }

            }; //!< Aggregates and sum energy terms

    }//namespace
}//namespace
