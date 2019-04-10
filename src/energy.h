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

#ifdef ENABLE_POWERSASA
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
                TimeRelativeOfTotal<std::chrono::microseconds> timer;
                virtual double energy(Change&)=0; //!< energy due to change
                virtual void to_json(json &j) const; //!< json output
                virtual void sync(Energybase*, Change&);
                virtual void init(); //!< reset and initialize
                virtual inline void force(std::vector<Point>&) {}; // update forces on all particles
                inline virtual ~Energybase() {};
        };

        void to_json(json &j, const Energybase &base); //!< Converts any energy class to json object

        /**
         * @brief Check for overlap between atoms and the simulation container
         *
         * If found infinite energy is returned. Not needed for cuboidal geometry
         * as there's nover any overlap due to PBC.
         */
        template<class Tspace>
            struct ContainerOverlap : public Energybase {
                const Tspace &spc;
                ContainerOverlap(const Tspace &spc) : spc(spc) {
                    name = "ContainerOverlap";
                }
                double energy(Change &change) override {
                    //if (spc.geo.type not_eq Geometry::Chameleon::CUBOID) // cuboid have PBC in all directions
                        if (change) {
                            // all groups have been updated
                            if (change.dV or change.all) {
                                for (auto &g : spc.groups) // loop over *all* groups in system
                                    for (auto &p : g) // loop over *all* active particles in group
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
                                }
                                else
                                    // only a subset of atoms were updated
                                    for (int i : d.atoms) // loop over specific atoms
                                        if (spc.geo.collision( (g.begin()+i)->pos ))
                                            return pc::infty;
                            }
                        }
                    return 0;
                }
            };

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

            void update(const Point &box);
        };

        void from_json(const json &j, EwaldData &d);

        void to_json(json &j, const EwaldData &d);

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
                        for (int k=0; k<d.Qion.size(); k++)
                            E += d.Aks[k] * std::norm( d.Qion[k] );
                    return 2 * pc::pi / spc->geo.getVolume() * E * d.lB;
                }
            };

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] Ewald - IonIonPolicy")
        {
            using doctest::Approx;
            typedef Space<Geometry::Chameleon, Particle<Charge,Dipole>> Tspace;

            Tspace spc;
            spc.p.resize(2);
            spc.geo  = R"( {"type": "cuboid", "length": 10} )"_json;
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
                    Tspace& spc;
                public:

                    Ewald(const json &j, Tspace &spc) : policy(spc), spc(spc) {
                        name = "ewald";
                        data = j;
                        init();
                    }

                    void init() override {
                        data.update( spc.geo.getLength() );
                        policy.updateComplex(data); // brute force. todo: be selective
                    }

                    double energy(Change &change) override {
                        double u=0;
                        if (change) {
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

                    void sync(Energybase *basePtr, Change&) override {
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
                        if (change.dV || change.all || change.dN) {
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

        /**
         * @brief Constrain system using reaction coordinates
         *
         * If outside specified `range`, infinity energy is returned, causing rejection.
         */
        class Constrain : public Energybase {
            private:
                std::string type;
                std::shared_ptr<ReactionCoordinate::ReactionCoordinateBase> rc=nullptr;
            public:
                template<class Tspace>
                    Constrain(const json &j, Tspace &spc) {
                        using namespace Faunus::ReactionCoordinate;
                        name = "constrain";
                        type = j.at("type").get<std::string>();
                        try {
                            if      (type=="atom")     rc = std::make_shared<AtomProperty>(j, spc);
                            else if (type=="molecule") rc = std::make_shared<MoleculeProperty>(j, spc);
                            else if (type=="system")   rc = std::make_shared<SystemProperty>(j, spc);
                            else if (type=="cmcm")     rc = std::make_shared<MassCenterSeparation>(j, spc);
                            if (rc==nullptr)
                                throw std::runtime_error("unknown coordinate type");

                        } catch (std::exception &e) {
                            throw std::runtime_error("error for reaction coordinate '"
                                    + type + "': " + e.what() + usageTip["coords=["+type+"]"]  );
                        }
                    }

                inline double energy(Change &change) override {
                    if (change) {
                        double val = (*rc)(); // calculate reaction coordinate
                        if (not rc->inRange(val)) // is it within allowed range?
                            return pc::infty;  // if not, return infinite energy
                    }
                    return 0;
                }

                inline void to_json(json &j) const override {
                    j = *rc;
                    j["type"] = type;
                    j.erase("resolution");
                }
        };

        /**
         * @brief Base class for external potentials
         *
         * This will apply an external energy to a defined
         * list of molecules, either acting on individual
         * atoms or the mass-center. The specific energy
         * function, `func` is injected in derived classes.
         */
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
                            if (molids.find(g.id) != molids.end()) {
                                if (COM and g.atomic==false) { // apply only to center of mass
                                    Tparticle cm; // fake particle representin molecule
                                    cm.charge = Geometry::monopoleMoment(g.begin(), g.end());
                                    cm.pos = g.cm;
                                    u = func(cm);
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

                    /*
                     * @todo The `dN` check is very inefficient
                     * as it calculates the external potential on *all*
                     * particles.
                     */
                    double energy(Change &change) override {
                        assert(func!=nullptr);
                        double u=0;
                        if (change.dV or change.all or change.dN) {
                            for (auto &g : spc.groups) { // check all groups
                                u += _energy(g);
                                if (std::isnan(u))
                                    break;
                            }
                        } else
                            for (auto &d : change.groups) {
                                auto &g = spc.groups.at(d.index); // check specified groups
                                if (d.all or COM)  // check all atoms in group
                                    u += _energy(g); // _energy also checks for molecule id
                                else {       // check only specified atoms in group
                                    if (molids.find(g.id) != molids.end())
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

        /**
         * @brief Custom external potential on molecules
         */
        template<typename Tspace, typename base=ExternalPotential<Tspace>>
            class CustomExternal : public base { 
                private:
                    ExprFunction<double> expr;
                    struct Data { // variables
                        double q=0, x=0, y=0, z=0;
                    };
                    Data d;
                    json jin; // initial json input

                 public:
                    CustomExternal(const json &j, Tspace &spc) : base(j,spc) {
                        base::name = "customexternal";
                        jin = j;
                        auto &_j = jin["constants"];
                        if (_j==nullptr)
                            _j = json::object();
                        _j["e0"] = pc::e0;
                        _j["kB"] = pc::kB;
                        _j["kT"] = pc::kT();
                        _j["Nav"] = pc::Nav;
                        _j["T"] = pc::temperature;
                        expr.set(jin, {
                                {"q",&d.q}, {"x",&d.x}, {"y",&d.y}, {"z",&d.z} } );
                        base::func = [&](const typename Tspace::Tparticle &a) {
                            d.x = a.pos.x();
                            d.y = a.pos.y();
                            d.z = a.pos.z();
                            d.q = a.charge;
                            return expr();
                        };
                    }

                    void to_json(json &j) const override {
                        j = jin;
                        base::to_json(j);
                    }
            };

        /**
         * @brief Mean field electric potential from outside rectangular simulation box.
         * @date Asljunga, December 2010.
         *
         * Calculates the average potential outside a simulation box due to ion
         * densities inside the box, extended to infinity.
         * The update() function calculates the average charge densities in slits in the xy-plane
         * and with a spacing dz. This is used to evaluate the electric potential along
         * the z-axis by using the method by Torbjorn and CO (Mol. Phys. 1996, 87:407). To avoid
         * energy drifts, update() returns the energy change brought about by updating the charge profile.
         */
        template<typename Tspace, typename base=ExternalPotential<Tspace>>
            class ExternalAkesson : public base { 
                private:
                    using base::spc;
                    std::string filename;                //!< File name for average charge
                    bool fixed;
                    unsigned int nstep=0;                //!< Internal between samples
                    unsigned int nphi=0;                 //!< Distance between phi updating
                    unsigned int updatecnt=0;            //!< Number of time rho has been updated
                    double epsr;                         //!< Relative dielectric constant
                    double dz;                           //!< z spacing between slits (A)
                    double lB;                           //!< Bjerrum length (A)
                    double halfz;                        //!< Half box length in z direction
                    Equidistant2DTable<double> Q;        //!< instantaneous net charge

                public:
                    unsigned int cnt=0;                  //!< Number of charge density updates
                    Equidistant2DTable<double,Average<double>> rho; //!< Charge density at z (unit A^-2)
                    Equidistant2DTable<double> phi;         //!< External potential at z (unit: beta*e)

                private:
                    void to_json(json &j) const override {
                        j = {
                            {"lB",lB}, {"dz",dz}, {"nphi",nphi},
                            {"epsr",epsr}, {"file",filename}, {"nstep",nstep},
                            {"Nupdates",updatecnt}, {"fixed",fixed}
                        };
                        base::to_json(j);
                        _roundjson(j,5);
                    }

                    void save() {
                        std::ofstream f(filename);
                        if (f) {
                            f.precision(16);
                            f << rho;
                        }
                        else throw std::runtime_error("cannot save file '"s + filename + "'");
                    }

                    void load() {
                        std::ifstream f(filename);
                        if (f) {
                            rho << f;
                            update_phi();
                        }
                        else std::cerr << "density file '" << filename << "' not loaded." << endl;
                     }

                    /*
                     * This is Eq. 15 of the mol. phys. 1996 paper by Greberg et al.
                     * (sign typo in manuscript: phi^infty(z) should be "-2*pi*z" on page 413, middle)
                     */
                    inline double phi_ext(double z, double a) const {
                        double a2 = a*a,
                               z2 = z*z;
                        return
                            - 2 * pc::pi * z
                            - 8 * a * std::log( ( std::sqrt( 2*a2 + z2 ) + a ) / std::sqrt( a2+z2 ) )
                            + 2 * z * (0.5*pc::pi
                                    +std::asin( ( a2*a2 - z2*z2 - 2*a2*z2) / std::pow(a2+z2,2) ) );
                    }

                    void sync(Energybase *basePtr, Change&) override {
                        if (not fixed) {
                            auto other = dynamic_cast<decltype(this)>(basePtr);
                            assert(other);
                            // only trial energy (new) require sync
                            if (other->key==Energybase::OLD)
                                if (cnt != other->cnt) {
                                    assert(cnt < other->cnt && "trial cnt's must be smaller");
                                    cnt = other->cnt;
                                    rho = other->rho;
                                    phi = other->phi;
                                }
                        }
                    }

                    // update average charge density
                    void update_rho() {
                        updatecnt++;
                        Point L = spc.geo.getLength();
                        double area = L.x()*L.y();
                        if (L.x() not_eq L.y() or 0.5*L.z()!=halfz)
                            throw std::runtime_error("Requires box Lx=Ly and Lz=const.");

                        Q.clear();
                        for (auto &g : spc.groups) // loop over all groups
                            for (auto &p : g)        // ...and their active particles
                                Q( p.pos.z() ) += p.charge;
                        for (double z=-halfz; z<=halfz; z+=dz)
                            rho(z) += Q(z) / area;
                    }

                    // update average external potential
                    void update_phi() {
                        Point L = spc.geo.getLength();
                        double a = 0.5*L.x();
                        for (double z=-halfz; z<=halfz; z+=dz) {
                            double s=0;
                            for (double zn=-halfz; zn<=halfz; zn+=dz)
                                if (rho(zn).cnt>0)
                                    s += rho(zn).avg() * phi_ext( std::fabs(z-zn), a );  // Eq. 14 in Greberg paper
                            phi(z) = lB*s;
                        }
                    }

                public:
                    ExternalAkesson(const json &j, Tspace &spc) : base(j,spc) {
                        base::name = "akesson";
                        base::cite = "doi:10/dhb9mj";

                        xjson _j = j; // json variant where items are deleted after access
                        _j.erase("com");
                        _j.erase("molecules");

                        nstep = _j.at("nstep").get<unsigned int>();
                        epsr = _j.at("epsr").get<double>();
                        fixed = _j.value("fixed", false);
                        nphi = _j.value("nphi", 10);

                        halfz = 0.5*spc.geo.getLength().z();
                        lB = pc::lB(epsr);

                        dz = _j.value("dz", 0.2); // read z resolution
                        Q.setResolution(dz, -halfz, halfz);
                        rho.setResolution(dz, -halfz, halfz);
                        phi.setResolution(dz, -halfz, halfz);

                        filename = _j.value("file", "mfcorr.dat"s);
                        load();

                        base::func = [&phi=phi](const typename Tspace::Tparticle &p) {
                            return p.charge * phi( p.pos.z() );
                        };

                        if (not _j.empty()) // throw exception of unused/unknown keys are passed
                            throw std::runtime_error("unused key(s) for '"s + base::name + "':\n" + _j.dump());
                    }

                    double energy(Change &change) override {
                        if (not fixed) // pho(z) unconverged, keep sampling
                            if (base::key==Energybase::OLD) { // only sample on accepted configs
                                cnt++;
                                if (cnt % nstep == 0)
                                    update_rho();
                                if (cnt % nstep*nphi == 0)
                                    update_phi();
                            }
                        return base::energy(change);
                    }

                    ~ExternalAkesson() {
                        // save only if still updating and if energy type is "OLD",
                        // that is, accepted configurations (not trial)
                        if (not fixed and base::key==Energybase::OLD)
                            save();
                    }
            };


        /**
         * @brief Confines molecules inside geometric shapes
         */
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
                        {"sphere", sphere},
                        {"cylinder", cylinder},
                        {"cuboid", cuboid}
                    };

                public:
                    Confine(const json &j, Tspace &spc) : base(j,spc) {
                        base::name = "confine";
                        k = value_inf(j, "k") * 1.0_kJmol; // get floating point; allow inf/-inf
                        type = m.at( j.at("type") );

                        if (type==sphere or type==cylinder) {
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
                                spc.scaleVolumeTriggers.push_back( [&radius=radius](Tspace&, double Vold, double Vnew) {
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
                        if (type==sphere or type==cylinder)
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
            typedef std::vector<std::shared_ptr<Potential::BondData>> BondVector;
            BondVector inter;  // inter-molecular bonds
            std::map<int, BondVector> intra; // intra-molecular bonds

          private:
            void update_intra() {
                using namespace Potential;
                intra.clear();
                for (size_t i = 0; i < spc.groups.size(); i++) {
                    auto& group = spc.groups[i];
                    for (auto& bond : molecules<Tpvec>.at(group.id).bonds) {
                        intra[i].push_back( bond->clone() ); // deep copy BondData from MoleculeData
                        intra[i].back()->shift( std::distance(spc.p.begin(), group.begin()) );
                        Potential::setBondEnergyFunction( intra[i].back(), spc.p );
                    }
                }
            } // finds and adds all intra-molecular bonds of active molecules

            double sum_energy(const BondVector &bonds) const {
                double energy = 0;
                for (auto& bond : bonds) {
                    assert(bond->hasEnergyFunction());
                    energy += bond->energy(spc.geo.getDistanceFunc());
                }
                return energy;
            } // sum energy in vector of BondData

            double sum_energy(const BondVector &bonds, int particle_ndx) const {
                double energy = 0;
                for (auto& bond : bonds) {
                    if (std::find(bond->index.begin(), bond->index.end(), particle_ndx) != bond->index.end()) {
                        assert(bond->hasEnergyFunction());
                        energy += bond->energy(spc.geo.getDistanceFunc());
                    }
                }
                return energy;
            } // sum energy in vector of BondData for matching particle index

          public:
            Bonded(const json &j, Tspace &spc) : spc(spc) {
                name = "bonded";
                update_intra();
                if (j.is_object())
                    if (j.count("bondlist")==1)
                        inter = j["bondlist"].get<BondVector>();
                for (auto& i : inter) // set all energy functions
                    Potential::setBondEnergyFunction( i, spc.p );
            }

            void to_json(json &j) const override {
                if (!inter.empty())
                    j["bondlist"] = inter;
                if (!intra.empty()) {
                    json& _j = j["bondlist-intramolecular"];
                    _j = json::array();
                    for (auto& i : intra)
                        for (auto& b : i.second)
                            _j.push_back(b);
                }
            }

            double energy(Change &change) override {
                double energy = 0;
                if (change) {
                    energy += sum_energy(inter); // energy of inter-molecular bonds

                    if (change.all || change.dV) { // compute all active groups
                        for (auto& i : intra) { // energies of intra-molecular bonds
                            if (! spc.groups[i.first].empty()) {// add only if group is active
                                energy += sum_energy(i.second);
                            }
                        }
                    } else { // compute only the affected groups
                        for (auto& group : change.groups) {
                            auto& intra_group = intra[group.index];
                            if (group.internal) {
                                if (group.all) { // all internal positions updated
                                    energy += sum_energy(intra_group);
                                } else { // only partial update
                                    // offset = index of first particle in group
                                    int offset = std::distance(spc.p.begin(), spc.groups[group.index].begin());
                                    for (int i : group.atoms) // d.atoms is relative to group
                                        energy += sum_energy(intra_group, i + offset);
                                }
                            }
                        } // for-loop over groups
                    }
                }
                return energy;
            }; // brute force -- refine this!
        };

        /**
         * @brief Nonbonded energy using a pair-potential
         */
        template<typename Tspace, typename Tpairpot>
            class Nonbonded : public Energybase {
                private:
                    double g2gcnt=0, g2gskip=0;
                    PairMatrix<double> cutoff2; // matrix w. group-to-group cutoff

                protected:
                    typedef typename Tspace::Tpvec Tpvec;
                    typedef typename Tspace::Tgroup Tgroup;
                    double Rc2_g2g=pc::infty;

                    // control of when OpenMP should be used
                    bool omp_enable=false;
                    bool omp_i2all=false;
                    bool omp_g2g=false;
                    bool omp_p2p=false;

                    void to_json(json &j) const override {
                        j["pairpot"] = pairpot;
                        if (omp_enable) {
                            json _a = json::array();
                            if (omp_p2p) _a.push_back("p2p");
                            if (omp_g2g) _a.push_back("g2g");
                            if (omp_i2all) _a.push_back("i2all");
                            j["openmp"] = _a;
                        }
                        j["cutoff_g2g"] = json::object();
                        auto &_j = j["cutoff_g2g"];
                        for (auto &a : Faunus::molecules<typename Tspace::Tpvec>)
                            for (auto &b : Faunus::molecules<typename Tspace::Tpvec>)
                                if (a.id()>=b.id())
                                    _j[a.name+" "+b.name] = sqrt( cutoff2(a.id(), b.id()) );
                    }

                    template<typename T>
                        inline bool cut(const T &g1, const T &g2) {
                            g2gcnt++;
                            if (g1.atomic || g2.atomic)
                                return false;
                            if ( spc.geo.sqdist(g1.cm, g2.cm) < cutoff2(g1.id, g2.id) )
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
                     * from zero) of changed particles within the group.
                     */
                    double g_internal(const Tgroup &g, const std::vector<int> &index=std::vector<int>()) {
                        using namespace ranges;
                        double u = 0;
                        if (index.empty() and not molecules<Tpvec>.at(g.id).rigid) // assume that all atoms have changed
                            for ( auto i = g.begin(); i != g.end(); ++i )
                                for ( auto j=i; ++j != g.end(); )
                                    u += i2i(*i, *j);
                        else { // only a subset has changed
                            auto fixed = view::ints( 0, int(g.size()) )
                                | view::remove_if(
                                        [&index](int i){return std::binary_search(index.begin(), index.end(), i);});
                            for (int i : index) { // moved<->static
                                for (int j : fixed ) {
                                    u += i2i( *(g.begin()+i), *(g.begin()+j));
                                }
                            }
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
#pragma omp parallel for reduction (+:u) if (omp_enable and omp_i2all)
                            for (size_t ig=0; ig<spc.groups.size(); ig++) {
                                auto &g = spc.groups[ig];
                                if (&g!=&(*it))        // avoid self-interaction
                                    if (not cut(g, *it))  // check g2g cut-off
                                        for (auto &j : g) // loop over particles in other group
                                            u += i2i(i,j);
                            }
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
                     * to the internal index (starting at zero) of the first group, `g1
                     * NOTE: the interpretation of this function is extended to also consider the mutual interactions
                     * of a subset of each group and in such case returns sub1 <-> 2 and !sub1<->sub2,
                     * hence excluding !sub1 <-> !sub2 in comparision to calling onconstrained g2g. In absence
                     * of sub1 any sub2 is ignored.
                     */
                    virtual double g2g(
                            const Tgroup &g1,
                            const Tgroup &g2,
                            const std::vector<int> &index=std::vector<int>(),
                            const std::vector<int> &jndex=std::vector<int>() )
                    {
                        using namespace ranges;
                        double u = 0;
                        if (not cut(g1,g2)) {
                            if ( index.empty() && jndex.empty() ) // if index is empty, assume all in g1 have changed
#pragma omp parallel for reduction (+:u) schedule (dynamic) if (omp_enable and omp_p2p)  
                                for (size_t i=0; i<g1.size(); i++)
                                    for (size_t j=0; j<g2.size(); j++)
                                        u += i2i( *(g1.begin()+i), *(g2.begin()+j) );
                            else { // only a subset of g1
                                for (auto i : index)
                                    for (auto j=g2.begin(); j!=g2.end(); ++j)
                                        u += i2i( *(g1.begin()+i), *j);
                                if ( not jndex.empty() ) {
                                    auto fixed = view::ints( 0, int(g1.size()) )
                                        | view::remove_if(
                                                [&index](int i){return std::binary_search(index.begin(), index.end(), i);});
                                    for (auto i : jndex) // moved2        <-|
                                        for (auto j : fixed) // static1   <-|
                                            u += i2i( *(g2.begin()+i), *(g1.begin()+j));
                                }
                            }
                        }
                        return u;
                    }

                public:
                    Tspace& spc;   //!< Space to operate on
                    Tpairpot pairpot; //!< Pair potential

                    Nonbonded(const json &j, Tspace &spc) : spc(spc) {
                        name="nonbonded";
                        pairpot = j;

                        // controls for OpenMP
                        auto it = j.find("openmp");
                        if (it != j.end())
                            if (it->is_array())
                                if (it->size()>0) {
                                    omp_enable=true;
                                    for (const std::string &k : *it)
                                        if (k=="g2g") omp_g2g=true;
                                        else if (k=="p2p") omp_p2p=true;
                                        else if (k=="i2all") omp_i2all=true;
#ifndef _OPENMP
                                    std::cerr << "warning: nonbonded requests unavailable OpenMP." << endl;
#endif
                                }

                        // disable all group-to-group cutoffs by setting infinity
                        for (auto &i : Faunus::molecules<typename Tspace::Tpvec>)
                            for (auto &j : Faunus::molecules<typename Tspace::Tpvec>)
                                cutoff2.set(i.id(), j.id(), pc::infty);

                        it = j.find("cutoff_g2g");
                        if (it != j.end()) {
                            // old style input w. only a single cutoff
                            if (it->is_number()) {
                                Rc2_g2g = std::pow( it->get<double>(), 2 );
                                for (auto &i : Faunus::molecules<typename Tspace::Tpvec>)
                                    for (auto &j : Faunus::molecules<typename Tspace::Tpvec>)
                                        cutoff2.set(i.id(), j.id(), Rc2_g2g);
                            }
                            // new style input w. multiple cutoffs between molecules
                            else if (it->is_object()) {
                                // ensure that there is a default, fallback cutoff
                                Rc2_g2g = std::pow( it->at("default").get<double>(), 2);
                                for (auto &i : Faunus::molecules<typename Tspace::Tpvec>)
                                    for (auto &j : Faunus::molecules<typename Tspace::Tpvec>)
                                        cutoff2.set(i.id(), j.id(), Rc2_g2g);
                                // loop for space separated molecule pairs in keys
                                for (auto& i : it->items()) {
                                    auto v = words2vec<std::string>( i.key() );
                                    if (v.size()==2) {
                                        int id1 = (*findName( Faunus::molecules<typename Tspace::Tpvec>, v[0])).id();
                                        int id2 = (*findName( Faunus::molecules<typename Tspace::Tpvec>, v[1])).id();
                                        cutoff2.set( id1, id2, std::pow(i.value().get<double>(),2) );
                                    }
                                }
                            }
                        }
                    }

                    void force(std::vector<Point> &forces) override {
                        auto &p = spc.p; // alias to particle vector (reference)
                        assert(forces.size() == p.size() && "the forces size must match the particle size");
                        for (size_t i=0; i<p.size()-1; i++)
                            for (size_t j=i+1; j<p.size(); j++) {
                                //Point r = spc.geo.vdist(p[i].pos, p[j].pos); // minimum distance vector 
                                Point f ;//= pairpot.force( p[i], p[j], r.squaredNorm(), r );
                                forces[i] += f;
                                forces[j] -= f;
                            }
                    }

                    double energy(Change &change) override {
                        using namespace ranges;
                        double u=0;

                        if (change) {

                            if (change.dV) {
#pragma omp parallel for reduction (+:u) schedule (dynamic) if (omp_enable and omp_g2g)  
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
#pragma omp parallel for reduction (+:u) schedule (dynamic) if (omp_enable and omp_g2g)
                                for ( auto i = spc.groups.begin(); i < spc.groups.end(); ++i ) {
                                    for ( auto j=i; ++j != spc.groups.end(); )
                                        u += g2g( *i, *j );
                                    u += g_internal(*i);
                                }
                                // more todo here...
                                return u;
                            }

                            // if exactly ONE molecule is changed
                            if (change.groups.size()==1 && not change.dN) {
                                auto& d = change.groups[0];
                                auto gindex = spc.groups.at(d.index).to_index(spc.p.begin()).first;

                                // exactly one atom has move
                                if (d.atoms.size()==1)
                                    return i2all(spc.p.at(gindex+d.atoms[0]));

                                // more atoms moved
                                auto& g1 = spc.groups.at(d.index);
                                //for (auto &g2 : spc.groups)
#pragma omp parallel for reduction (+:u) schedule (dynamic) if (omp_enable and omp_g2g)
                                for (size_t i=0; i<spc.groups.size(); i++) {
                                    auto &g2 = spc.groups[i];
                                    if (&g1 != &g2)
                                        u += g2g(g1, g2, d.atoms);
                                }
                                if (d.internal)
                                    u += g_internal(g1, d.atoms);
                                return u;
                            }

                            auto moved = change.touchedGroupIndex(); // index of moved groups
                            auto fixed = view::ints( 0, int(spc.groups.size()) )
                                | view::remove_if(
                                        [&moved](int i){return std::binary_search(moved.begin(), moved.end(), i);}
                                        ); // index of static groups

                            if (change.dN) {
                                /*auto moved = change.touchedGroupIndex(); // index of moved groups
                                std::vector<int> Moved;
                                for (auto i: moved) {
                                    Moved.push_back(i);
                                }
                                std::sort( Moved.begin(), Moved.end() );
                                auto fixed = view::ints( 0, int(spc.groups.size()) )
                                    | view::remove_if(
                                            [&Moved](int i){return std::binary_search(Moved.begin(), Moved.end(), i);}
                                            ); // index of static groups*/
                                for ( auto cg1 = change.groups.begin(); cg1 < change.groups.end() ; ++cg1 ) { // Loop over all changed groups
                                    std::vector<int> ifiltered, jfiltered; // Active atoms
                                    for (auto i: cg1->atoms) {
                                        if ( i < spc.groups.at(cg1->index).size() )
                                            ifiltered.push_back(i);
                                    }
                                    // Skip if the group is empty
                                    if ( not ifiltered.empty() )
                                        for ( auto j : fixed )
                                            u += g2g( spc.groups.at(cg1->index), spc.groups[j], ifiltered, jfiltered );

                                    for ( auto cg2 = cg1; ++cg2 != change.groups.end(); ) {
                                        for (auto i: cg2->atoms)
                                            if ( i < spc.groups.at(cg2->index).size() )
                                                jfiltered.push_back(i);
                                        // Skip if both groups are empty
                                        if ( not (ifiltered.empty() && jfiltered.empty()) )
                                            u += g2g( spc.groups.at(cg1->index),  spc.groups.at(cg2->index), ifiltered, jfiltered );
                                        jfiltered.clear();
                                    }
                                    if ( not ifiltered.empty() && cg1->dNatomic )
                                        u += g_internal( spc.groups.at( cg1->index ), ifiltered );
                                }
                                return u;
                            }

                            // moved<->moved
                            if (change.moved2moved) {
                                for ( auto i = moved.begin(); i != moved.end(); ++i )
                                    for ( auto j=i; ++j != moved.end(); )
                                        u += g2g( spc.groups[*i], spc.groups[*j] );
                            }

                            // moved<->static
                            if (omp_enable and omp_g2g) {
                                std::vector<std::pair<int,int>> pairs( size(moved) * size(fixed) );
                                size_t cnt=0;
                                for (auto i : moved)
                                    for( auto j : fixed)
                                        pairs[cnt++] = {i,j};
#pragma omp parallel for reduction (+:u) schedule (dynamic) if (omp_enable and omp_g2g)
                                for (size_t i=0; i<pairs.size(); i++)
                                    u += g2g(spc.groups[pairs[i].first], spc.groups[pairs[i].second]);
                            } else
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
                    Tspace &spc;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
                    double g2g(const Tgroup &g1, const Tgroup &g2, const std::vector<int> &index=std::vector<int>(), const std::vector<int> &jndex=std::vector<int>()) override {
#pragma GCC diagnostic pop
                        int i = &g1 - &base::spc.groups.front();
                        int j = &g2 - &base::spc.groups.front();
                        if (j<i)
                            std::swap(i,j);
                        if (base::key==Energybase::NEW) {        // if this is from the trial system,
                            double u = 0;
                            if (not base::cut(g1,g2)) {
                                for (auto &i : g1)
                                    for (auto &j : g2)
                                        u += base::i2i(i,j);
                            }
                            cache(i,j) = u;
                        }
                        return cache(i,j);                     // return (cached) value
                    }

                public:
                    NonbondedCached(const json &j, Tspace &spc) : base(j,spc), spc(spc) {
                        base::name += "EM";
                        init();
                    }

                    void init() override {
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
                    } //!< Cache pair interactions in matrix

                    double energy(Change &change) override {
                        using namespace ranges;
                        double u=0;

                        if (change) {

                            if (change.all || change.dV) {
#pragma omp parallel for reduction (+:u) schedule (dynamic) if (this->omp_enable)
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

#pragma omp parallel for reduction (+:u) schedule (dynamic) if (this->omp_enable and this->omp_g2g)
                                for (size_t i=0; i<spc.groups.size(); i++) {
                                    auto &g2 = spc.groups[i];
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
                            if (change.moved2moved)
                                for ( auto i = moved.begin(); i != moved.end(); ++i )
                                    for ( auto j=i; ++j != moved.end(); )
                                        u += g2g( base::spc.groups[*i], base::spc.groups[*j] );
                            // moved<->static
                            if (this->omp_enable and this->omp_g2g) {
                                std::vector<std::pair<int,int>> pairs( size(moved) * size(fixed) );
                                size_t cnt=0;
                                for (auto i : moved)
                                    for( auto j : fixed)
                                        pairs[cnt++] = {i,j};
#pragma omp parallel for reduction (+:u) schedule (dynamic) if (this->omp_enable and this->omp_g2g)
                                for (size_t i=0; i<pairs.size(); i++)
                                    u += g2g(spc.groups[pairs[i].first], spc.groups[pairs[i].second]);
                            } else
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
                                for (int i=0; i<d.index; i++)
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
                    bool overwrite_penalty=true; // overwrites the input penalty function
                    bool nodrift;       // avoid energy drift when upgrading penalty function
                    bool quiet;         // hold mund
                    size_t dim=0;       // number of reaction coordinate
                    size_t cnt=0;       // number of calls to `sync()`
                    size_t nupdate;     // update frequency [steps]
                    size_t samplings;
                    size_t nconv=0;
                    double udelta=0;    // total energy change of updating penalty function
                    double scale;       // scaling factor for f0
                    double f0;          // penalty increment
                    std::string file, hisfile;
                    std::vector<Tcoord> rcvec; // vector of reaction coordinate functions (length = 1 or 2)
                    std::vector<double> coord; // latest reaction coordinate (length = 1 or 2)

                    Table<int> histo;      // sampling along reaction coordinates
                    Table<double> penalty; // penalty function

                public:
                    Penalty(const json &j, Tspace &spc) : spc(spc) {
                        using namespace ReactionCoordinate;
                        name = "penalty";
                        overwrite_penalty = j.value("overwrite", true);
                        f0 = j.at("f0").get<double>();
                        scale = j.at("scale").get<double>();
                        quiet = j.value("quiet", true);
                        nupdate = j.at("update").get<double>();
                        samplings = j.value("samplings", 1);
                        nodrift = j.value("nodrift", true);
                        file = j.at("file").get<std::string>();
                        hisfile = j.value("histogram", "penalty-histogram.dat");
                        std::vector<double> binwidth, min, max;

                        if (scale<0 or scale>1)
                            throw std::runtime_error("`scale` must be in the interval [0:1]");

                        for (auto &i : j.at("coords"))
                            if (i.is_object())
                                if (i.size()==1) {
                                    std::shared_ptr<ReactionCoordinate::ReactionCoordinateBase> rc=nullptr;
                                    for (auto it=i.begin(); it!=i.end(); ++it) {

                                        try {
                                            if (it.key()=="atom")
                                                rc = std::make_shared<AtomProperty>(it.value(), spc);
                                            else if (it.key()=="molecule")
                                                rc = std::make_shared<MoleculeProperty>(it.value(), spc);
                                            else if (it.key()=="system")
                                                rc = std::make_shared<SystemProperty>(it.value(), spc);
                                            else if (it.key()=="cmcm")
                                                rc = std::make_shared<MassCenterSeparation>(it.value(), spc);
                                            if (rc!=nullptr) {
                                                if (rc->min>=rc->max || rc->binwidth<=0)
                                                    throw std::runtime_error("min<max and binwidth>0 required");
                                                rcvec.push_back(rc);
                                                binwidth.push_back( rc->binwidth );
                                                min.push_back( rc->min );
                                                max.push_back( rc->max );
                                            } else
                                                throw std::runtime_error("unknown coordinate type");

                                        } catch (std::exception &e) {
                                            throw std::runtime_error("error for reaction coordinate '"
                                                    + it.key() + "': " + e.what() + usageTip["coords=["+it.key()+"]"]  );
                                        }
                                    }
                                }
                        dim = binwidth.size();
                        if (dim<1 or dim>2)
                            throw std::runtime_error("exactly one or two coordinates required");

                        coord.resize(rcvec.size(), 0);
                        histo.reInitializer(binwidth, min, max);
                        penalty.reInitializer(binwidth, min, max);

                        std::ifstream f(MPI::prefix+file);
                        if (f) {
                            cout << "Loading penalty function '" << MPI::prefix+file << "'" << endl;
                            std::string hash;
                            f >> hash >> f0 >> samplings >> nconv;
                            for (int row=0; row<penalty.rows(); row++)
                                for (int col=0; col<penalty.cols(); col++)
                                    if (not f.eof())
                                        f >> penalty(row,col);
                                    else
                                        throw std::runtime_error("penalty file dimension mismatch");
                        }
                    }

                    virtual ~Penalty() {
                        if (overwrite_penalty) {
                            std::ofstream f(MPI::prefix + file);
                            if (f) {
                                f.precision(16);
                                f << "# " << f0 << " " << samplings << " " << nconv << "\n"
                                    << penalty.array() - penalty.minCoeff() << "\n";
                                f.close();
                            }
                        }

                        std::ofstream f2(MPI::prefix + hisfile);
                        if (f2)
                            f2 << histo << "\n";
                        // add function to save to numpy-friendly file...
                    }

                    void to_json(json &j) const override {
                        j["file"] = file;
                        j["scale"] = scale;
                        j["update"] = nupdate;
                        j["nodrift"] = nodrift;
                        j["histogram"] = hisfile;
                        j["f0_final"] = f0;
                        j["overwrite"] = overwrite_penalty;
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
                        if (change) {
                            for (size_t i=0; i<rcvec.size(); i++) {
                                coord.at(i) = rcvec[i]->operator()();
                                if ( not rcvec[i]->inRange(coord[i]) )
                                    return pc::infty;
                            }
                            penalty.to_index(coord);
                            u = penalty[coord];
                        }
                        // reaching here, `coord` always reflects
                        // the current reaction coordinate
                        return (nodrift) ? u - udelta : u;
                    }

                    /*
                     * @todo: If this is called before `energy()`, the coord
                     * is never calculated and causes undefined behavior
                     */
                    virtual void update(const std::vector<double> &c) {
                        if (++cnt % nupdate == 0 and f0>0) {
                            bool b = histo.minCoeff() >= (int)samplings;
                            if (b) {
                                double min = penalty.minCoeff(); // define minimun penalty energy
                                penalty = penalty.array() - min; // ...to zero
                                if (not quiet)
                                    cout << "Barriers/kT: penalty = " << penalty.maxCoeff()
                                        << " histogram = " << std::log(double(histo.maxCoeff())/histo.minCoeff()) << endl;
                                f0 = f0 * scale; // reduce penalty energy
                                samplings = std::ceil( samplings / scale );
                                histo.setZero();
                                udelta += -min;
                            }
                        }
                        coord = c;
                        histo[coord]++;
                        penalty[coord] += f0;
                        udelta += f0;
                    }

                    void sync(Energybase *basePtr, Change&) override {
                        // this function is called when a move is accepted
                        // or rejected, as well as when initializing the system
                        auto other = dynamic_cast<decltype(this)>(basePtr);
                        assert(other);
                        update(other->coord);
                        other->update(other->coord); // this is to keep cnt and samplings in sync

                        // some assertions...
                        assert( samplings == other->samplings );
                        assert( coord == other->coord );
                        assert( cnt == other->cnt );
                        assert( udelta == other->udelta) ;
                    } // @todo: this doubles the MPI communication
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
                using Base::quiet;

                Eigen::VectorXi weights; // array w. mininum histogram counts
                Eigen::VectorXd buffer; // receive buffer for penalty functions

                PenaltyMPI(const json &j, Tspace &spc) : Base(j,spc) {
                    weights.resize( MPI::mpi.nproc() );
                    buffer.resize( penalty.size()*MPI::mpi.nproc() ); // recieve buffer for penalty func
                }

                void update(const std::vector<double> &c) override {
                    using namespace Faunus::MPI;
                    double uold = penalty[c];
                    if (++cnt % this->nupdate == 0 and f0>0) {

                        int min = histo.minCoeff(); // if min>0 --> all RC's visited
                        MPI_Barrier(mpi.comm); // wait for all walkers to reach here
                        MPI_Allgather(&min, 1, MPI_INT, weights.data(), 1, MPI_INT, mpi.comm);

                        // if at least one walker has sampled full RC space at least `samplings` times
                        if ( weights.maxCoeff() > samplings ) { // change to minCoeff()?
                            // master collects penalty from all slaves
                            MPI_Gather(penalty.data(), penalty.size(), MPI_DOUBLE,
                                    buffer.data(), penalty.size(), MPI_DOUBLE, 0, mpi.comm);

                            // master performs the average
                            if (mpi.isMaster()) {
                                penalty.setZero();
                                for (int i=0; i<mpi.nproc(); i++)
                                    penalty += Eigen::Map<Eigen::MatrixXd>(
                                            buffer.data()+i*penalty.size(), penalty.rows(), penalty.cols() );
                                penalty = ( penalty.array() - penalty.minCoeff() ) / double(mpi.nproc());
                            }

                            // master sends the averaged penalty function to all slaves
                            MPI_Bcast(penalty.data(), penalty.size(), MPI_DOUBLE, 0, mpi.comm);
                            nconv += 1;

                            // at this point, *all* penalty functions shall be identical

                            // save penalty function to disk
                            if (mpi.isMaster()) {
                                std::ofstream f(file + ".walkersync" + std::to_string(nconv));
                                if (f) {
                                    f.precision(16);
                                    f << "# " << f0 << " " << samplings << " " << nconv << "\n"
                                        << penalty.array() << endl;
                                }
                            }

                            // save histogram to disk
                            std::ofstream f(MPI::prefix + hisfile + ".walkersync" + std::to_string(nconv));
                            if (f) {
                                f << histo << endl;
                                f.close();
                            }

                            // print information to console
                            if (min>0 and not quiet) {
                                cout << "Barriers/kT: penalty = " << penalty.maxCoeff()
                                    << " histogram = " << std::log(double(histo.maxCoeff())/histo.minCoeff()) << endl;
                            }

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

#ifdef ENABLE_POWERSASA
        /*
         * @todo:
         * - can only a subset of sasa be calculated? Note that it's the
         *   `update_coord()` function that takes up most time.
         * - delegate to GPU? In the PowerSasa paper this is mentioned
         */
        template<class Tspace>
            class SASAEnergy : public Energybase {
                public:
                    std::vector<float> sasa, radii;
                private:
                    typedef typename Tspace::Tparticle Tparticle;
                    typedef typename Tspace::Tpvec Tpvec;
                    Tspace& spc;
                    double probe; // sasa probe radius (angstrom)
                    double conc=0;// co-solute concentration (mol/l)
                    Average<double> avgArea; // average surface area
                    std::shared_ptr<POWERSASA::PowerSasa<float,Point>> ps=nullptr;

                    void updateSASA(const Tpvec &p) {
                        assert(ps != nullptr);
                        radii.resize(p.size());
                        std::transform(p.begin(), p.end(), radii.begin(),
                                [this](auto &a){ return atoms[a.id].sigma*0.5 + this->probe;});

                        ps->update_coords(spc.positions(), radii); // slowest step!

                        for (size_t i=0; i<p.size(); i++) {
                            auto &a = atoms[p[i].id];
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

                    /*
                     * @note
                     * This is not enough as the PowerSasa object contains data
                     * that also need syncing. It works due to the `update` (expensive!)
                     * call in `energy`.
                     */
                    void sync(Energybase *basePtr, Change &c) override {
                        auto other = dynamic_cast<decltype(this)>(basePtr);
                        if (other) {
                            if (c.all || c.dV) {
                                radii = other->radii;
                                sasa = other->sasa;
                            } else {
                                for (auto &d : c.groups) {
                                    int offset = std::distance(spc.p.begin(), spc.groups.at(d.index).begin());
                                    for (int j : d.atoms) {
                                        int i = j + offset;
                                        radii[i] = other->radii[i];
                                        sasa[i] = other->sasa[i];
                                    }
                                }
                            }
                        }
                    }

                public:
                    SASAEnergy(const json &j, Tspace &spc) : spc(spc) {
                        name = "sasa";
                        cite = "doi:10.1002/jcc.21844";
                        probe = j.value("radius", 1.4) * 1.0_angstrom;
                        conc = j.at("molarity").get<double>() * 1.0_molar;
                        init();
                    }

                    void init() override {
                        radii.resize( spc.p.size() );
                        std::transform( spc.p.begin(), spc.p.end(), radii.begin(),
                                [this](auto &a){ return atoms[a.id].sigma*0.5 + this->probe;} );

                        if (ps==nullptr)
                            ps = std::make_shared<POWERSASA::PowerSasa<float,Point>>(spc.positions(),radii);
                        updateSASA(spc.p);
                    }

                    double energy(Change&) override {
                        double u=0, A=0;
                        /*
                         * ideally we want to call `update` only if `key==NEW` but
                         * syncronising the PowerSasa object is difficult since it's
                         * non-copyable.
                         */
                        updateSASA(spc.p); // ideally we want 
                        for (size_t i=0; i<spc.p.size(); ++i) {
                            auto &a = atoms[ spc.p[i].id ];
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
                Example2D(const json&, Tspace &spc): i(spc.p.at(0).pos) { name = "Example2D"; }
            double energy(Change &change) override;
        };

        template<typename Tspace>
            class Hamiltonian : public Energybase, public BasePointerVector<Energybase> {
                protected:
                    double maxenergy=pc::infty; //!< Maximum allowed energy change
                    typedef typename Tspace::Tparticle Tparticle;
                    void to_json(json &j) const override {
                        for (auto i : this->vec)
                            j.push_back(*i);
                    }

                    void addEwald(const json &j, Tspace &spc) {
                        if (j.count("coulomb")==1)
                            if (j["coulomb"].count("type")==1)
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
                        typedef CombinedPairPotential<Coulomb,HardSphere<Tparticle>> PrimitiveModel;

                        Energybase::name="hamiltonian";

                        // add container overlap energy for non-cuboidal geometries
                        if (spc.geo.type not_eq Geometry::Chameleon::CUBOID)
                            push_back<Energy::ContainerOverlap<Tspace>>(spc);

                        for (auto &m : j.at("energy")) {// loop over move list
                            size_t oldsize = vec.size();
                            for (auto it=m.begin(); it!=m.end(); ++it) {
                                try {
                                    if (it.key()=="nonbonded_coulomblj")
                                        push_back<Energy::Nonbonded<Tspace,CoulombLJ>>(it.value(), spc);

                                    if (it.key()=="nonbonded_coulomblj_EM")
                                        push_back<Energy::NonbondedCached<Tspace,CoulombLJ>>(it.value(), spc);

                                    if (it.key()=="nonbonded")
                                        push_back<Energy::Nonbonded<Tspace,TabulatedPotential<typename Tspace::Tparticle>>>(it.value(), spc);
                                    
                                    if (it.key()=="nonbonded_exact")
                                        push_back<Energy::Nonbonded<Tspace,FunctorPotential<typename Tspace::Tparticle>>>(it.value(), spc);

                                    if (it.key()=="nonbonded_cached")
                                        push_back<Energy::NonbondedCached<Tspace,TabulatedPotential<typename Tspace::Tparticle>>>(it.value(), spc);

                                    if (it.key()=="nonbonded_coulombwca")
                                        push_back<Energy::Nonbonded<Tspace,CoulombWCA>>(it.value(), spc);

                                    if (it.key()=="nonbonded_pm" or it.key()=="nonbonded_coulombhs")
                                        push_back<Energy::Nonbonded<Tspace,PrimitiveModel>>(it.value(), spc);

                                    if (it.key()=="nonbonded_pmwca")
                                        push_back<Energy::Nonbonded<Tspace,PrimitiveModelWCA>>(it.value(), spc);

                                    if (it.key()=="bonded")
                                        push_back<Energy::Bonded<Tspace>>(it.value(), spc);

                                    if (it.key()=="customexternal")
                                        push_back<Energy::CustomExternal<Tspace>>(it.value(), spc);

                                    if (it.key()=="akesson")
                                        push_back<Energy::ExternalAkesson<Tspace>>(it.value(), spc);

                                    if (it.key()=="confine")
                                        push_back<Energy::Confine<Tspace>>(it.value(), spc);

                                    if (it.key()=="constrain")
                                        push_back<Energy::Constrain>(it.value(), spc);

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

                                    if (it.key()=="maxenergy") {
                                        maxenergy = it.value().get<double>();
                                        continue;
                                    }

                                    if (vec.size()==oldsize)
                                        throw std::runtime_error("unknown term");

                                } catch (std::exception &e) {
                                    throw std::runtime_error("Error adding energy '" + it.key() + "': " + e.what() + usageTip[it.key()]);
                                }
                            }
                        }
                    }

                    double energy(Change &change) override {
                        double du=0;
                        for (auto i : this->vec) {
                            i->key=key;
                            i->timer.start();
                            du += i->energy(change);
                            i->timer.stop();
                            if (du>=maxenergy)
                                break; // stop summing energies
                        }
                        return du;
                    } //!< Energy due to changes

                    void init() override {
                        for (auto i : this->vec)
                            i->init();
                    }

                    void sync(Energybase* basePtr, Change &change) override {
                        auto other = dynamic_cast<decltype(this)>(basePtr);
                        if (other)
                            if (other->size()==size()) {
                                for (size_t i=0; i<size(); i++)
                                    this->vec[i]->sync( other->vec[i].get(), change );
                                return;
                            }
                        throw std::runtime_error("hamiltonian mismatch");
                    }

            }; //!< Aggregates and sum energy terms

    }//namespace
}//namespace
