#pragma once
#include "core.h"
#include "units.h"

namespace Faunus {
    /**
     * @brief General properties for atoms
     */
    struct AtomData {
        public:
            int _id=-1;
            std::string name;     //!< Name
            double eps=0;         //!< LJ epsilon [kJ/mol] (pair potentials should convert to kT)
            double activity=0;    //!< Chemical activity [mol/l]
            double alphax=0;      //!< Excess polarisability (unit-less)
            double charge=0;      //!< Particle charge [e]
            double dp=0;          //!< Translational displacement parameter [angstrom]
            double dprot=0;       //!< Rotational displacement parameter [degrees]
            Point  mu={1,0,0};    //!< Dipole moment unit vector
            double mulen=0;       //!< Dipole moment scalar [eÃ]
            double mw=1;          //!< Weight
            Point  scdir={1,0,0}; //!< Sphero-cylinder direction
            double sclen=0;       //!< Sphere-cylinder length [angstrom]
            double sigma=0;       //!< Diameter for e.g Lennard-Jones etc. [angstrom]
            double hdr=0;         //!< Hydrodynamic radius [angstrom]
            double tension=0;     //!< Surface tension [kT/Å^2]
            double tfe=0;         //!< Transfer free energy [J/mol/angstrom^2/M]
            double squarewell_depth;     //!< Depth of square-well potential [kJ/mol] (pair potentials should convert to kT)
            double hz_eps;               //!< Strength of Hertz potential [kJ/mol] (pair potentials should convert to kT)
            double squarewell_threshold; //!< Threshold for square-well potential [angstrom]
            int& id(); //!< Type id
            const int& id() const; //!< Type id
            bool hydrophobic=false;  //!< Is the particle hydrophobic?
            bool implicit=false; //!< Is the particle implicit (e.g. proton)?
    };

    void to_json(json& j, const AtomData &a);
    void from_json(const json& j, AtomData& a);

    /**
     * @brief Construct vector of atoms from json array
     *
     * Items are added to existing items while if an item
     * already exists, it will be overwritten.
     */
    void from_json(const json& j, std::vector<AtomData> &v);

    extern std::vector<AtomData> atoms; //!< Global instance of atom list

    template<class Trange>
        auto findName(Trange &rng, const std::string &name) {
            return std::find_if( rng.begin(), rng.end(), [&name](auto &i){ return i.name==name; });
        } //!< Returns iterator to first element with member `name` matching input

    template<class Trange>
        std::vector<int> names2ids(Trange &rng, const std::vector<std::string> &names) {
            std::vector<int> index;
            index.reserve(names.size());
            for (auto &n : names) {
                // wildcard selecting all id's
                if (n=="*") {
                    index.resize( rng.size() );
                    std::iota(index.begin(), index.end(), 0);
                    return index;
                }
                auto it = findName(rng, n);
                if (it!=rng.end())
                    index.push_back(it->id());
                else
                    throw std::runtime_error("name '" + n + "' not found");
            }
            return index;
        } //!< Convert vector of names into vector of id's from Trange (exception if not found)

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] AtomData") {
        using doctest::Approx;

        json j = R"({ "atomlist" : [
             { "A": { "r":1.1, "pactivity":2 } },
             { "B": { "activity":0.2, "eps":0.05, "dp":9.8, "dprot":3.14, "mw":1.1, "tfe":0.98, "tension":0.023 } }
             ]})"_json;

        pc::temperature = 298.15;

        atoms = j["atomlist"].get<decltype(atoms)>();
        auto &v = atoms; // alias to global atom list

        CHECK(v.size()==2);
        CHECK(v.front().id()==0);
        CHECK(v.front().name=="A"); // alphabetic order in std::map
        CHECK(v.front().sigma==Approx(2*1.1e-10_m));
        CHECK(v.front().activity==Approx(0.01_molar));
        CHECK(v.back().tfe==Approx(0.98_kJmol/(1.0_angstrom*1.0_angstrom*1.0_molar)));

        AtomData a = json(v.back()); // AtomData -> JSON -> AtomData

        CHECK(a.name=="B");
        CHECK(a.id()==1);
        CHECK(a.activity==Approx(0.2_molar));
        CHECK(a.eps==Approx(0.05_kJmol));
        CHECK(a.dp==Approx(9.8));
        CHECK(a.dprot==Approx(3.14));
        CHECK(a.mw==Approx(1.1));
        CHECK(a.tfe==Approx(0.98_kJmol/1.0_angstrom/1.0_angstrom/1.0_molar));
        CHECK(a.tension==Approx(0.023_kJmol/1.0_angstrom/1.0_angstrom));

        auto it = findName(v, "B");
        CHECK( it->id() == 1 );
        it = findName(v, "unknown atom");
        CHECK( it==v.end() );
    }
#endif

} // end of Faunus namespace
