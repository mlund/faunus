#pragma once
#include "core.h"

namespace Faunus {

    /**
     * @brief General properties for molecules
     */
    template<class Tpvec>
        class MoleculeData {
            private:
                int _id=-1;
            public:
                int& id() { return _id; } //!< Type id
                const int& id() const { return _id; } //!< Type id

                int Ninit = 0;             //!< Number of initial molecules
                std::string name;          //!< Molecule name
                std::string structure;     //!< Structure file (pqr|aam|xyz)
                bool atomic=false;         //!< True if atomic group (salt etc.)
                bool rotate=true;          //!< True if molecule should be rotated upon insertion
                bool keeppos=false;        //!< Keep original positions of `structure`
                double activity=0;         //!< Chemical activity (mol/l)
                Point insdir = {1,1,1};    //!< Insertion directions
                Point insoffset = {0,0,0}; //!< Insertion offset

                std::vector<int> atoms;    //!< Sequence of atoms in molecule (atom id's)
                std::vector<Tpvec> conformations;           //!< Conformations of molecule
                std::discrete_distribution<> confDist;      //!< Weight of conformations

                /**
                 * @brief Store a single conformation
                 * @param vec Vector of particles
                 * @param weight Relative weight of conformation (default: 1)
                 */
                void addConformation( const Tpvec &vec, double weight = 1 )
                {
                    if ( !conformations.empty())
                    {     // resize weights
                        auto w = confDist.probabilities();// (defaults to 1)
                        w.push_back(weight);
                        confDist = std::discrete_distribution<>(w.begin(), w.end());
                    }
                    conformations.push_back(vec);
                    assert(confDist.probabilities().size() == conformations.size());
                }

                void loadConformation(const std::string &file)
                {
                    Tpvec v;
                    std::string suffix = structure.substr(file.find_last_of(".") + 1);
                    if ( suffix == "pqr" )
                        FormatPQR::load(structure, v);
                    //if ( suffix == "xyz" )
                    //FormatXYZ::load(structure, v);
                    if ( !v.empty())
                    {
                        if ( keeppos == false )
                            Geometry::cm2origo(
                                    Geometry::Sphere(1e20), v); // move to origo
                        pushConformation(v);        // add conformation
                        for ( auto &p : v )           // add atoms to atomlist
                            atoms.push_back(p.id);
                    }
                    if ( v.empty() )
                        throw std::runtime_error("Structure " + structure + " not loaded. Filetype must be .aam/.pqr/.xyz");
                }
        }; // end of class

    template<class T>
        void to_json(json& j, const MoleculeData<T> &a) {
            auto& _j = j[a.name];
            _j["activity"] = a.activity / 1.0_molar;
            _j["atomic"] = a.atomic;
            _j["id"] = a.id();
            _j["insdir"] = a.insdir;
            _j["insoffset"] = a.insoffset;
            _j["keeppos"] = a.keeppos;
        }

    template<class T>
        void from_json(const json& j, MoleculeData<T>& a) {
            if (j.is_object()==false || j.size()!=1)
                throw std::runtime_error("Invalid JSON data for AtomData");
            for (auto it=j.begin(); it!=j.end(); ++it) {
                a.name = it.key();
                auto& val = it.value();
                a.activity = val.value("activity", a.activity) * 1.0_molar;
                a.atomic = val.value("atomic", a.atomic);
                a.id() = val.value("id", a.id());
                a.insdir = val.value("insdir", a.insdir);
                a.insoffset = val.value("insoffset", a.insoffset);
                a.keeppos = val.value("keeppos", a.keeppos);
            }
        }

    template<typename Tpvec>
        static std::vector<MoleculeData<Tpvec>> molecules = {}; //!< Global instance of molecule list

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] MoleculeData") {
        using doctest::Approx;

        json j = {
            { "moleculelist",
                {
                    { "B",
                        {
                            {"activity",0.2}, {"atomic",true},
                            {"insdir",{0.5,0,0}}, {"insoffset",{-1.1, 0.5, 10}}
                        }
                    },
                    { "A", { {"atomic",false} } }
                }
            }
        };
        typedef Particle<Radius, Charge, Dipole, Cigar> T;
        typedef std::vector<T> Tpvec;
        typedef MoleculeData<Tpvec> Tmoldata;

        molecules<Tpvec> = j["moleculelist"].get<decltype(molecules<Tpvec>)>(); // fill global instance
        auto &v = molecules<Tpvec>; // reference to global molecule vector

        //std::vector<Tmoldata> v = j["moleculelist"];

        CHECK(v.size()==2);
        CHECK(v.front().id()==0);
        CHECK(v.front().name=="A"); // alphabetic order in std::map
        CHECK(v.front().atomic==false);

        MoleculeData<T> m = json(v.back()); // moldata --> json --> moldata

        CHECK(m.name=="B");
        CHECK(m.id()==1);
        CHECK(m.activity==Approx(0.2_molar));
        CHECK(m.atomic==true);
        CHECK(m.insdir==Point(0.5,0,0));
        CHECK(m.insoffset==Point(-1.1,0.5,10));
    }
#endif

}//namespace

