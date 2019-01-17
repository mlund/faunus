#include "atomdata.h"

namespace Faunus {

    int &AtomData::id() { return _id; }

    const int &AtomData::id() const { return _id; }

    void to_json(json &j, const AtomData &a) {
        auto& _j = j[a.name];
        _j = {
            {"activity", a.activity / 1.0_molar},
            {"alphax", a.alphax}, {"q", a.charge},
            {"dp", a.dp / 1.0_angstrom}, {"dprot", a.dprot / 1.0_rad},
            {"eps", a.eps / 1.0_kJmol}, {"mw", a.mw},
            {"sigma", a.sigma / 1.0_angstrom},
            {"tension", a.tension * 1.0_angstrom*1.0_angstrom / 1.0_kJmol},
            {"tfe", a.tfe * 1.0_angstrom*1.0_angstrom * 1.0_molar / 1.0_kJmol},
            {"mu", a.mu}, {"mulen",a.mulen},
            {"scdir", a.scdir}, {"sclen", a.sclen},
            {"id", a.id()}
        };
        if (a.hydrophobic)
            _j["hydrophobic"] = a.hydrophobic;
    }

    void from_json(const json &j, AtomData &a) {
        if (j.is_object()==false || j.size()!=1)
            throw std::runtime_error("Invalid JSON data for AtomData");
        for (auto it : j.items()) {
            a.name = it.key();
            xjson val = it.value();
            a.activity = val.value("activity", a.activity) * 1.0_molar;
            a.alphax   = val.value("alphax", a.alphax);
            a.charge   = val.value("q", a.charge);
            a.dp       = val.value("dp", a.dp) * 1.0_angstrom;
            a.dprot    = val.value("dprot", a.dprot) * 1.0_rad;
            a.eps      = val.value("eps", a.eps) * 1.0_kJmol;
            a.id()     = val.value("id", a.id());
            a.mu       = val.value("mu", a.mu);
            a.mulen    = val.value("mulen", a.mulen);
            a.scdir    = val.value("scdir", a.scdir);
            a.sclen    = val.value("sclen", a.sclen);
            a.mw       = val.value("mw", a.mw);
            a.sigma    = val.value("sigma", 0.0) * 1.0_angstrom;
            if (fabs(a.sigma)<1e-20)
                a.sigma = 2.0*val.value("r", 0.0) * 1.0_angstrom;
            a.tension  = val.value("tension", a.tension) * 1.0_kJmol / (1.0_angstrom*1.0_angstrom);
            a.tfe      = val.value("tfe", a.tfe) * 1.0_kJmol / (1.0_angstrom*1.0_angstrom*1.0_molar);
            a.hydrophobic = val.value("hydrophobic", false);
            if (not val.empty()) // throw exception of unused/unknown keys are passed
                throw std::runtime_error("unused key(s) for atom '"s + a.name + "':\n"
                        + val.dump() + usageTip["atomlist"]);
        }
    }

    void from_json(const json &j, std::vector<AtomData> &v) {
        auto it = j.find("atomlist");
        json _j = ( it==j.end() ) ? j : *it;
        v.reserve( v.size() + _j.size() );

        for ( auto &i : _j ) {
            if ( i.is_string() ) // treat ax external file to load
                from_json( openjson(i.get<std::string>()), v );
            else if ( i.is_object() ) {
                AtomData a = i;
                auto it = findName( v, a.name );
                if ( it==v.end() )
                    v.push_back( a ); // add new atom
                else
                    *it = a;
            }
        }
        for (size_t i=0; i<v.size(); i++)
            v[i].id() = i; // id must match position in vector
    }

    std::vector<AtomData> atoms;

} // end of Faunus
