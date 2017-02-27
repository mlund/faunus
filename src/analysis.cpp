#include <faunus/analysis.h>
#include <faunus/species.h>
#include <faunus/energy.h>
#include <faunus/space.h>
#include <faunus/inputfile.h>
#include <faunus/geometry.h>
#include <faunus/textio.h>

namespace Faunus
{

  namespace Analysis
  {
    AnalysisBase::AnalysisBase() : w(30), cnt(0)
    {
        stepcnt = 0;
    }

    AnalysisBase::AnalysisBase( Tmjson &j, string name ) : w(30), cnt(0), name(name)
    {
        if (!j.is_object())
            std::runtime_error("Analysis JSON entry must be of type object");
        steps = j.value("nstep", 0);
        stepcnt = 0;
    }

    AnalysisBase::~AnalysisBase() {}

    string AnalysisBase::_info() { return string(); }

    Tmjson AnalysisBase::_json() { return Tmjson(); }

    void AnalysisBase::_sample()
    {
        assert(!"We should never reach here -- implement _sample() function");
        /* make pure virtual! */
    }

    void AnalysisBase::sample()
    {
        stepcnt++;
        if ( stepcnt == steps )
        {
            cnt++;
            stepcnt = 0;
            timer.start();
            _sample();
            timer.stop();
        }
    }

    void AnalysisBase::_test( UnitTest &t ) {}

    void AnalysisBase::test( UnitTest &t ) { _test(t); }

    string AnalysisBase::info()
    {
        using namespace textio;
        std::ostringstream o;
        if ( cnt > 0 )
        {
            o << header("Analysis: " + name);
            if ( !cite.empty())
                o << pad(SUB, w, "Reference:") << cite << "\n";

            o << pad(SUB, w, "Sample interval") << steps << "\n";

            if ( cnt > 0 )
            {
                o << pad(SUB, w, "Number of sample events") << cnt << "\n";
                double time = timer.result();
                if ( time > 1e-3 )
                    o << pad(SUB, w, "Relative time") << time << "\n";
            }
            o << _info();
        }
        return o.str();
    }

    Tmjson AnalysisBase::json()
    {
        Tmjson j;
        if ( !name.empty())
            if ( cnt > 0 )
            {
                j[name] = {
                    {"nstep", steps},
                    {"samples", cnt},
                    {"relative time", timer.result()}
                };
                if ( !cite.empty())
                {
                    j[name]["citation"] = cite;
                }
                j = merge(j, _json());
            }
        return j;
    }

    void BilayerStructure::_test( UnitTest &t )
    {
        t("bilayer_order", S.avg());
        t("bilayer_area", A.avg());
    }

    Tmjson MeanForce::_json()
    {
        if ( mf1.cnt>0 && mf2.cnt>0 )
            return {
                { name,
                    {
                        { "groups", {g1, g2} },
                        { "meanforce", { mf1.avg(), mf2.avg() } },
                        { "forceunit", "kT/angstrom" }
                    }
                }
            };
        return Tmjson();
    }

    void MeanForce::_sample() {
        func();
        cout << mf1.cnt << " " << mf2.cnt << "\n";
    }

    void SystemEnergy::_sample() {
        f << energy() << "\n"; 
    }

    void PairFunctionBase::_sample()
    {
        for (auto &d : datavec)
            update(d);
    }

    Tmjson PairFunctionBase::_json()
    {
        Tmjson j;
        auto &_j = j[name];
        for (auto &d : datavec)
            _j[ d.name1+"-"+d.name2 ] = {
                { "dr", d.dr },
                { "file", d.file },
                { "dim", d.dim },
		{ "Rhyper", d.Rhypersphere }
            };
        return j;
    }

    PairFunctionBase::PairFunctionBase( Tmjson j, string name ) : AnalysisBase(j, name) {
        try {
            for (auto &i : j.at("pairs"))
                if (i.is_object())
                {
                    data d;
                    d.file = i.at("file");
                    d.name1 = i.at("name1");
                    d.name2 = i.at("name2");
                    d.dim = i.value("dim", 3);
                    d.dr = i.value("dr", 0.1);
                    d.hist.setResolution(d.dr);
		    d.Rhypersphere = i.value("Rhyper", -1.0);
                    datavec.push_back( d );
                }
        }
        catch(std::exception& e) {
            throw std::runtime_error(name + ": " + e.what());
        }

        if (datavec.empty())
            std::cerr << name + ": no sample sets defined for analysis\n";
    }

    PairFunctionBase::~PairFunctionBase()
    {
        for (auto &d : datavec)
            d.hist.save( d.file );
    }

    void CombinedAnalysis::sample()
    {
        cnt++;
        for ( auto i : v )
            i->sample();
    }

    string CombinedAnalysis::info()
    {
        std::ostringstream o;
        for ( auto i : v )
            o << i->info();
        return o.str();
    }

    string CombinedAnalysis::_info() { return string(); }

    void CombinedAnalysis::_sample() {}

    void CombinedAnalysis::test( UnitTest &test )
    {
        for ( auto i : v )
            i->test(test);
    }

    Tmjson CombinedAnalysis::json()
    {
        Tmjson js;
        for ( auto i : v )
            js = merge(js, i->json());
        return js;
    }

    CombinedAnalysis::~CombinedAnalysis()
    {
        if (cnt>0) {
            std::ofstream f(jsonfile);
            if ( f )
                f << std::setw(4) << json() << std::endl;
        }
    }

  }//Analysis namespace
}//Faunus namespace
