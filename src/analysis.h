#pragma once
#include "space.h"

namespace Faunus {

    namespace Analysis {

        class AnalysisBase {
            private:
                inline virtual void _to_json(json &j) const {};
                inline virtual void _from_json(const json &j) {};
                virtual void _sample()=0;
                int stepcnt=0;
                TimeRelativeOfTotal<std::chrono::microseconds> timer;

            protected:
                int steps=0; //!< Sample interval (do not modify)
                int cnt=0;   //!< number of samples

            public:
                std::string name; //!< descriptive name
                std::string cite; //!< reference, url, doi etc. describing the analysis

                inline void to_json(json &j) const {
                    assert( !name.empty() );
                    auto& _j = j[name];
                    _to_json(_j);
                    _j["nstep"] = steps;
                    _j["samples"] = cnt;
                    if (!cite.empty())
                        _j["citations"] = cite;
                    _j["relative time"] = timer.result();
                } //!< JSON report w. statistics, output etc.

                inline void from_json(const json &j) {
                    auto &_j = j.at(name);
                    steps = _j.value("nstep", 0);
                    _from_json(_j);
                } //!< configure from json object

                inline void sample()
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

        };

        inline void to_json(json &j, const AnalysisBase &b) { b.to_json(j); }
        inline void from_json(const json &j, AnalysisBase &b) { b.from_json(j); }

        class SystemEnergy : public AnalysisBase {
            private:
                std::string file;
                std::ofstream f;
                std::function<double()> energyFunc;
                inline void _sample() override { f << cnt*steps << " " << energyFunc() << "\n"; }
                inline void _from_json(const json &j) override {
                    file = j.at("file");
                    if (f)
                        f.close();
                    f.open(file);
                    if (!f)
                        throw std::runtime_error(name + ": cannot open output file " + file);
                }

            public:
                template<class Tenergy>
                    SystemEnergy( const json &j, Tenergy &pot ) {
                        name = "systemenergy";
                        from_json(j);
                        energyFunc = [&pot]() {
                            Change c;
                            c.all=true;
                            return pot.energy(c);
                        };
                    }
        }; //!< Save system energy to disk. Keywords: `nstep`, `file`.

        /** @brief Example analysis */
        template<class T, class Enable = void>
            struct _analyse {
                void sample(T &p) {
                    std::cout << "not a dipole!" << std::endl;
                } //!< Sample
            }; // primary template

        /** @brief Example analysis */
        template<class T>
            struct _analyse<T, typename std::enable_if<std::is_base_of<Dipole, T>::value>::type> {
                void sample(T &p) {
                    std::cout << "dipole!" << std::endl;
                } //!< Sample
            }; // specialized template

    }//namespace

}//namespace
