#pragma once
#include "space.h"

namespace Faunus {

    namespace Analysis {

        class Analysisbase {
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
                    auto &_j = j[name];
                    if (cnt>0) {
                        _j["relative time"] = timer.result();
                        _j["nstep"] = steps;
                        _j["samples"] = cnt;
                    }
                    if (!cite.empty())
                        _j["citation"] = cite;
                    _to_json(_j);
                } //!< JSON report w. statistics, output etc.

                inline void from_json(const json &j) {
                    steps = j.value("nstep", 0);
                    _from_json(j);
                } //!< configure from json object

                inline virtual void sample()
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

        inline void to_json(json &j, const Analysisbase &base) {
            base.to_json( j );
        }

        class SystemEnergy : public Analysisbase {
            private:
                std::string file;
                std::ofstream f;
                std::function<double()> energyFunc;
                Average<double> uavg; //!< mean energy
                double uinit;

                inline void _sample() override {
                    double u = energyFunc();
                    uavg+=u;
                    f << cnt*steps << " " << u << "\n";
                }

                inline void _to_json(json &j) const override {
                    j["file"] = file;
                    j["init"] = uinit;
                    j["final"] = energyFunc();
                    if (cnt>0)
                        j["mean"] = uavg.avg();
                }

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
                            Change change;
                            change.all = true;
                            return pot.energy(change);
                        };
                        uinit = energyFunc(); // initial energy
                    }
        }; //!< Save system energy to disk. Keywords: `nstep`, `file`.

        class SaveState : public Analysisbase {
            private:
                std::function<void(std::string)> writeFunc = nullptr;
                std::string file;
                inline void _to_json(json &j) const override {
                    j["file"] = file;
                }
            public:
                template<class Tspace>
                    SaveState(const json &j, Tspace &spc) {
                        using std::ref;
                        using std::placeholders::_1;

                        name = "savestate";
                        from_json(j);
                        steps = j.value("nstep", -1);
                        file = j.at("file");
                        std::string suffix = file.substr(file.find_last_of(".") + 1);
                        if ( suffix == "aam" )
                            writeFunc = std::bind(
                                    []( std::string file, Tspace &s ) { FormatAAM::save(file, s.p); },
                                    _1, std::ref(spc));
                        if ( suffix == "pqr" )
                            writeFunc = std::bind(
                                    []( std::string file, Tspace &s ) { FormatPQR::save(file, s.p, s.geo.getLength()); },
                                    _1, std::ref(spc));
                        if ( suffix == "state" )
                            writeFunc = [&spc](const std::string &file) {
                                std::ofstream f(file);
                                if (f)
                                    f << std::setw(4) << json(spc);
                            };
                    }

                ~SaveState() {
                    if (steps==-1)
                        _sample();
                }

                void _sample() override {
                    writeFunc(file);
                }

        };

        class CombinedAnalysis : public BasePointerVector<Analysisbase> {
            public:
                template<class Tspace, class Tenergy>
                    CombinedAnalysis(const json &j, Tspace &spc, Tenergy &pot) {
                        for (auto &m : j.at("analysis")) // loop over move list
                            for (auto it=m.begin(); it!=m.end(); ++it) {
                                if (it.key()=="systemenergy")
                                    push_back<SystemEnergy>(it.value(), pot);
                                if (it.key()=="savestate")
                                    push_back<SaveState>(it.value(), spc);
                                // additional analysis go here...
                            }
                    }

                inline void sample() {
                    for (auto i : this->vec)
                        i->sample();
                }

        }; //!< Aggregates analysis

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
