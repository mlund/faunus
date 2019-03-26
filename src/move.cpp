#include "core.h"
#include "move.h"

namespace Faunus {
    namespace Move {

        Random Movebase::slump; // static instance of Random (shared for all moves)

        void Movebase::from_json(const json &j) {
            auto it = j.find("repeat");
            if (it!=j.end()) {
                if (it->is_number())
                    repeat = it->get<double>();
                else
                if (it->is_string())
                    if (it->get<std::string>()=="N")
                        repeat = -1;
            }
            _from_json(j);
            if (repeat<0)
                repeat=0;
        }

        void Movebase::to_json(json &j) const {
            _to_json(j);
            if (timer_move.result() > 0.01) // only print if more than 1% of the time
                j["relative time (without energy calc)"] = timer_move.result();
            if (timer.result() > 0.01) // only print if more than 1% of the time
                j["relative time"] = timer.result();
            j["acceptance"] = double(accepted)/cnt;
            j["repeat"] = repeat;
            j["moves"] = cnt;
            if (!cite.empty())
                j["cite"] = cite;
            _roundjson(j, 3);
        }

        void Movebase::move(Change &change) {
            timer.start();
            timer_move.start();
            cnt++;
            change.clear();
            _move(change);
            if (change.empty())
                timer.stop();
            timer_move.stop();
        }

        void Movebase::accept(Change &c) {
            accepted++;
            _accept(c);
            timer.stop();
        }

        void Movebase::reject(Change &c) {
            rejected++;
            _reject(c);
            timer.stop();
        }

        double Movebase::bias(Change&, double, double) {
            return 0; // du
        }

        void Movebase::_accept(Change&) {}

        void Movebase::_reject(Change&) {}

        void from_json(const json &j, Movebase &m) {
            m.from_json( j );
        }

        void to_json(json &j, const Movebase &m) {
            assert( !m.name.empty() );
            m.to_json(j[m.name]);
        }

        void ChainRotationMovebase::_from_json(const json &j) {
            molname = j.at("molecule");
            dprot = j.at("dprot");
            allow_small_box = j.value("skiplarge", true); // todo rename the json attribute and make false default
        }

        void ChainRotationMovebase::_to_json(json &j) const {
            using namespace u8;
            j = {
                    {"molecule", molname}, {"dprot", dprot},
                    {u8::rootof + u8::bracket("r_cm" + u8::squared), std::sqrt(msqdispl.avg())}
            };
            if(small_box_encountered > 0) {
                j["skipped"] = double(small_box_encountered) / cnt;  // todo rename the json attribute
            }
            _roundjson(j,3);
        }

        void ChainRotationMovebase::_move(Change &change) {
            permit_move = true;
            sqdispl = 0;
            if (std::fabs(dprot) > 1e-9) {
                if (select_segment() > 0) {
                    double angle = dprot * (slump() - 0.5);
                    rotate_segment(angle);
                    store_change(change);
                }
            }
        }

        void ChainRotationMovebase::_accept(Change&) { msqdispl += sqdispl; }
        void ChainRotationMovebase::_reject(Change&) { msqdispl += 0; }
        double ChainRotationMovebase::bias(Change&, double, double) { return permit_move ? 0 : pc::infty; }
    }
}

