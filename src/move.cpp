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
            cnt++;
            change.clear();
            _move(change);
            if (change.empty())
                timer.stop();
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

        double Movebase::bias(Change &c, double uold, double unew) {
            return 0; // du
        }

        void Movebase::_accept(Change &) {}

        void Movebase::_reject(Change &) {}

        void from_json(const json &j, Movebase &m) {
            m.from_json( j );
        }

        void to_json(json &j, const Movebase &m) {
            assert( !m.name.empty() );
            m.to_json(j[m.name]);
        }
    }
}
