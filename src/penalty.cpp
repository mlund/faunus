#include "penalty.h"

namespace Faunus {
    namespace ReactionCoordinate {

        void ReactionCoordinateBase::_to_json(json&) const {}

        double ReactionCoordinateBase::normalize(double) const { return 1.; }

        double ReactionCoordinateBase::operator()() {
            assert(f!=nullptr);
            return f();
        }

        bool ReactionCoordinateBase::inRange(double coord) const {
            return (coord>=min && coord<=max);
        }

        void to_json(json &j, const ReactionCoordinateBase &r) {
            j = { {"range",{r.min,r.max}}, {"resolution",r.binwidth} };
            r._to_json(j);
        }

        void from_json(const json &j, ReactionCoordinateBase &r) {
            r.binwidth = j.value("resolution", 0.5);
            auto range = j.value("range", std::vector<double>({0,0}));
            if (range.size()==2)
                if (range[0]<=range[1]) {
                    r.min = range[0];
                    r.max = range[1];
                    return;
                }
            throw std::runtime_error(r.name + ": 'range' require two numbers: [min, max>=min]");
        }

        void SystemProperty::_to_json(json &j) const {
            j["property"] = property;
        }

        void AtomProperty::_to_json(json &j) const {
            j["property"] = property;
            j["index"] = index;
            if (dir.squaredNorm()>1e-9)
                j["dir"] = dir;
        }

        void MoleculeProperty::_to_json(json &j) const {
            j["property"] = property;
            j["index"] = index;
            if (dir.squaredNorm()>1e-9)
                j["dir"] = dir;
            if (indexes.size()>=2)
                j["indexes"] = indexes;
        }

        double MassCenterSeparation::normalize(double coord) const {
            int dim=dir.sum();
            if (dim==2) return 1/(2*pc::pi*coord);
            if (dim==3) return 1/(4*pc::pi*coord*coord);
            return 1.0;
        }

        void MassCenterSeparation::_to_json(json &j) const {
            j["dir"] = dir;
            j["indexes"] = indexes;
            j["type"] = type;
        }

    }
}

