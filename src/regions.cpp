#include "regions.h"
#include "space.h"

namespace Faunus {
namespace Region {

RegionBase::RegionBase(RegionType type) : type(type){};

/**
 * The opposite of `createRegion()`
 *
 * Having this in place, std::vector of shared pointers can
 * automatically be serialized to json
 */
void to_json(json &j, const std::shared_ptr<RegionBase> &region) {
    j[region->name] = json::object();
    region->to_json(j[region->name]);
}

/**
 * Expects an object where KEY is an arbitrary, user-defined name and
 * the VALUE is another object defining rhe region type and other
 * required properties. Example:
 *
 * ~~~ yaml
 * mysubspace:
 *     type: within
 *     threshold: 7
 *     indexes: [1,2]
 *     com: true
 * ~~~
 */
std::shared_ptr<RegionBase> createRegion(const json &j, Space &spc) {
    std::shared_ptr<RegionBase> region;
    if (j.is_object() and j.size() == 1) { // object with exactly one key expected
        for (auto &it : j.items()) {       // this is just to get the key and value
            const std::string &type = it.value().at("type");
            if (type == WithinGroups::type_name)
                region = std::make_shared<WithinGroups>(it.value(), spc);

            // add `else if` statements for future regions here

            if (region) // good, region properly constructed
                region->name = it.key();
            else // not so good
                throw std::runtime_error("unknown region for "s + it.key());
        }
    }
    assert(region); // temporary; implement exceptions
    return region;
}

void WithinGroups::to_json(json &j) const {
    j = {{"type", type_name}, {"threshold", sqrt(threshold2)}, {"com", com}, {"indexes", indexes}};
}

const std::string WithinGroups::type_name = "within";

WithinGroups::WithinGroups(const json &j, Space &spc) : RegionBase(WITHIN), spc(spc) {
    threshold2 = std::pow(j.at("threshold").get<double>(), 2);
    com = j.value("com", false);

    // the following checks first for `indexes` (list of molecule index) and,
    // if not found, for `molecules` (list of molecule names). An exception
    // is thrown if:
    // - indexes are out of range
    // - if unknown molecules are given
    // - neither `indexes` nor `molecules` are found
    // - if com=true and one of the selected molecules are atomic
    auto it = j.find("indexes");
    if (it != j.end()) {
        indexes = it->get<decltype(indexes)>();
        if (not indexes.empty())
            if (*std::max_element(indexes.begin(), indexes.end()) >= spc.groups.size())
                throw std::runtime_error("indexes out of range");
    } else { // no indexes given, now check for molecules
        auto molids = Faunus::names2ids(Faunus::molecules, j.at("molecules").get<std::vector<std::string>>());
        for (size_t i = 0; i < spc.groups.size(); i++) {
            if (std::find(molids.begin(), molids.end(), spc.groups[i].id) != molids.end())
                indexes.push_back(i);
        }
    }

    // if COM requested, atomic groups are forbidden
    if (com)
        for (size_t i : indexes)
            if (spc.groups.at(i).atomic)
                throw std::runtime_error("com cannot be used with atomic groups");

    if (indexes.empty())
        std::cerr << "warning: no molecules selected for region" << endl;
}

bool WithinGroups::isInside(const Point &a) const {
    for (size_t i : indexes) {       // loop over user defined group index
        auto &g = spc.groups.at(i);  // ref. to current group
        if (com and not g.empty()) { // check only with mass-center
            assert(g.atomic == false);
            if (spc.geo.sqdist(a, g.cm) < threshold2)
                return true;
        } else {
            for (Particle &b : g) // loop over active particles in group
                if (spc.geo.sqdist(a, b.pos) < threshold2)
                    return true;
        }
    }
    return false;
}

double WithinGroups::volume() const { return (com) ? 4 * pc::pi / 3 * std::pow(threshold2, 1.5) : -1; }

} // namespace Region
} // namespace Faunus
