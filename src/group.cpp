#include "group.h"
#include <nlohmann/json.hpp>
#include <Eigen/Geometry>
#include "aux/eigensupport.h"

namespace Faunus {

template <class T> Group<T>::Group(Group<T> &o) : base(o.begin(), o.trueend()) { *this = operator=(o); }

template <class T> Group<T>::Group(const Group<T> &o) : base(o.begin(), o.trueend()) { *this = operator=(o); }

template <class T> Group<T>::Group(Group<T>::iter begin, Group<T>::iter end) : base(begin,end) {}

template <class T> Group<T> &Group<T>::operator=(const Group<T> &o) {
    if (&o == this)
        return *this;
    shallowcopy(o);
    if (o.begin()!=begin())
        std::copy(o.begin(), o.trueend(), begin()); // copy all particle data
    return *this;
}

template <class T> Group<T> &Group<T>::shallowcopy(const Group<T> &o) {
    if (&o != this) {
        if (this->capacity() != o.capacity())
            throw std::runtime_error("Group::shallowcopy: capacity mismatch");
        this->resize(o.size());
        id = o.id;
        atomic = o.atomic;
        compressible = o.compressible;
        cm = o.cm;
        confid = o.confid;
    }
    return *this;
}

template <class T> bool Group<T>::contains(const T &a, bool include_inactive) const {
    int size = (include_inactive ? this->capacity() : this->size());
    if (size > 0) {
        int d = &a - &(*(this->begin()));
        if (d >= 0 and d < size) {
            return true;
        }
    }
    return false;
}

template <class T> double Group<T>::mass() const {
    return std::accumulate(begin(), end(), 0.0, [](double sum, auto &i) { return sum + i.traits().mw; });
}

template <class T> std::vector<std::reference_wrapper<Point>> Group<T>::positions() const {
    return Faunus::member_view(begin(), end(), &Particle::pos);
    //return ranges::view::transform(*this, [](auto &i) -> Point& {return i.pos;});
}

template <class T> void Group<T>::wrap(Geometry::BoundaryFunction boundary) {
    boundary(cm);
    for (auto &i : *this)
        boundary(i.pos);
}

template <class T> void Group<T>::rotate(const Eigen::Quaterniond &Q, Geometry::BoundaryFunction boundary) {
    Geometry::rotate(begin(), end(), Q, boundary, -cm);
}

template <class T> void Group<T>::translate(const Point &d, Geometry::BoundaryFunction boundary) {
    cm += d;
    boundary(cm);
    for (auto &i : *this) {
        i.pos += d;
        boundary(i.pos);
    }
}

/**
 * @param mask Bitmask based on enum `Group::Selectors`
 * @return Lambda function that returns true if group matches mask
 *
 * ~~~ cpp
 * auto filter = Group::getSelectionFilter(Select::ACTIVE | Select::NEUTRAL);
 * bool b = filter(mygroup); // true if mygroup is active and uncharged
 * ~~~
 */
// constexpr std::function<bool(const Group<Particle> &)> getGroupFilter(unsigned int mask) {
//    return [mask = mask](const Group<Particle> &g) { return g.match<mask>(); };
//}

template class Group<Particle>;

void to_json(json &j, const Group<Particle> &g) {
    j = {
        {"id", g.id}, {"cm", g.cm}, {"atomic", g.atomic}, {"compressible", g.compressible},  {"size", g.size()}
    };
    if (g.capacity()>g.size())
        j["capacity"] = g.capacity();
    if (g.confid!=0)
        j["confid"] = g.confid;
}

void from_json(const json &j, Group<Particle> &g) {
    g.resize( j.at("size").get<int>() );
    g.trueend() = g.begin() + j.value("capacity", g.size());
    g.id = j.at("id").get<unsigned int>();
    g.cm = j.at("cm").get<Point>();
    g.atomic = j.at("atomic").template get<bool>();
    g.compressible = j.value("compressible", false);
    g.confid = j.value("confid", 0);
}

} // namespace Faunus
