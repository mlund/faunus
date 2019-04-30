#include "space.h"

namespace Faunus {

bool Change::data::operator<(const Faunus::Change::data &a) const { return index < a.index; }

void Change::clear() {
    dV = false;
    all = false;
    dN = false;
    moved2moved = true;
    groups.clear();
    assert(empty());
}
bool Change::empty() const {
    if (dV==false)
        if (all==false)
            if (groups.empty())
                if (dN==false)
                    return true;
    return false;
}
Change::operator bool() const {
    return not empty();
}
} // namespace