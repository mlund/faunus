#include "core.h"

using namespace Faunus;
using namespace std;

int main() {

    Random r;

    vector<int> v = {1,2,3,4,5,6};

    auto i = v | ranges::view::sample(2, r.engine);
    cout << i << endl;
}