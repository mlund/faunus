#include <faunus/faunus.h>
using namespace Faunus;

int main(int argc, char** argv) {
  if (argc==4 || argc==5) {
    string file = argv[1];
    double prope = std::stod(argv[2]);
    int it = std::stoi(argv[3]);
    cout << "file       = " << file << "\n"
      << "prope      = " << prope << " aa\n"
      << "iterations = " << it << "\n";

    // PQR or AAM molecular file format?
    std::vector<PointParticle> p;
    if (file.find(".pqr")!=std::string::npos)
      FormatPQR::load(file,p);
    else
      FormatAAM::load(file,p);
    double V=Geometry::calcVolume(p,it,prope);

    cout << "volume     = " << V << " aa^3 = " << V/1e3 << " nm^3\n";
    cout << "radius     = " << cbrt(3*V/4/acos(-1)) << " aa\n";

    // if molecular weight is given
    if (argc==5) {
      double mw = std::stod(argv[4]);
      double rho = mw/6.022e23/V * 1e24;
      cout << "density    = " << rho << " g/ml\n";
      cout << "specific volume = " << 1/rho << " ml/g\n";
    }

  } else
    std::cerr << "usage: " << argv[0]
      << " file(.pqr|.aam) properadius iterations [Mw]\n"; 
}
