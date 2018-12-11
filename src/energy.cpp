
#include "energy.h"

void Faunus::Energy::Energybase::to_json(json &j) const {}

void Faunus::Energy::Energybase::sync(Energybase *, Change &) {}

void Faunus::Energy::Energybase::init() {}

void Faunus::Energy::to_json(json &j, const Energybase &base) {
    assert(not base.name.empty());
    if (base.timer) 
        j["relative time"] = base.timer.result();
    if (not base.cite.empty())
        j[base.name]["reference"] = base.cite;
    base.to_json( j[base.name] );
}

void Faunus::Energy::EwaldData::update(const Point &box) {
    L = box;
    int kcc = std::ceil(kc);
    check_k2_zero = 0.1*std::pow(2*pc::pi/L.maxCoeff(), 2);
    int kVectorsLength = (2*kcc+1) * (2*kcc+1) * (2*kcc+1) - 1;
    if (kVectorsLength == 0) {
        kVectors.resize(3,1);
        Aks.resize(1);
        kVectors.col(0) = Point(1,0,0); // Just so it is not the zero-vector
        Aks[0] = 0;
        kVectorsInUse = 1;
        Qion.resize(1);
        Qdip.resize(1);
    } else {
        double kc2 = kc*kc;
        kVectors.resize(3, kVectorsLength);
        Aks.resize(kVectorsLength);
        kVectorsInUse = 0;
        kVectors.setZero();
        Aks.setZero();
        int startValue = 1 - int(ipbc);
        for (int kx = 0; kx <= kcc; kx++) {
            double dkx2 = double(kx*kx);
            for (int ky = -kcc*startValue; ky <= kcc; ky++) {
                double dky2 = double(ky*ky);
                for (int kz = -kcc*startValue; kz <= kcc; kz++) {
                    double factor = 1.0;
                    if(kx > 0)
                        factor *= 2;
                    if(ky > 0 && ipbc)
                        factor *= 2;
                    if(kz > 0 && ipbc)
                        factor *= 2;
                    double dkz2 = double(kz*kz);
                    Point kv = 2*pc::pi*Point(kx/L.x(),ky/L.y(),kz/L.z());
                    double k2 = kv.dot(kv);
                    if (k2 < check_k2_zero) // Check if k2 != 0
                        continue;
                    if (spherical_sum)
                        if( (dkx2/kc2) + (dky2/kc2) + (dkz2/kc2) > 1)
                            continue;
                    kVectors.col(kVectorsInUse) = kv;
                    Aks[kVectorsInUse] = factor*std::exp(-k2/(4*alpha*alpha))/k2;
                    kVectorsInUse++;
                }
            }
        }
        Qion.resize(kVectorsInUse);
        Qdip.resize(kVectorsInUse);
        Aks.conservativeResize(kVectorsInUse);
        kVectors.conservativeResize(3,kVectorsInUse);
    }
}

void Faunus::Energy::from_json(const json &j, EwaldData &d) {
    d.alpha = j.at("alpha");
    d.rc = j.at("cutoff");
    d.kc = j.at("kcutoff");
    d.ipbc = j.value("ipbc", false);
    d.spherical_sum = j.value("spherical_sum", true);
    d.lB = pc::lB( j.at("epsr") );
    d.eps_surf = j.value("epss", 0.0);
    d.const_inf = (d.eps_surf < 1) ? 0 : 1; // if unphysical (<1) use epsr infinity for surrounding medium
}

void Faunus::Energy::to_json(json &j, const EwaldData &d) {
    j = {{"lB", d.lB}, {"ipbc", d.ipbc}, {"epss", d.eps_surf},
        {"alpha", d.alpha}, {"cutoff", d.rc}, {"kcutoff", d.kc},
        {"wavefunctions", d.kVectors.cols()}, {"spherical_sum", d.spherical_sum}};
}

double Faunus::Energy::Example2D::energy(Change &change) {
    double s=1+std::sin(2*pc::pi*i.x())+std::cos(2*pc::pi*i.y());
    if (i.x()>=-2.00 && i.x()<=-1.25) return 1*s;
    if (i.x()>=-1.25 && i.x()<=-0.25) return 2*s;
    if (i.x()>=-0.25 && i.x()<= 0.75) return 3*s;
    if (i.x()>= 0.75 && i.x()<= 1.75) return 4*s;
    if (i.x()>= 1.75 && i.x()<= 2.00) return 5*s;
    return 1e10;
}
