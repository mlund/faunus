#include "core.h"

double Faunus::_round(double x, int n) {
    std::stringstream o;
    o << std::setprecision(n) << x;
    return std::stod(o.str());
}

void Faunus::_roundjson(Faunus::json &j, int n) {
    if (j.is_object())
        for (auto &i : j)
            if (i.is_number_float())
                i = _round(i,n);
}

double Faunus::value_inf(const Faunus::json &j, const std::string &key) {
    auto it = j.find(key);
    if (it==j.end())
        throw std::runtime_error("unknown json key '" + key + "'");
    else
    if (it->is_string()) {
        if (*it=="inf")
            return std::numeric_limits<double>::infinity();
        if (*it=="-inf")
            return -std::numeric_limits<double>::infinity();
        throw std::runtime_error("value must be number or 'inf'");
    }
    return double(*it);
}

void Faunus::to_json(nlohmann::json &j, const Faunus::Tensor &t) {
    j = { t(0,0), t(0,1), t(0,2), t(1,1), t(1,2), t(2,2) };
}

void Faunus::from_json(const nlohmann::json &j, Faunus::Tensor &t) {
    if ( j.size()!=6 || !j.is_array() )
        throw std::runtime_error("Json->Tensor: array w. exactly six coefficients expected.");
    t = Tensor(j[0],j[1],j[2],j[3],j[4],j[5]);
}

void Faunus::from_json(const Faunus::json &j, Faunus::Random &r) {
    if (j.is_object()) {
        auto seed = j.value("seed", std::string());
        try {
            if (seed=="default" || seed=="fixed")
                return;
            if (seed=="hardware")
                r.engine = decltype(r.engine)(std::random_device()());
            else if (!seed.empty()) {
                std::stringstream s(seed);
                s.exceptions( std::ios::badbit | std::ios::failbit );
                s >> r.engine;
            }
        }
        catch (std::exception &e) {
            std::cerr << "error initializing random from json: " << e.what();
            throw;
        }
    }
}

void Faunus::to_json(Faunus::json &j, const Faunus::Random &r) {
    std::ostringstream o;
    o << r.engine;
    j["seed"] = o.str();
}

Faunus::Point Faunus::ranunit_neuman(Faunus::Random &rand) {
    double r2;
    Point p;
    do {
        p = {rand()-0.5, rand()-0.5, rand()-0.5};
        r2 = p.squaredNorm();
    } while ( r2 > 0.25 );
    return p / std::sqrt(r2);
}

Faunus::Point Faunus::ranunit_polar(Faunus::Random &rand) {
    return rtp2xyz( {1, 2*pc::pi*rand(), std::acos(2*rand()-1)} );
}

Faunus::Point Faunus::xyz2rtp(const Faunus::Point &p, const Faunus::Point &origin) {
    Point xyz = p - origin;
    double radius = xyz.norm();
    return {
            radius,
            std::atan2( xyz.y(), xyz.x() ),
            std::acos( xyz.z()/radius) };
}

Faunus::Point Faunus::rtp2xyz(const Faunus::Point &rtp, const Faunus::Point &origin) {
    return origin + rtp.x() * Point(
            std::cos(rtp.y()) * std::sin(rtp.z()),
            std::sin(rtp.y()) * std::sin(rtp.z()),
            std::cos(rtp.z()) );
}

Faunus::json Faunus::merge(const Faunus::json &a, const Faunus::json &b) {
    json result = a.flatten();
    json tmp = b.flatten();
    for ( auto it = tmp.begin(); it != tmp.end(); ++it )
        result[it.key()] = it.value();
    return result.unflatten();
}

Faunus::json Faunus::openjson(const std::string &file) {
    json js;
    std::ifstream f (file );
    if ( f ) {
        try {
            f >> js;
        }
        catch(std::exception& e) {
            throw std::runtime_error("Syntax error in JSON file " + file + ": " + e.what());
        }
    }
    else
        throw std::runtime_error("Cannot find or read JSON file " + file);
    return js;
}

Faunus::Tensor::Tensor(double xx, double xy, double xz, double yy, double yz, double zz) {
    (*this) << xx, xy, xz, xy, yy, yz, xz, yz, zz;
}

void Faunus::Tensor::rotate(const Faunus::Tensor::base &m) {
    (*this) = m * (*this) * m.transpose();
}

void Faunus::Tensor::eye() { *this = base::Identity(3, 3); }

Faunus::QuaternionRotate::QuaternionRotate() {}

Faunus::QuaternionRotate::QuaternionRotate(double angle, Faunus::Point u) { set(angle,u); }

void Faunus::QuaternionRotate::set(double angle, Faunus::Point u) {
    this->angle = angle;
    u.normalize();
    first = Eigen::AngleAxisd(angle, u);
    second << 0, -u.z(), u.y(), u.z(), 0, -u.x(), -u.y(), u.x(), 0;
    second =
            Eigen::Matrix3d::Identity() + second * std::sin(angle)
            + second * second * (1 - std::cos(angle));

    // Quaternion can be converted to rotation matrix:
    // second = first.toRotationMatrix()
}

Faunus::Point Faunus::QuaternionRotate::operator()(Faunus::Point a, std::function<void(Faunus::Point &)> boundary,
                                                   const Faunus::Point &shift) const {
    a = a - shift;
    boundary(a);
    a = first * a + shift;
    boundary(a);
    return a;
    // https://www.cc.gatech.edu/classes/AY2015/cs4496_spring/Eigen.html
}

auto Faunus::QuaternionRotate::operator()(const Eigen::Matrix3d &a) const {
    return second * a * second.transpose();
}

void Faunus::Radius::to_json(Faunus::json &j) const { j["r"] = radius; }

void Faunus::Radius::from_json(const Faunus::json &j) { radius = j.value("r", 0.0); }

void Faunus::Charge::to_json(Faunus::json &j) const { j["q"] = charge; }

void Faunus::Charge::from_json(const Faunus::json &j) { charge = j.value("q", 0.0); }

void Faunus::Dipole::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &) {
    mu = q * mu;
}

void Faunus::Dipole::to_json(Faunus::json &j) const {
    j["mulen"] = mulen ;
    j["mu"] = mu;
}

void Faunus::Dipole::from_json(const Faunus::json &j) {
    mulen = j.value("mulen", 0.0);
    mu = j.value("mu", Point(1,0,0) );
}

void Faunus::Quadrupole::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m) { Q.rotate(m); }

void Faunus::Quadrupole::to_json(Faunus::json &j) const { j["Q"] = Q; }

void Faunus::Quadrupole::from_json(const Faunus::json &j) { Q = j.value("Q", Q); }

void Faunus::Cigar::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &) {
    scdir = q * scdir;
}

void Faunus::Cigar::to_json(Faunus::json &j) const {
    j["sclen"] = sclen;
    j["scdir"] = scdir;
}

void Faunus::Random::seed() { engine = std::mt19937(std::random_device()()); }

Faunus::Random::Random() : dist01(0,1) {}

double Faunus::Random::operator()() { return dist01(engine); }

int Faunus::Random::range(int min, int max) {
    std::uniform_int_distribution<int> d(min, max);
    return d(engine);
}

std::string Faunus::u8::bracket(const std::string &s) {
    return "\u27e8" + s + "\u27e9";
}

void Faunus::ParticlePropertyBase::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &) {}
