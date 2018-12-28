#include "core.h"

namespace Faunus {

    double _round(double x, int n) {
        std::stringstream o;
        o << std::setprecision(n) << x;
        return std::stod(o.str());
    }

    void _roundjson(json &j, int n) {
        if (j.is_object())
            for (auto &i : j)
                if (i.is_number_float())
                    i = _round(i,n);
    }

    double value_inf(const json &j, const std::string &key) {
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

    void to_json(nlohmann::json &j, const Tensor &t) {
        j = { t(0,0), t(0,1), t(0,2), t(1,1), t(1,2), t(2,2) };
    }

    void from_json(const nlohmann::json &j, Tensor &t) {
        if ( j.size()!=6 || !j.is_array() )
            throw std::runtime_error("Json->Tensor: array w. exactly six coefficients expected.");
        t = Tensor(j[0],j[1],j[2],j[3],j[4],j[5]);
    }

    json merge(const json &a, const json &b) {
        json result = a.flatten();
        json tmp = b.flatten();
        for ( auto it = tmp.begin(); it != tmp.end(); ++it )
            result[it.key()] = it.value();
        return result.unflatten();
    }

    json openjson(const std::string &file, bool throw_if_file_not_found) {
        json js;
        std::ifstream f(file);
        if ( f ) {
            try {
                f >> js;
            }
            catch(std::exception& e) {
                throw std::runtime_error("Syntax error in JSON file " + file + ": " + e.what());
            }
        }
        else if (throw_if_file_not_found)
            throw std::runtime_error("Cannot find or read JSON file " + file);
        return js;
    }

    bool assertKeys(const json &j, const std::vector<std::string> &okkeys, bool exception) {
        assert(j.is_object());
        for (auto it : j.items()) {
            if (std::find(okkeys.begin(), okkeys.end(), it.key()) == okkeys.end() ) {
                if (exception)
                    throw std::runtime_error("unknown key: "s + it.key());
                else return false;
            }
        }
        return true;
    }

    void TipFromTheManual::load(const std::vector<std::string> &files) {
        slump.seed();
        // try loading `files`; stop of not empty
        for (auto &i : files) {
            db = openjson(i, false); // allow for file not found
            if (not db.empty())
                break;
        }
    }

    std::string TipFromTheManual::operator[](const std::string &key) {
        std::string t;
        if (not tip_already_given) {
            auto it = db.find(key);
            if (it!=db.end()) {
                t = "\n\nNeed help, my young Padawan?\n\n" + it->get<std::string>();
                if (key=="coulomb")
                    t += db.at("coulomb types").get<std::string>();
                auto ascii = db["ascii"].get<std::vector<std::string>>();
                if (not ascii.empty()) {
                    t += *(slump.sample(ascii.begin(), ascii.end())) + "\n";
                    tip_already_given = true;
                }
            }
        }
        return t; // empty string of no tips available
    }

    TipFromTheManual usageTip; // Global instance

    Tensor::Tensor(double xx, double xy, double xz, double yy, double yz, double zz) {
        (*this) << xx, xy, xz, xy, yy, yz, xz, yz, zz;
    }

    void Tensor::rotate(const Tensor::base &m) {
        (*this) = m * (*this) * m.transpose();
    }

    void Tensor::eye() { *this = base::Identity(3, 3); }

    QuaternionRotate::QuaternionRotate() {}

    QuaternionRotate::QuaternionRotate(double angle, Point u) { set(angle,u); }

    void QuaternionRotate::set(double angle, Point u) {
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

    Point QuaternionRotate::operator()(Point a, std::function<void(Point &)> boundary,
            const Point &shift) const {
        a = a - shift;
        boundary(a);
        a = first * a + shift;
        boundary(a);
        return a;
        // https://www.cc.gatech.edu/classes/AY2015/cs4496_spring/Eigen.html
    }

    auto QuaternionRotate::operator()(const Eigen::Matrix3d &a) const {
        return second * a * second.transpose();
    }

    std::string addGrowingSuffix(const std::string& file) {
        //using std::experimental::filesystem; // exp. c++17 feature, not available on MacOS (Dec. 2018)
        size_t cnt=0;
        std::string newfile;
        auto exists = [&]() {
            std::ifstream f(newfile.c_str());
            return f.good();
        };

        do {
            newfile = file + "." + std::to_string(cnt++);
        } while (not exists());
        return newfile;
    }

} // end of namespace
