#include "core.h"
#include "random.h"
#include "units.h"
#include <iomanip>
#include <fstream>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/null_sink.h>

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

    TipFromTheManual::TipFromTheManual() {
        random = std::make_shared<Random>();
        random->seed();
    }

    /**
     * @brief Load JSON tips database
     * @param files vector of file names
     *
     * Iterate over `files` until file is successfully opened. The JSON
     * file must consist of objects with `key`s and markdown formatted
     * tips. If no file can be opened, the database is left empty.
     */
    void TipFromTheManual::load(const std::vector<std::string> &files) {
        // try loading `files`; stop of not empty
        for (auto &i : files) {
            db = openjson(i, false); // allow for file not found
            if (not db.empty())
                break;
        }
    }

    /**
     * @brief If possible, give help based on short keys/tags
     */
    std::string TipFromTheManual::operator[](const std::string &key) {
        std::string t;
        if (not tip_already_given) {
            // look for help for the given `key`
            auto it = db.find(key);
            if (it!=db.end()) {
                t = "\nNeed help, my young apprentice?\n\n" + it->get<std::string>();

                // for the Coulomb potential, add additional table w. types
                if (key=="coulomb")
                    t += "\n" + db.at("coulomb types").get<std::string>();

                // for the custom potential, add also list of symbols
                if (key=="custom")
                    t += "\n" + db.at("symbol").get<std::string>();

                tip_already_given = true;

                // add ascii art
                if (asciiart) {
                    it = db.find("ascii");
                    if (it != db.end())
                        if (not it->empty() and it->is_array())
                            t += random->sample(it->begin(), it->end())->get<std::string>() + "\n";
                }
            }
            buffer = t;
        }
        return (quiet) ? std::string() : t;
    }

    TipFromTheManual usageTip; // Global instance

    // global loggers as a dummy instance
    // they should be replaced with proper instances in faunus, pyfaunus and unittests if desired
    std::shared_ptr<spdlog::logger> faunus_logger = spdlog::create<spdlog::sinks::null_sink_st>("null");
    std::shared_ptr<spdlog::logger> mcloop_logger = faunus_logger;

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

    json::size_type SingleUseJSON::count(const std::string &key) const { return json::count(key); }

    bool SingleUseJSON::empty() const { return json::empty(); }

    SingleUseJSON::SingleUseJSON(const json &j) : json(j) {}

    std::string SingleUseJSON::dump(int w) const { return json::dump(w); }

    void SingleUseJSON::clear() { json::clear(); }

    json SingleUseJSON::at(const std::string &key) {
        json val = json::at(key);
        json::erase(key);
        return val;
    }

    json SingleUseJSON::operator[](const std::string &key) { return at(key); }

    void SingleUseJSON::erase(const std::string &key) { json::erase(key); }
    bool SingleUseJSON::is_object() const { return json::is_object(); }

    Point xyz2rth(const Point &p, const Point &origin, const Point &dir, const Point &dir2) {
        assert(fabs(dir.norm() - 1.0) < 1e-6);
        assert(fabs(dir2.norm() - 1.0) < 1e-6);
        assert(fabs(dir.dot(dir2)) < 1e-6); // check if unit-vectors are perpendicular
        Point xyz = p - origin;
        double h = xyz.dot(dir);
        Point xy = xyz - dir * h;
        double x = xy.dot(dir2);
        Point y = xy - dir2 * x;
        double theta = std::atan2(y.norm(), x);
        double radius = xy.norm();
        return {radius, theta, h};
    }

    Point xyz2rtp(const Point &p, const Point &origin) {
        Point xyz = p - origin;
        double radius = xyz.norm();
        return {radius, std::atan2(xyz.y(), xyz.x()), std::acos(xyz.z() / radius)};
    }

    Point rtp2xyz(const Point &rtp, const Point &origin) {
        return origin + rtp.x() * Point(std::cos(rtp.y()) * std::sin(rtp.z()), std::sin(rtp.y()) * std::sin(rtp.z()),
                                        std::cos(rtp.z()));
    }
    Point ranunit(Random &rand, const Point &dir) {
        double r2;
        Point p;
        do {
            for (size_t i = 0; i < 3; i++)
                p[i] = (rand() - 0.5) * dir[i];
            r2 = p.squaredNorm();
        } while (r2 > 0.25);
        return p / std::sqrt(r2);
    }

    Point ranunit_polar(Random &rand) { return rtp2xyz({1, 2 * pc::pi * rand(), std::acos(2 * rand() - 1)}); }

} // end of Faunus namespace

template class nlohmann::basic_json<>;

