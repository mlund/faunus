#include <doctest/doctest.h>
#include "core.h"
#include "aux/iteratorsupport.h"
#include "random.h"
#include "units.h"
#include "particle.h"
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/null_sink.h>
#include <range/v3/view/filter.hpp>

namespace Faunus {

using doctest::Approx;

double _round(double x, int n) {
    std::stringstream o;
    o << std::setprecision(n) << x;
    return std::stod(o.str());
}

void _roundjson(json& j, int n) {
    if (j.is_object()) {
        for (auto& i : j) {
            if (i.is_number_float()) {
                i = _round(i, n);
            }
        }
    }
}

double value_inf(const json& j, const std::string& key) {
    auto it = j.find(key);
    if (it == j.end())
        throw std::runtime_error("unknown json key '" + key + "'");
    else if (it->is_string()) {
        if (*it == "inf")
            return std::numeric_limits<double>::infinity();
        if (*it == "-inf")
            return -std::numeric_limits<double>::infinity();
        throw std::runtime_error("value must be number or 'inf'");
    }
    return double(*it);
}

/**
 * @brief Reads content of a JSON file into a JSON object
 * @param filename
 * @return json object from the file
 * @throw IOError when file cannot be read
 * @throw json::parse_error when json ill-formatted
 */
json openjson(const std::string& filename) {
    json j;
    if (std::ifstream stream(filename); stream) {
        stream >> j;
    } else {
        throw IOError("cannot open file '{}'", filename);
    }
    return j;
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
void TipFromTheManual::load(const std::vector<std::string>& files) {
    // try loading `files`; stop if not empty
    for (const auto& file : files) {
        try {
            database = openjson(file); // allow for file not found
            if (!database.empty()) {
                break;
            }
        } catch (...) {
            // ignore
        }
    }
}

/**
 * @brief If possible, give help based on short keys/tags
 */
std::string TipFromTheManual::operator[](const std::string& key) {
    std::string t;
    if (not tip_already_given) {
        // look for help for the given `key`
        auto it = database.find(key);
        if (it != database.end()) {
            t = "\nNeed help, my young apprentice?\n\n" + it->get<std::string>();

            // for the Coulomb potential, add additional table w. types
            if (key == "coulomb") {
                t += "\n" + database.at("coulomb types").get<std::string>();
            }

            // for the custom potential, add also list of symbols
            if (key == "custom") {
                t += "\n" + database.at("symbol").get<std::string>();
            }

            tip_already_given = true;

            // add ascii art
            if (asciiart) {
                it = database.find("ascii");
                if (it != database.end()) {
                    if (not it->empty() and it->is_array()) {
                        t += random->sample(it->begin(), it->end())->get<std::string>() + "\n";
                    }
                }
            }
        }
        buffer = t;
    }
    return (quiet) ? std::string() : t;
}

void TipFromTheManual::pick(const std::string& key) { operator[](key); }

TipFromTheManual usageTip; // Global instance

// global loggers as a dummy instance
// they should be replaced with proper instances in faunus, pyfaunus and unittests if desired
std::shared_ptr<spdlog::logger> faunus_logger = spdlog::create<spdlog::sinks::null_sink_st>("null");
std::shared_ptr<spdlog::logger> mcloop_logger = faunus_logger;

std::string addGrowingSuffix(const std::string& file) {
    // using std::experimental::filesystem; // exp. c++17 feature, not available on MacOS (Dec. 2018)
    size_t cnt = 0;
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

std::tuple<const std::string&, const json&> jsonSingleItem(const json& j) {
    if (!j.is_object() || j.size() != 1) {
        throw std::runtime_error("invalid data: single key expected");
    }
    const auto j_it = j.cbegin();
    return {j_it.key(), j_it.value()};
}

json::size_type SingleUseJSON::count(const std::string& key) const { return json::count(key); }

bool SingleUseJSON::empty() const { return json::empty(); }

SingleUseJSON::SingleUseJSON(const json& j) : json(j) {}

std::string SingleUseJSON::dump(int w) const { return json::dump(w); }

void SingleUseJSON::clear() { json::clear(); }

json SingleUseJSON::at(const std::string& key) {
    json val = json::at(key);
    json::erase(key);
    return val;
}

json SingleUseJSON::operator[](const std::string& key) { return at(key); }

void SingleUseJSON::erase(const std::string& key) { json::erase(key); }
bool SingleUseJSON::is_object() const { return json::is_object(); }

Point xyz2rth(const Point& p, const Point& origin, const Point& dir, const Point& dir2) {
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

Point xyz2rtp(const Point& p, const Point& origin) {
    Point xyz = p - origin;
    double radius = xyz.norm();
    return {radius, std::atan2(xyz.y(), xyz.x()), std::acos(xyz.z() / radius)};
}

Point rtp2xyz(const Point& rtp, const Point& origin) {
    return origin + rtp.x() * Point(std::cos(rtp.y()) * std::sin(rtp.z()), std::sin(rtp.y()) * std::sin(rtp.z()),
                                    std::cos(rtp.z()));
}
Point ranunit(Random& rand, const Point& dir) {
    double r2;
    Point p;
    do {
        for (size_t i = 0; i < 3; i++) {
            p[i] = (rand() - 0.5) * dir[i];
        }
        r2 = p.squaredNorm();
    } while (r2 > 0.25);
    return p / std::sqrt(r2);
}

TEST_CASE("[Faunus] ranunit") {
    Random r;
    int n = 4e5;
    Point rtp(0, 0, 0);
    for (int i = 0; i < n; i++) {
        rtp += xyz2rtp(ranunit(r));
    }
    rtp = rtp / n;
    CHECK(rtp.x() == doctest::Approx(1));
    CHECK(rtp.y() == doctest::Approx(0).epsilon(0.005));          // theta [-pi:pi] --> <theta>=0
    CHECK(rtp.z() == doctest::Approx(pc::pi / 2).epsilon(0.005)); // phi [0:pi] --> <phi>=pi/2
}

Point ranunit_polar(Random& rand) { return rtp2xyz({1, 2 * pc::pi * rand(), std::acos(2 * rand() - 1)}); }

TEST_CASE("[Faunus] ranunit_polar") {
    Random r;
    int n = 2e5;
    Point rtp(0, 0, 0);
    for (int i = 0; i < n; i++) {
        rtp += xyz2rtp(ranunit_polar(r));
    }
    rtp = rtp / n;
    CHECK(rtp.x() == doctest::Approx(1));
    CHECK(rtp.y() == doctest::Approx(0).epsilon(0.005));          // theta [-pi:pi] --> <theta>=0
    CHECK(rtp.z() == doctest::Approx(pc::pi / 2).epsilon(0.005)); // phi [0:pi] --> <phi>=pi/2
}

GenericError::GenericError(const std::exception& e) : GenericError(e.what()) {}
GenericError::GenericError(const std::runtime_error& e) : std::runtime_error(e) {}
GenericError::GenericError(const std::string& msg) : std::runtime_error(msg) {}
GenericError::GenericError(const char* msg) : std::runtime_error(msg) {}

const json& ConfigurationError::attachedJson() const { return attached_json; }

ConfigurationError& ConfigurationError::attachJson(const json j) {
    attached_json = j;
    return *this;
}

void displayError(spdlog::logger& logger, const std::exception& e, int level) {
    const std::string padding = level > 0 ? "... " : "";

    logger.error(padding + e.what());

    // ConfigurationError can carry a JSON snippet which should be shown for debugging.
    if (auto config_error = dynamic_cast<const ConfigurationError*>(&e);
        config_error != nullptr && !config_error->attachedJson().empty()) {
        logger.debug(padding + "JSON snippet:\n{}", config_error->attachedJson().dump(4));
    }

    // Process nested exceptions in a tail recursion.
    try {
        std::rethrow_if_nested(e);
    } catch (const std::exception& e) {
        displayError(logger, e, level + 1);
    }
}

TEST_SUITE_BEGIN("Core");

TEST_CASE("[Faunus] infinite/nan") {
    CHECK(std::isnan(0.0 / 0.0));
    CHECK(std::isnan(0.0 / 0.0 * 1.0));
    CHECK(std::isinf(pc::infty));
    CHECK(std::isinf(std::log(0.0)));
    CHECK(std::isnan(std::sqrt(-1.0)));
    CHECK(std::numeric_limits<double>::has_signaling_NaN);
    CHECK(std::isnan(std::numeric_limits<double>::signaling_NaN()));
    CHECK(std::isnan(std::numeric_limits<double>::signaling_NaN() * 1.0));
}

TEST_CASE("[Faunus] distance") {
    std::vector<long long int> v = {10, 20, 30, 40, 30};
    auto rng = v | ranges::cpp20::views::filter([](int i) { return i == 30; });
    CHECK(Faunus::distance(v.begin(), rng.begin()) == 2);
    auto it = rng.begin();
    CHECK(Faunus::distance(v.begin(), ++it) == 4);
}

TEST_CASE("[Faunus] member_view") {
    struct data {
        int N = 0;
    };
    std::vector<data> v(2);                                // original vector
    auto Nvec = member_view(v.begin(), v.end(), &data::N); // ref. to all N's in vec
    int cnt = 1;
    for (int& i : Nvec) { // modify N values
        i = cnt++;
    }
    CHECK(v[0].N == 1); // original vector was changed
    CHECK(v[1].N == 2); // original vector was changed
}

TEST_CASE("[Faunus] asEigenMatrix") {
    using doctest::Approx;
    std::vector<Particle> v(4);
    v[0].pos.x() = 5;
    v[1].pos.y() = 10;
    v[2].pos.z() = 2;
    auto m = asEigenMatrix(v.begin(), v.end(), &Particle::pos);

    CHECK(m.cols() == 3);
    CHECK(m.rows() == 4);
    CHECK(m.row(0).x() == 5);
    CHECK(m.row(1).y() == 10);
    CHECK(m.row(2).z() == 2);
    CHECK(m.sum() == 17);
    m.row(0).z() += 0.5;
    CHECK(v[0].pos.z() == Approx(0.5));

    v[2].charge = 2;
    v[3].charge = -12;
    auto m2 = asEigenVector(v.begin() + 1, v.end(), &Particle::charge);
    CHECK(m2.cols() == 1);
    CHECK(m2.rows() == 3);
    CHECK(m2.col(0).sum() == Approx(-10));
}

TEST_SUITE_END();

} // namespace Faunus

template class nlohmann::basic_json<>;
