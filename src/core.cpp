#define _USE_MATH_DEFINES
#include <doctest/doctest.h>
#include "core.h"
#include "units.h"
#include "aux/iteratorsupport.h"
#include "random.h"
#include "particle.h"
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/null_sink.h>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/numeric/accumulate.hpp>

namespace Faunus {

double roundValue(double value, const int number_of_digits) {
    std::stringstream o;
    o << std::setprecision(number_of_digits) << value;
    return std::stod(o.str());
}

void roundJSON(json& j, int number_of_digits) {
    if (j.is_object()) {
        for (auto& value : j) {
            if (value.is_number_float()) {
                value = roundValue(value, number_of_digits);
            } else if (value.is_object() && !value.empty()) {
                roundJSON(value, number_of_digits);
            }
        }
    }
}

double getValueInfinity(const json& j, const std::string& key) {
    auto value = j.find(key);
    if (value == j.end()) {
        throw std::runtime_error("unknown json key '" + key + "'");
    }
    if (value->is_string()) {
        if (*value == "inf" || *value == "oo") {
            return std::numeric_limits<double>::infinity();
        }
        if (*value == "-inf" || *value == "-oo") {
            return -std::numeric_limits<double>::infinity();
        }
        throw std::runtime_error("value must be number or 'inf'");
    }
    return static_cast<double>(*value);
}

/**
 * @brief Reads content of a JSON file into a JSON object
 * @param filename
 * @return json object from the file
 * @throw IOError when file cannot be read
 * @throw json::parse_error when json ill-formatted
 */
json loadJSON(const std::string& filename) {
    if (std::ifstream stream(filename); stream) {
        json j;
        stream >> j;
        return j;
    }
    throw IOError("cannot open file '{}'", filename);
}

TipFromTheManual::TipFromTheManual() : random(std::make_shared<Random>()) { random->seed(); }

/**
 * @brief Load JSON tips database
 * @param files vector of file names
 *
 * Iterate over `files` until file is successfully opened. The JSON
 * file must consist of objects with `key`s and markdown formatted
 * tips. If no file can be opened, the database is left empty.
 */
void TipFromTheManual::load(const std::vector<std::string>& files) {
    for (const auto& file : files) {
        try {
            if (database = loadJSON(file); !database.empty()) {
                break;
            }
        } catch (...) {}
    }
}

/**
 * @brief If possible, give help based on short keys/tags
 */
std::string TipFromTheManual::operator[](std::string_view key) {
    std::string tip;
    if (not tip_already_given) {
        // look for help for the given `key`
        if (auto it = database.find(key); it != database.end()) {
            tip = "\nNeed help, my young apprentice?\n\n" + it->get<std::string>();
            if (key == "coulomb") { // for the Coulomb potential, add additional table w. types
                tip += "\n" + database.at("coulomb types").get<std::string>();
            } else if (key == "custom") { // for the custom potential, add also list of symbols
                tip += "\n" + database.at("symbol").get<std::string>();
            }
            tip_already_given = true;
            if (asciiart) { // add ascii art?
                if (it = database.find("ascii"); it != database.end() && !it->empty() && it->is_array()) {
                    tip += random->sample(it->begin(), it->end())->get<std::string>() + "\n";
                }
            }
        }
        buffer = tip;
    }
    return (quiet) ? std::string() : tip;
}

void TipFromTheManual::pick(const std::string& key) { operator[](key); }

TipFromTheManual usageTip; // Global instance

// global loggers as a dummy instance
// they should be replaced with proper instances in faunus, pyfaunus and unittests if desired
std::shared_ptr<spdlog::logger> faunus_logger = spdlog::create<spdlog::sinks::null_sink_st>("null");
std::shared_ptr<spdlog::logger> mcloop_logger = faunus_logger;

std::string addGrowingSuffix(const std::string& file) {
    // using std::experimental::filesystem; // exp. c++17 feature, not available on MacOS (Dec. 2018)
    int cnt = 0;
    std::string newfile;
    auto exists = [&]() { return std::ifstream(newfile).good(); };
    do {
        newfile = fmt::format("{}.{}", file, cnt++);
    } while (not exists());
    return newfile;
}

std::tuple<const std::string&, const json&> jsonSingleItem(const json& j) {
    if (j.is_object() && j.size() == 1) {
        const auto it = j.cbegin();
        return {it.key(), it.value()};
    }
    throw std::runtime_error("invalid data: single key expected");
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
    auto h = xyz.dot(dir);
    Point xy = xyz - dir * h;
    auto x = xy.dot(dir2);
    Point y = xy - dir2 * x;
    auto theta = std::atan2(y.norm(), x);
    auto radius = xy.norm();
    return {radius, theta, h};
}

Point xyz2rtp(const Point& p, const Point& origin) {
    Point xyz = p - origin;
    const auto radius = xyz.norm();
    return {radius, std::atan2(xyz.y(), xyz.x()), std::acos(xyz.z() / radius)};
}

Point rtp2xyz(const Point& rtp, const Point& origin) {
    return origin + rtp.x() * Point(std::cos(rtp.y()) * std::sin(rtp.z()), std::sin(rtp.y()) * std::sin(rtp.z()),
                                    std::cos(rtp.z()));
}

/**
 * Generates random points on a cube with side-lengths 1 and origin at [0,0,0].
 * Keep going until inside embedded sphere of radius 0.5, then return normalized vector.
 */
Point randomUnitVector(Random& rand, const Point& directions) {
    constexpr double squared_radius = 0.25; // 0.5^2
    double squared_norm;
    Point position;
    do {
        for (int i = 0; i < 3; i++) {
            position[i] = (rand() - 0.5) * directions[i];
        }
        squared_norm = position.squaredNorm();
    } while (squared_norm > squared_radius);
    return position / std::sqrt(squared_norm);
}

TEST_CASE("[Faunus] randomUnitVector") {
    Random r;
    int n = 4e5;
    Point rtp(0, 0, 0);
    for (int i = 0; i < n; i++) {
        rtp += xyz2rtp(randomUnitVector(r));
    }
    rtp = rtp / n;
    CHECK(rtp.x() == doctest::Approx(1));
    CHECK(rtp.y() == doctest::Approx(0).epsilon(0.005));          // theta [-pi:pi] --> <theta>=0
    CHECK(rtp.z() == doctest::Approx(M_PI / 2.0).epsilon(0.005)); // phi [0:pi] --> <phi>=pi/2
}

Point randomUnitVectorPolar(Random& rand) { return rtp2xyz({1.0, 2.0 * M_PI * rand(), std::acos(2.0 * rand() - 1.0)}); }

TEST_CASE("[Faunus] randomUnitVectorPolar") {
    Random r;
    int n = 2e5;
    Point rtp(0, 0, 0);
    for (int i = 0; i < n; i++) {
        rtp += xyz2rtp(randomUnitVectorPolar(r));
    }
    rtp = rtp / n;
    CHECK(rtp.x() == doctest::Approx(1));
    CHECK(rtp.y() == doctest::Approx(0).epsilon(0.005));          // theta [-pi:pi] --> <theta>=0
    CHECK(rtp.z() == doctest::Approx(0.5 * M_PI).epsilon(0.005)); // phi [0:pi] --> <phi>=pi/2
}

GenericError::GenericError(const std::exception& e) : GenericError(e.what()) {}
GenericError::GenericError(const std::runtime_error& e) : std::runtime_error(e) {}
GenericError::GenericError(const std::string& msg) : std::runtime_error(msg) {}
GenericError::GenericError(const char* msg) : std::runtime_error(msg) {}

const json& ConfigurationError::attachedJson() const { return attached_json; }

ConfigurationError& ConfigurationError::attachJson(const json& j) {
    attached_json = j;
    return *this;
}

void displayError(spdlog::logger& logger, const std::exception& e, int level) {
    const std::string padding = level > 0 ? "... " : "";
    logger.error(padding + e.what());

    // ConfigurationError can carry a JSON snippet which should be shown for debugging.
    if (const auto* config_error = dynamic_cast<const ConfigurationError*>(&e);
        config_error != nullptr && !config_error->attachedJson().empty()) {
        if (level > 0) {
            logger.debug("... JSON snippet:\n{}", config_error->attachedJson().dump(4));
        } else {
            logger.debug("JSON snippet:\n{}", config_error->attachedJson().dump(4));
        }
    }

    // Process nested exceptions in a tail recursion.
    try {
        std::rethrow_if_nested(e);
    } catch (const std::exception& e) { displayError(logger, e, level + 1); }
}

TEST_SUITE_BEGIN("Core");

TEST_CASE("[Faunus] infinite/nan") {
    CHECK(std::isnan(0.0 / 0.0));
    CHECK(std::isnan(0.0 / 0.0 * 1.0));
    CHECK(std::isinf(std::numeric_limits<double>::infinity()));
    CHECK(std::isinf(std::log(0.0)));
    CHECK(std::isnan(std::sqrt(-1.0)));
    CHECK(std::numeric_limits<double>::has_signaling_NaN);
    CHECK(std::isnan(std::numeric_limits<double>::signaling_NaN()));
    CHECK(std::isnan(std::numeric_limits<double>::signaling_NaN() * 1.0));
}

TEST_CASE("[Faunus] distance") {
    std::vector<long long int> v = {10, 20, 30, 40, 30};
    auto rng = v | ranges::cpp20::views::filter([](auto i) { return i == 30; });
    CHECK(Faunus::distance(v.begin(), rng.begin()) == 2);
    auto it = rng.begin();
    CHECK(Faunus::distance(v.begin(), ++it) == 4);
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

/**
 * @param molarity Molar salt concentration
 * @param valencies valencies for participating ions {1,-1} ~ NaCl, {2,-1} ~ MgCl2, {1,3,-2} ~ KAl(SO4)2
 * @throw If stoichiometry cannot be resolved
 */
Electrolyte::Electrolyte(const double molarity, const std::vector<int>& valencies)
    : molarity(molarity)
    , valencies(valencies) {
    namespace rv = ranges::cpp20::views;
    const auto sum_positive = ranges::accumulate(valencies | rv::filter([](auto v) { return v > 0; }), 0);
    const auto sum_negative = ranges::accumulate(valencies | rv::filter([](auto v) { return v < 0; }), 0);
    const auto gcd = std::gcd(sum_positive, sum_negative);
    if (sum_positive == 0 || sum_negative == 0 || gcd == 0) {
        throw std::runtime_error("cannot resolve stoichiometry; did you provide both + and - ions?");
    }
    auto nu_times_squared_valency = valencies | rv::transform([&](const auto valency) {
                                        const auto nu = (valency > 0 ? -sum_negative : sum_positive) / gcd;
                                        return nu * valency * valency;
                                    });
    ionic_strength = 0.5 * molarity * static_cast<double>(ranges::accumulate(nu_times_squared_valency, 0));
    faunus_logger->debug("salt molarity {} --> ionic strength = {:.3f} mol/l", molarity, ionic_strength);
}

/**
 * Back calculates the molarity, assuming 1:-1 salt charges. Note that the all stored properties
 * are temperature dependent which is why the bjerrum_length is required.
 */
Electrolyte::Electrolyte(double debye_length, double bjerrum_length) {
    valencies = {1, -1};
    ionic_strength = molarity =
        std::pow(1.0 / debye_length, 2) / (8.0 * pc::pi * bjerrum_length * 1.0_angstrom * 1.0_molar);
    faunus_logger->debug("debyelength {} Å --> 1:1 salt molarity = {:.3f}", debye_length, molarity);
}

/**
 * The salt composition is automatically resolved, and the ionic strength, I, is calculated according to
 * I = 0.5 * concentration * sum( nu_i * valency_i^2 )
 * where `nu` are the minimum stoichiometric coefficients, deduced by assuming a net-neutral salt.
 */
double Electrolyte::ionicStrength() const { return ionic_strength; }

double Electrolyte::getMolarity() const { return molarity; }

const std::vector<int>& Electrolyte::getValencies() const { return valencies; }

/**
 * @param bjerrum_length Bjerrum length in Angstrom
 * @return Debye screening length in Angstrom
 */
double Electrolyte::debyeLength(const double bjerrum_length) const {
    return 1.0 / std::sqrt(8.0 * pc::pi * bjerrum_length * 1.0_angstrom * ionic_strength * 1.0_molar);
}

/**
 * Create a salt object using either:
 *
 * 1. `debyelength` and `epsr` pair (throws if the latter is absent)
 * 2. `molarity` and `valencies` pair (the latter has default value of [1,-1])
 *
 * ~~~ yaml
 *     debyelength: 30, epsr: 80
 *     molarity: 0.1, valencies: [1:-2]
 *     salt: 0.1, valencies: [1:-2]    # ok, but deprecated
 *     debyelength: 30                 # throws due to missing 'epsr'
 *     molarity: 0.1, valencies: [1:2] # throws due to violated electroneutrality
 * ~~~
 *
 * @throw if `debyelength` is found but no dielectric constant, `epsr`.
 */
std::optional<Electrolyte> makeElectrolyte(const json& j) {
    if (auto it = j.find("debyelength"); it != j.end()) {
        const auto debye_length = it->get<double>() * 1.0_angstrom;
        const auto relative_dielectric_constant = j.at("epsr").get<double>(); // may throw!
        const auto bjerrum_length = pc::bjerrumLength(relative_dielectric_constant);
        return Electrolyte(debye_length, bjerrum_length);
    }
    auto molarity = j.value("molarity", j.value("salt", 0.0));
    if (molarity > 0.0) {
        auto valencies = j.value("valencies", std::vector<int>{1, -1});
        return Electrolyte(molarity, valencies);
    }
    return std::nullopt;
}

TEST_CASE("[Faunus] Electrolyte") {
    using doctest::Approx;
    CHECK(Electrolyte(0.1, {1, -1}).ionicStrength() == Approx(0.1));                       // NaCl
    CHECK(Electrolyte(0.1, {2, -2}).ionicStrength() == Approx(0.5 * (0.1 * 4 + 0.1 * 4))); // CaSO₄
    CHECK(Electrolyte(0.1, {2, -1}).ionicStrength() == Approx(0.5 * (0.1 * 4 + 0.2)));     // CaCl₂
    CHECK(Electrolyte(0.1, {1, -2}).ionicStrength() == Approx(0.5 * (0.2 + 0.1 * 4)));     // K₂SO₄
    CHECK(Electrolyte(0.1, {1, -3}).ionicStrength() == Approx(0.5 * (0.3 + 0.1 * 9)));     // Na₃Cit
    CHECK(Electrolyte(0.1, {3, -1}).ionicStrength() == Approx(0.5 * (0.3 + 0.1 * 9)));     // LaCl₃
    CHECK(Electrolyte(0.1, {2, -3}).ionicStrength() == Approx(0.5 * (0.3 * 4 + 0.2 * 9))); // Ca₃(PO₄)₂
    CHECK(Electrolyte(0.1, {1, 3, -2}).ionicStrength() == Approx(0.5 * (0.1 * 1 + 0.1 * 9 + 0.1 * 2 * 4))); // KAl(SO₄)₂
    CHECK(Electrolyte(1.0, {2, 3, -2}).ionicStrength() == Approx(0.5 * (2 * 4 + 2 * 9 + 5 * 4))); // Ca₂Al₂(SO₄)₅
    CHECK_THROWS(Electrolyte(0.1, {1, 1}));
    CHECK_THROWS(Electrolyte(0.1, {-1, -1}));
    CHECK_THROWS(Electrolyte(0.1, {0, 0}));

    SUBCASE("debyeLength") {
        CHECK(Electrolyte(0.03, {1, -1}).debyeLength(7.0) == Approx(17.7376102214));
        CHECK_THROWS(makeElectrolyte(R"({"debyelength": 30.0"})"_json)); // 'epsr' is missing
        const auto bjerrum_length = pc::bjerrumLength(80);
        CHECK(makeElectrolyte(R"({"debyelength": 30.0, "epsr": 80})"_json).value().debyeLength(bjerrum_length) ==
              Approx(30));
    }

    SUBCASE("debye length input") {
        const auto bjerrum_length = 7.0;
        CHECK(Electrolyte(30, bjerrum_length).debyeLength(bjerrum_length) == Approx(30));
        CHECK(Electrolyte(30, bjerrum_length).getMolarity() == Approx(0.0104874272));
        CHECK(Electrolyte(30, bjerrum_length).ionicStrength() == Approx(0.0104874272));
    }
}

void to_json(json& j, const Electrolyte& electrolyte) {
    const auto bjerrum_length = pc::bjerrumLength(pc::T());
    j = {{"molarity", electrolyte.getMolarity()},
         {"valencies", electrolyte.getValencies()},
         {"molar ionic strength", electrolyte.ionicStrength()},
         {"debyelength", electrolyte.debyeLength(bjerrum_length)}};
}

} // namespace Faunus

template class nlohmann::basic_json<>;
