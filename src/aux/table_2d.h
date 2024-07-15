#pragma once
#include <map>
#include <vector>
#include <string>
#include <average.h>
#include <Eigen/Core>
#include <fstream>

namespace Faunus {
/**
 * @brief General class for handling 2D tables - xy data, for example.
 * @date Lund 2011
 * @note `Tx` is used as the `std::map` key and which may be
 * problematic due to direct floating point comparison (== operator).
 * We have not experienced any issues with this, though. This uses
 * `std::map` and table lookup is of complexity logarithmic with N.
 */
template <typename Tx, typename Ty> class Table2D
{
  protected:
    typedef std::map<Tx, Ty> Tmap;

    Ty count()
    {
        Ty cnt = 0;
        for (auto& m : map)
            cnt += m.second;
        return cnt;
    }

    Tx dx;
    Tmap map;
    std::string name;

  private:
    Tx round(Tx x) { return (x >= 0) ? int(x / dx + 0.5) * dx : int(x / dx - 0.5) * dx; }

    double get(Tx x) { return operator()(x); }

  public:
    enum type
    {
        HISTOGRAM,
        XYDATA
    };

    type tabletype;

    /** @brief Sum of all y values (same as `count()`) */
    Ty sumy() const
    {
        Ty sum = 0;
        for (auto& m : map)
            sum += m.second;
        return sum;
    }

    /**
     * @brief Constructor
     * @param resolution Resolution of the x axis
     * @param key Table type: HISTOGRAM or XYDATA
     */
    Table2D(Tx resolution = 0.2, type key = XYDATA)
    {
        tabletype = key;
        setResolution(resolution);
    }

    /** @brief Convert to map */
    std::map<std::string, std::vector<double>> to_map()
    {
        std::map<std::string, std::vector<double>> m;
        m["x"].reserve(map.size());
        m["y"].reserve(map.size());
        for (auto& i : map) {
            m["x"].push_back(i.first);
            m["y"].push_back(get(i.first));
        }
        return m;
    }

    void clear() { map.clear(); }

    void setResolution(Tx resolution)
    {
        assert(resolution > 0);
        dx = resolution;
        map.clear();
    }

    void setResolution(std::vector<Tx>& resolution)
    {
        assert(resolution[0] > 0);
        dx = resolution[0];
        map.clear();
    }

    /** @brief Access operator - returns reference to y(x) */
    Ty& operator()(Tx x) { return map[round(x)]; }

    /** @brief Access operator - returns reference to y(x) */
    Ty& operator()(std::vector<Tx>& x) { return map[round(x[0])]; }

    /** @brief Find key and return corresponding value otherwise zero*/
    Ty find(std::vector<Tx>& x)
    {
        Ty value = 0;
        auto it = map.find(round(x[0]));
        if (it != map.end())
            value = it->second;
        return value;
    }

    /** @brief Save table to disk */
    template <class T = double> void save(const std::string& filename, T scale = 1, T translate = 0)
    {
        if (tabletype == HISTOGRAM) {
            if (!map.empty())
                map.begin()->second *= 2; // compensate for half bin width
            if (map.size() > 1)
                (--map.end())->second *= 2; // -//-
        }

        if (!map.empty()) {
            std::ofstream f(filename.c_str());
            f.precision(10);
            if (f) {
                for (auto& m : map)
                    f << m.first << " " << (m.second + translate) * scale << "\n";
            }
        }

        if (tabletype == HISTOGRAM) {
            if (!map.empty())
                map.begin()->second /= 2; // restore half bin width
            if (map.size() > 1)
                (--map.end())->second /= 2; // -//-
        }
    }

    /** @brief Save normalized table to disk */
    template <class T = double> void normSave(const std::string& filename)
    {
        if (tabletype == HISTOGRAM) {
            if (!map.empty())
                map.begin()->second *= 2; // compensate for half bin width
            if (map.size() > 1)
                (--map.end())->second *= 2; // -//-
        }

        if (!map.empty()) {
            std::ofstream f(filename.c_str());
            f.precision(10);
            Ty cnt = count() * dx;
            if (f) {
                for (auto& m : map)
                    f << m.first << " " << m.second / cnt << "\n";
            }
        }

        if (tabletype == HISTOGRAM) {
            if (!map.empty())
                map.begin()->second /= 2; // restore half bin width
            if (map.size() > 1)
                (--map.end())->second /= 2; // -//-
        }
    }

    /** @brief Sums up all previous elements and saves table to disk */
    template <class T = double> void sumSave(std::string filename, T scale = 1)
    {
        if (tabletype == HISTOGRAM) {
            if (!map.empty())
                map.begin()->second *= 2; // compensate for half bin width
            if (map.size() > 1)
                (--map.end())->second *= 2; // -//-
        }

        if (!map.empty()) {
            std::ofstream f(filename.c_str());
            f.precision(10);
            if (f) {
                Ty sum_t = 0.0;
                for (auto& m : map) {
                    sum_t += m.second;
                    f << m.first << " " << sum_t * scale << "\n";
                }
            }
        }

        if (tabletype == HISTOGRAM) {
            if (!map.empty())
                map.begin()->second /= 2; // restore half bin width
            if (map.size() > 1)
                (--map.end())->second /= 2; // -//-
        }
    }

    const Tmap& getMap() const { return map; }

    Tmap& getMap() { return map; }

    Tx getResolution() { return dx; }

    /*! Returns average */
    Tx mean()
    {
        assert(!map.empty());
        Tx avg = 0;
        for (auto& m : map)
            avg += m.first * m.second;
        return avg / count();
    }

    /*! Returns standard deviation */
    Tx std()
    {
        assert(!map.empty());
        Tx std2 = 0;
        Tx avg = mean();
        for (auto& m : map)
            std2 += m.second * (m.first - avg) * (m.first - avg);
        return sqrt(std2 / count());
    }

    /*! Returns iterator of minumum y */
    typename Tmap::const_iterator min()
    {
        assert(!map.empty());
        Ty min = std::numeric_limits<Ty>::max();
        typename Tmap::const_iterator it;
        for (auto m = map.begin(); m != map.end(); ++m)
            if (m->second < min) {
                min = m->second;
                it = m;
            }
        return it;
    }

    /*! Returns iterator of maximum y */
    typename Tmap::const_iterator max()
    {
        assert(!map.empty());
        Ty max = std::numeric_limits<Ty>::min();
        typename Tmap::const_iterator it;
        for (auto m = map.begin(); m != map.end(); ++m)
            if (m->second > max) {
                max = m->second;
                it = m;
            }
        return it;
    }

    /*! Returns x at minumum x */
    Tx minx()
    {
        assert(!map.empty());
        Tx x = 0;
        for (auto& m : map) {
            x = m.first;
            break;
        }
        return x;
    }

    /*! Returns average in interval */
    Ty avg(const std::vector<Tx>& limits)
    {
        Ty avg = 0;
        int cnt = 0;
        assert(!map.empty());
        for (auto& m : map) {
            if (m.first >= limits[0] && m.first <= limits[1]) {
                avg += m.second;
                ++cnt;
            }
        }
        if (cnt > 0)
            avg /= cnt;
        return avg;
    }

    /**
     * @brief Convert table2D to vector of floats
     */
    std::vector<double> hist2buf(int& size)
    {
        std::vector<double> sendBuf;
        assert(!map.empty());
        for (auto& m : map) {
            sendBuf.push_back(m.first);
            sendBuf.push_back(m.second);
        }
        sendBuf.resize(size, -1);
        return sendBuf;
    }

    /**
     * @brief Convert vector of floats to table2D
     */
    void buf2hist(std::vector<double>& v)
    {
        this->clear();
        assert(!v.empty());
        std::map<double, Average<double>> all;
        for (int i = 0; i < int(v.size()) - 1; i += 2)
            if (v.at(i + 1) != -1)
                all[v.at(i)] += v.at(i + 1);
        for (auto& m : all)
            this->operator()(m.first) = m.second.avg();
    }

    /**
     * @brief Load table from disk
     * @note The first line - used for comments - is ignored.
     * @todo Implement end bin compensation as in the save()
     * function when loading HISTOGRAMs
     */
    bool load(const std::string& filename)
    {
        std::ifstream f(filename.c_str());
        if (f) {
            map.clear();
            while (!f.eof()) {
                Tx x;
                double y;
                f >> x >> y;
                operator()(x) = y;
            }
            if (tabletype == HISTOGRAM) {
                if (!map.empty())
                    map.begin()->second /= 2; // restore half bin width
                if (map.size() > 1)
                    (--map.end())->second /= 2; // -//-
            }
            return true;
        }
        return false;
    }

    /**
     * @brief Convert table to matrix
     */
    Eigen::MatrixXd tableToMatrix()
    {
        assert(!this->map.empty() && "Map is empty!");
        Eigen::MatrixXd table(2, map.size());
        table.setZero();
        int I = 0;
        for (auto& m : this->map) {
            table(0, I) = m.first;
            table(1, I) = m.second;
            I++;
        }
        return table;
    }
};

/**
 * @brief Subtract two tables
 */
template <class Tx, class Ty> Table2D<Tx, Ty> operator-(Table2D<Tx, Ty>& a, Table2D<Tx, Ty>& b)
{
    assert(a.tabletype == b.tabletype && "Table a and b needs to be of same type");
    Table2D<Tx, Ty> c(std::min(a.getResolution(), b.getResolution()), a.tabletype);
    auto a_map = a.getMap();
    auto b_map = b.getMap();

    if (a.tabletype == Table2D<Tx, Ty>::HISTOGRAM) {
        if (!a_map.empty())
            a_map.begin()->second *= 2; // compensate for half bin width
        if (a_map.size() > 1)
            (--a_map.end())->second *= 2; // -//-
        if (!b_map.empty())
            b_map.begin()->second *= 2; // compensate for half bin width
        if (b_map.size() > 1)
            (--b_map.end())->second *= 2; // -//-
    }

    for (auto& m1 : a_map) {
        for (auto& m2 : b_map) {
            c(m1.first) = m1.second - m2.second;
            break;
        }
    }

    if (a.tabletype == Table2D<Tx, Ty>::HISTOGRAM) {
        if (!a_map.empty())
            a_map.begin()->second /= 2; // compensate for half bin width
        if (a_map.size() > 1)
            (--a_map.end())->second /= 2; // -//-
        if (!b_map.empty())
            b_map.begin()->second /= 2; // compensate for half bin width
        if (b_map.size() > 1)
            (--b_map.end())->second /= 2; // -//-
    }
    return c;
}

/**
 * @brief Addition two tables
 */
template <class Tx, class Ty> Table2D<Tx, Ty> operator+(Table2D<Tx, Ty>& a, Table2D<Tx, Ty>& b)
{
    assert(a.tabletype == b.tabletype && "Table a and b needs to be of same type");
    Table2D<Tx, Ty> c(std::min(a.getResolution(), b.getResolution()), a.tabletype);
    auto a_map = a.getMap();
    auto b_map = b.getMap();

    if (a.tabletype == Table2D<Tx, Ty>::HISTOGRAM) {
        if (!a_map.empty())
            a_map.begin()->second *= 2; // compensate for half bin width
        if (a_map.size() > 1)
            (--a_map.end())->second *= 2; // -//-
        if (!b_map.empty())
            b_map.begin()->second *= 2; // compensate for half bin width
        if (b_map.size() > 1)
            (--b_map.end())->second *= 2; // -//-
    }

    for (auto& m : a_map) {
        c(m.first) += m.second;
    }
    for (auto& m : b_map) {
        c(m.first) += m.second;
    }

    if (a.tabletype == Table2D<Tx, Ty>::HISTOGRAM) {
        if (!a_map.empty())
            a_map.begin()->second /= 2; // compensate for half bin width
        if (a_map.size() > 1)
            (--a_map.end())->second /= 2; // -//-
        if (!b_map.empty())
            b_map.begin()->second /= 2; // compensate for half bin width
        if (b_map.size() > 1)
            (--b_map.end())->second /= 2; // -//-
    }

    return c;
}

} // namespace Faunus