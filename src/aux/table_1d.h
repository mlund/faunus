#pragma once

#include <Eigen/Core>
#include <vector>
#include <fstream>

namespace Faunus {

/**
 * @brief Dynamic table for 1d data
 * @todo This really needs documentation and general refactoring!
 */
template <typename Tcoeff = double, typename base = Eigen::Matrix<Tcoeff, Eigen::Dynamic, Eigen::Dynamic>>
class Table : public base {
  private:
    using Tvec = std::vector<double>;
    Tvec _bw, _lo, _hi;
    using index_t = typename base::Index;
    static_assert(std::is_integral_v<index_t> && std::is_signed_v<index_t>, "signed integral value expected");
    index_t _rows;
    index_t _cols;

  public:
    explicit Table(const Tvec& bw = {1, 1}, const Tvec& lo = {0, 0}, const Tvec& hi = {2, 2})
    {
        reInitializer(bw, lo, hi);
    }

    // required for assignment from Eigen::Matrix and Eigen::Array objects
    template <typename OtherDerived>
    explicit Table(const Eigen::MatrixBase<OtherDerived>& other)
        : base(other)
    {
    }
    template <typename OtherDerived> Table& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->base::operator=(other);
        return *this;
    }
    template <typename OtherDerived>
    explicit Table(const Eigen::ArrayBase<OtherDerived>& other)
        : base(other)
    {
    }
    template <typename OtherDerived> Table& operator=(const Eigen::ArrayBase<OtherDerived>& other)
    {
        this->base::operator=(other);
        return *this;
    }

    void reInitializer(const Tvec& bw, const Tvec& lo, const Tvec& hi)
    {
        assert(bw.size() == 1 || bw.size() == 2);
        assert(bw.size() == lo.size() && lo.size() == hi.size());
        _bw = bw;
        _lo = lo;
        _hi = hi;
        _rows = (_hi[0] - _lo[0]) / _bw[0] + 1.;
        if (bw.size() == 2) {
            _cols = (_hi[1] - _lo[1]) / _bw[1] + 1;
        }
        else {
            _cols = 1;
        }
        base::resize(_rows, _cols);
        base::setZero();
    }

    void round(Tvec& v) const
    {
        for (Tvec::size_type i = 0; i != v.size(); ++i) {
            v[i] = (v[i] >= 0) ? int(v[i] / _bw[i] + 0.5) * _bw[i] : int(v[i] / _bw[i] - 0.5) * _bw[i];
        }
    }

    void to_index(Tvec& v) const
    {
        for (Tvec::size_type i = 0; i != v.size(); ++i) {
            v[i] = (v[i] >= 0) ? int(v[i] / _bw[i] + 0.5) : int(v[i] / _bw[i] - 0.5);
            v[i] = v[i] - _lo[i] / _bw[i];
        }
        v.resize(2, 0);
    }

    Tcoeff& operator[](const Tvec& v)
    {
        return base::operator()(static_cast<index_t>(v[0]), static_cast<index_t>(v[1]));
    }

    bool isInRange(const Tvec& v) const
    {
        bool in_range = true;
        for (Tvec::size_type i = 0; i != v.size(); ++i) {
            in_range = in_range && v[i] >= _lo[i] && v[i] <= _hi[i];
        }
        return in_range;
    }

    Tvec hist2buf(int) const
    {
        Tvec sendBuf;
        for (index_t i = 0; i < _cols; ++i) {
            for (index_t j = 0; j < _rows; ++j) {
                sendBuf.push_back(base::operator()(j, i));
            }
        }
        return sendBuf;
    }

    void buf2hist(const Tvec& v)
    {
        assert(!v.empty());
        base::setZero();
        auto p = static_cast<int>(v.size()) / this->size();
        Tvec::size_type n = 0;
        double nproc = p;
        while (p-- > 0) {
            for (index_t i = 0; i < _cols; ++i) {
                for (index_t j = 0; j < _rows; ++j) {
                    base::operator()(j, i) += v.at(n++) / nproc;
                }
            }
        }
    }

    base getBlock(const Tvec& slice)
    { // {xmin,xmax} or {xmin,xmax,ymin,ymax}
        Tvec w(4, 0);
        switch (slice.size()) {
        case 1:
            w[0] = w[1] = (slice[0] - _lo[0]) / _bw[0];
            w[3] = _cols - 1;
            break;
        case 2:
            w[0] = (slice[0] - _lo[0]) / _bw[0];
            w[1] = (slice[1] - _lo[0]) / _bw[0];
            break;
        case 3:
            w[0] = w[1] = _rows - 1;
            w[2] = w[3] = _cols - 1;
            break;
        case 4:
            w[0] = (slice[0] - _lo[0]) / _bw[0];
            w[1] = (slice[1] - _lo[0]) / _bw[0];
            w[2] = (slice[2] - _lo[1]) / _bw[1];
            w[3] = (slice[3] - _lo[1]) / _bw[1];
        }
        return this->block(w[0], w[2], w[1] - w[0] + 1, w[3] - w[2] + 1); // xmin,ymin,rows,cols
    }

    Tcoeff avg(const Tvec& v) const { return this->getBlock(v).mean(); }

    void save(const std::string& filename, Tcoeff scale = 1, Tcoeff translate = 0) const
    {
        Eigen::VectorXd v1(_cols + 1);
        Eigen::VectorXd v2(_rows + 1);
        v1(0) = v2(0) = base::size();
        for (index_t i = 1; i != _cols + 1; ++i) {
            v1(i) = (i - 1) * _bw[1] + _lo[1];
        }
        for (index_t i = 1; i != _rows + 1; ++i) {
            v2(i) = (i - 1) * _bw[0] + _lo[0];
        }
        base m(_rows + 1, _cols + 1);
        m.leftCols(1) = v2;
        m.topRows(1) = v1.transpose();
        m.bottomRightCorner(_rows, _cols) = *this;
        if (scale != 1) {
            m.bottomRightCorner(_rows, _cols) *= scale;
        }
        if (translate != 0) {
            m.bottomRightCorner(_rows, _cols) += base::Constant(_rows, _cols, translate);
        }
        std::ofstream stream(filename.c_str());
        if (stream) {
            stream.precision(10);
            if (_cols == 1) {
                stream << "#";
            }
            stream << m;
        }
    }

    void saveRow(const std::string& filename, const Tvec& v, Tcoeff scale = 1, Tcoeff translate = 0)
    {
        if (!this->isInRange(v)) {
            return;
        }
        auto b = this->getBlock(v);
        auto size = b.size();
        Eigen::VectorXd w(size);
        for (index_t i = 0; i != size; ++i) {
            w(i) = i * _bw[1] + _lo[1];
        }
        base m(size, 2);
        m.leftCols(1) = w;
        m.bottomRightCorner(size, 1) = b.transpose();
        if (scale != 1) {
            m.bottomRightCorner(size, 1) *= scale;
        }
        if (translate != 0) {
            m.bottomRightCorner(size, 1) += base::Constant(size, 1, translate);
        }
        std::ofstream stream(filename.c_str());
        if (stream) {
            stream.precision(10);
            stream << m;
        }
    }

    void load(const std::string& filename)
    {
        if (std::ifstream f(filename.c_str()); f) {
            index_t i = 0;
            index_t j = -1;
            std::string line;
            getline(f, line);
            while (getline(f, line)) {
                if (i > _rows - 1 || j > _cols - 1) {
                    throw std::runtime_error("file larger than expected:"s + filename);
                }
                j = -1;
                const std::istringstream iss(line);
                Tcoeff a;
                Tcoeff b;
                iss >> a;
                while (iss >> b) {
                    base::operator()(i, ++j) = b;
                }
                ++i;
            }
            if (i != _rows || j != _cols - 1) {
                throw std::runtime_error("file larger than expected:"s + filename);
            }
        }
    }
};

} // namespace Faunus