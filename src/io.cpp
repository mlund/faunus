#include <fstream>
#include "io.h"
#include "units.h"
#include "random.h"

namespace Faunus {

bool IO::readFile(const std::string &file, std::vector<std::string> &v) {
    std::ifstream f(file);
    if (f) {
        std::string s;
        while (getline(f, s))
            v.push_back(s);
        f.close();
        return true;
    }
    std::cerr << "# WARNING! FILE " << file << " NOT READ!\n";
    return false;
}

bool IO::writeFile(const std::string &file, const std::string &s, std::ios_base::openmode mode) {
    std::ofstream f(file, mode);
    if (f) {
        f << s;
        return true;
    }
    return false;
}

void IO::strip(std::vector<std::string> &v, const std::string &pat) {
    auto iter = v.begin();
    while (iter < v.end())
        if ((*iter).find(pat) != std::string::npos)
            v.erase(iter);
        else
            ++iter;
}

int FormatXTC::getNumAtoms() { return natoms_xtc; }

bool FormatAAM::keepcharges = true;

Particle &FormatAAM::s2p(const std::string &s, Particle &a) {
    std::stringstream o;
    std::string name;
    double radius, mw;
    int num;
    o << s;
    o >> name;
    auto it = findName(atoms, name);
    if (it == atoms.end())
        throw std::runtime_error("AAM load error: unknown atom name '" + name + "'.");

    a = *it;
    o >> num >> a.pos.x() >> a.pos.y() >> a.pos.z() >> a.charge >> mw >> radius;

    // does charge match AtomData?
    if (std::fabs(it->charge - a.charge) > pc::epsilon_dbl) {
        std::cerr << "charge mismatch on atom " << num << name << ": " << ((keepcharges) ? "ignoring" : "using")
                  << " `atomlist` value" << endl;
        if (not keepcharges)
            a.charge = it->charge; // let's use atomdata charge
    }

    // does radius match AtomData?
    if (std::fabs(it->sigma - 2 * radius) > pc::epsilon_dbl)
        std::cerr << "radius mismatch on atom " << num << name << ": using `atomlist` value" << endl;

    return a;
}
bool FormatAAM::load(const std::string &file, FormatAAM::Tpvec &target, bool _keepcharges) {
    keepcharges = _keepcharges;
    std::vector<std::string> v;
    target.clear();
    if (IO::readFile(file, v)) {
        IO::strip(v, "#");
        unsigned int n = atoi(v[0].c_str());
        target.resize(n);
        for (unsigned int i = 1; i <= n; i++) {
            assert(i < v.size());
            s2p(v.at(i), target.at(i - 1));
        }
        return true;
    }
    return false;
}
bool FormatAAM::save(const std::string &file, const FormatAAM::Tpvec &pv) {
    std::ostringstream o;
    o << pv.size() << endl;
    for (size_t i = 0; i < pv.size(); i++)
        o << p2s(pv[i], i);
    return IO::writeFile(file, o.str());
}

bool FormatXTC::open(std::string s) {
    if (xd != NULL)
        close();
    xd = xdrfile_open(&s[0], "r");
    if (xd != NULL) {
        int rc = read_xtc_natoms(&s[0], &natoms_xtc); // get number of atoms
        if (rc == exdrOK) {
            x_xtc = new rvec[natoms_xtc]; // resize coordinate array
            return true;
        }
    } else
        std::cerr << "# ioxtc error: xtc file could not be opened." << endl;
    return false;
}

void FormatXTC::close() {
    xdrfile_close(xd);
    xd = NULL;
    delete[] x_xtc;
}

FormatXTC::FormatXTC(double len) {
    prec_xtc = 1000.;
    time_xtc = step_xtc = 0;
    setbox(len);
    xd = NULL;
    x_xtc = NULL;
}

FormatXTC::~FormatXTC() { close(); }

void FormatXTC::setbox(double x, double y, double z) {
    assert(x > 0 && y > 0 && z > 0);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            xdbox[i][j] = 0;
    xdbox[0][0] = 0.1 * x; // corners of the
    xdbox[1][1] = 0.1 * y; // rectangular box
    xdbox[2][2] = 0.1 * z; // in nanometers!
}

void FormatXTC::setbox(double len) { setbox(len, len, len); }

void FormatXTC::setbox(const Point &p) { setbox(p.x(), p.y(), p.z()); }

std::string FormatPQR::writeCryst1(const Point &len, const Point &angle) {
    char buf[500];
    sprintf(buf, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n", len.x(), len.y(), len.z(), angle.x(),
            angle.y(), angle.z());
    return std::string(buf);
}

Point FormatPQR::load(const std::string &file, FormatPQR::Tpvec &p, bool keepcharges) {
    Point len(0, 0, 0);
    std::ifstream in(file);
    if (in) {
        Particle a;
        int iatom, ires;
        std::string line, key, aname, rname;
        while (std::getline(in, line)) {
            std::stringstream o(line);
            while (o >> key)
                if (key == "ATOM" or key == "HETATM") {
                    double radius;
                    o >> iatom >> aname;
                    auto it = findName(atoms, aname);
                    if (it == atoms.end())
                        throw std::runtime_error("PQR load error: unknown atom name '" + aname + "'.");
                    a = *it;
                    o >> rname >> ires >> a.pos.x() >> a.pos.y() >> a.pos.z() >> a.charge >> radius;

                    // does charge match AtomData?
                    if (std::fabs(it->charge - a.charge) > pc::epsilon_dbl) {
                        std::cerr << "charge mismatch on atom " << aname << " " << ires;
                        if (keepcharges)
                            std::cerr << "; using " << a.charge << " instead of `atomlist`s " << it->charge << "."
                                      << endl;
                        else {
                            std::cerr << "; using `atomlist`s " << it->charge << " over " << a.charge << "." << endl;
                            a.charge = it->charge;
                        }
                    }

                    // does radius match AtomData?
                    if (std::fabs(it->sigma - 2 * radius) > pc::epsilon_dbl)
                        std::cerr << "radius mismatch on atom " << aname << " " << ires << "; using `atomdata`s "
                                  << it->sigma / 2 << " over " << radius << "." << endl;

                    p.push_back(a);

                } else if (key == "CRYST1")
                    o >> len.x() >> len.y() >> len.z();
        }
    }
    return len;
}

void FormatPQR::load(const std::string &file, std::vector<FormatPQR::Tpvec> &v) {
    std::ifstream in(file);
    if (in) {
        Particle a;
        Tpvec p;
        int iatom, ires;
        std::string line, key, aname, rname;
        while (std::getline(in, line)) {
            std::stringstream o(line);
            while (o >> key)
                if (key == "ATOM" or key == "HETATM") {
                    double radius = 0;
                    o >> iatom >> aname;
                    auto it = findName(atoms, aname);
                    if (it == atoms.end())
                        throw std::runtime_error("PQR load error: unknown atom name '" + aname + "'.");
                    a = *it;
                    o >> rname >> ires >> a.pos.x() >> a.pos.y() >> a.pos.z() >> a.charge >> radius;
                    p.push_back(a);
                } else if (key == "END") {
                    v.push_back(p);
                    p.clear();
                    p.reserve(v.back().size());
                }
        }
    } else
        throw std::runtime_error("Error loading PQR trajectory " + file);
}

bool FormatPQR::save(const std::string &file, const FormatPQR::Tpvec &p, Point len, int n) {
    int nres = 1, natom = 1;
    char buf[500];
    std::ostringstream o;
    if (len.norm() > 1e-6)
        o << writeCryst1(len);
    for (auto &p_i : p) {
        auto &prop = atoms.at(p_i.id);
        std::string name = prop.name;
        double radius = prop.sigma / 2;
        sprintf(buf, "ATOM  %5d %-4s %-4s%5d    %8.3f %8.3f %8.3f %.3f %.3f\n", natom++, name.c_str(), name.c_str(),
                nres, (p_i.pos + len / 2).x(), (p_i.pos + len / 2).y(), (p_i.pos + len / 2).z(), p_i.charge,
                radius); // move particles inside the sim. box
        o << buf;
        if (atoms.at(p_i.id).name == "CTR")
            nres++;
        else if (natom % n == 0)
            nres++;
    }
    o << "END\n";
    return IO::writeFile(file, o.str());
}

bool FormatXYZ::save(const std::string &file, const FormatXYZ::Tpvec &p, const Point &box) {
    std::ostringstream o;
    o << p.size() << "\n" << box.transpose() << "\n";
    for (auto &i : p)
        o << atoms.at(i.id).name << " " << i.pos.transpose() << "\n";
    return IO::writeFile(file, o.str());
}

bool FormatXYZ::load(const std::string &file, FormatXYZ::Tpvec &p, bool append) {
    std::ifstream f(file);
    if (f) {
        if (append == false)
            p.clear();
        Particle a;
        std::string comment, name;
        size_t n;
        f >> n;
        p.reserve(p.size() + n);
        std::getline(f, comment); // ">>" token doesn't gobble new line
        std::getline(f, comment); // read comment line
        for (size_t i = 0; i < n; i++) {
            f >> name;
            auto it = findName(atoms, name);
            if (it == atoms.end())
                throw std::runtime_error("XYZ load error: unknown atom name '" + name + "'.");
            a = *it;
            f >> a.pos.x() >> a.pos.y() >> a.pos.z();
            p.push_back(a);
        }
        if (!p.empty())
            return true;
    }
    return false;
}

std::string FormatMXYZ::p2s(const Particle &, int) {
    std::ostringstream o;
    o.precision(5);
    // o << a.transpose() << " " << a.dir.transpose() << " " << a.patchdir.transpose() << "\n";
    return o.str();
}

Particle &FormatMXYZ::s2p(const std::string &s, Particle &a) {
    std::stringstream o;
    o << s;
    // o >> a.x() >> a.y() >> a.z() >> a.dir.x() >> a.dir.y() >> a.dir.z() >> a.patchdir.x() >> a.patchdir.y() >>
    //    a.patchdir.z();
    // a.init();
    return a;
}

bool FormatMXYZ::save(const std::string &file, const FormatMXYZ::Tpvec &p, const Point &len, int time) {
    std::ostringstream o;
    o << p.size() << "\n"
      << "sweep " << time << "; box " << len.transpose() << "\n";
    for (size_t i = 0; i < p.size(); i++)
        o << p2s(p[i], i);
    return IO::writeFile(file, o.str(), std::ios_base::app);
}

bool FormatMXYZ::load(const std::string &file, FormatMXYZ::Tpvec &p, Point &len) {
    std::stringstream o;
    std::vector<std::string> v;
    if (IO::readFile(file, v)) {
        IO::strip(v, "#");
        size_t n = atoi(v[0].c_str());
        if (p.size() != n)
            std::cerr << "# mxyz load error: number of particles in xyz file " << n << " does not match input file ("
                      << p.size() << ")!" << endl;
        o << v[1].erase(0, v[1].find_last_of("x") + 1);
        o >> len.x() >> len.y() >> len.z();
        for (size_t i = 2; i < n + 2; i++)
            s2p(v.at(i), p.at(i - 2));
        return true;
    }
    return false;
}

void FormatGRO::s2p(const std::string &s, Particle &dst) {
    std::stringstream o;
    std::string name;
    o << s.substr(10, 5) << s.substr(20, 8) << s.substr(28, 8) << s.substr(36, 8);
    o >> name;
    auto it = findName(atoms, name);
    if (it != atoms.end()) {
        o >> dst.pos.x() >> dst.pos.y() >> dst.pos.z();
        dst = *it;
    } else
        throw std::runtime_error("gro: unknown atom name");
    dst.pos = 10 * dst.pos; // nm->angstrom
}

bool FormatGRO::load(const std::string &file, FormatGRO::Tpvec &p) {
    p.clear();
    v.resize(0);
    if (IO::readFile(file, v)) {
        int last = atoi(v[1].c_str()) + 1;
        for (int i = 2; i <= last; i++) {
            Particle a;
            s2p(v[i], a);
            p.push_back(a);
        }
        return true;
    }
    return false;
}

bool FormatGRO::save(const std::string &file, FormatGRO::Tpvec &p, std::string mode) {
    int nres = 1, natom = 1;
    char buf[79];
    double halflen = len / 20; // halflen in nm
    std::ostringstream o;
    o << "Generated by Faunus -- http://faunus.sourceforge.net" << std::endl << p.size() << std::endl;
    for (auto &pi : p) {
        std::string name = atoms.at(pi.id).name;
        sprintf(buf, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n", nres, name.c_str(), name.c_str(), natom++,
                pi.pos.x() / 10 + halflen, pi.pos.y() / 10 + halflen, pi.pos.z() / 10 + halflen);
        o << buf;
        if (atoms[pi.id].name == "CTR")
            nres++;
    }
    if (len > 0)
        o << len / 10 << " " << len / 10 << " " << len / 10 << std::endl; // box side length in nm
    if (mode == "append")
        return IO::writeFile(file, o.str(), std::ios_base::app);
    else
        return IO::writeFile(file, o.str(), std::ios_base::out);
}

std::vector<Particle> fastaToParticles(const std::string &fasta, double spacing, const Point &origin) {
    std::vector<Particle> vec; // particle vector
    Particle p;                // single particle
    p.pos = origin;            // first atom place here
    auto ids = fastaToAtomIds(fasta);
    for (auto i : ids) { // loop over atom ids
        p.id = i;
        p.charge = atoms[i].charge;
        if (not vec.empty())
            p.pos = vec.back().pos + ranunit(random) * spacing;
        vec.push_back(p);
    }
    return vec;
}

bool loadStructure(const std::string &file, std::vector<Particle> &dst, bool append, bool keepcharges) {
    if (append == false)
        dst.clear();
    std::string suffix = file.substr(file.find_last_of(".") + 1);
    if (suffix == "aam")
        FormatAAM::load(file, dst, keepcharges);
    if (suffix == "pqr")
        FormatPQR::load(file, dst, keepcharges);
    if (suffix == "xyz")
        FormatXYZ::load(file, dst);
    if (not dst.empty())
        return true;
    return false;
}

} // namespace Faunus
