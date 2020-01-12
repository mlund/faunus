#include "io.h"
#include "units.h"
#include "random.h"
#include "spdlog/spdlog.h"
#include <fstream>
#include <iostream>

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
    faunus_logger->warn("cannot read file: {0}", file);
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
        faunus_logger->warn("charge mismatch on atom {0} {1}: {2} atomlist value", num, name,
                            (keepcharges) ? "ignoring" : "using");
        if (not keepcharges)
            a.charge = it->charge; // let's use atomdata charge
    }

    // does radius match AtomData?
    if (std::fabs(it->sigma - 2 * radius) > pc::epsilon_dbl)
        faunus_logger->warn("radius mismatch on atom {0} {1}: using atomlist value", num, name);
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
    o << pv.size() << std::endl;
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
        faunus_logger->warn("xtc file could not be opened");
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

void FormatXTC::setbox(double box_length) { setbox(box_length, box_length, box_length); }

void FormatXTC::setbox(const Point &box_length) { setbox(box_length.x(), box_length.y(), box_length.z()); }

std::string FormatPQR::writeCryst1(const Point &box_length, const Point &angle) {
    return fmt::format("CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n", box_length.x(),
                       box_length.y(), box_length.z(), angle.x(), angle.y(), angle.z());
}

bool FormatPQR::readCrystalRecord(const std::string &record, Point &box_length) {
    std::stringstream o(record);
    std::string key;
    o >> key;
    if (key == "CRYST1") {
        o >> box_length.x() >> box_length.y() >> box_length.z();
        return true;
    } else
        return false;
}

bool FormatPQR::readAtomRecord(const std::string &record, Particle &particle, double &radius) {
    std::stringstream o(record);
    std::string key;
    o >> key;
    if (key == "ATOM" or key == "HETATM") {
        int atom_index, res_index;
        std::string atom_name, res_name;
        o >> atom_index >> atom_name;
        auto it = findName(Faunus::atoms, atom_name);
        if (it == atoms.end())
            throw std::runtime_error("PQR load error: unknown atom name '" + atom_name + "'.");
        particle = *it;
        o >> res_name >> res_index >> particle.pos.x() >> particle.pos.y() >> particle.pos.z() >> particle.charge >>
            radius;
        return true;
    }
    return false;
}

Point FormatPQR::load(const std::string &filename, FormatPQR::Tpvec &destination, bool keep_charges) {
    Point box_length(0, 0, 0);
    std::ifstream file(filename);
    if (file) {
        bool log_charge = true, log_radius = true; // emit warnings only once
        std::string record;
        int atom_index = 0;
        while (std::getline(file, record)) {
            double radius;
            Particle particle;
            if (readAtomRecord(record, particle, radius)) { // read ATOM or HETATOM record
                atom_index++;
                AtomData &atom_data = Faunus::atoms.at(particle.id);

                // does charge match AtomData?
                if (std::fabs(atom_data.charge - particle.charge) > pc::epsilon_dbl) {
                    if (log_charge) {
                        faunus_logger->warn("charge mismatches in PQR file {0}: using {1} values", filename,
                                            keep_charges ? "PQR" : "atomlist");
                        log_charge = false;
                    }
                    if (keep_charges)
                        faunus_logger->debug("charge mismatch on atom {0} {1}: using {2} (PQR) over {3} (atomlist)",
                                             atom_data.name, atom_index, particle.charge, atom_data.charge);
                    else {
                        faunus_logger->debug("charge mismatch on atom {0} {1}: using {2} (atomlist) over {3} (PQR)",
                                             atom_data.name, atom_index, atom_data.charge, particle.charge);
                        particle.charge = atom_data.charge;
                    }
                }

                // does radius match AtomData?
                if (std::fabs(atom_data.sigma - 2 * radius) > pc::epsilon_dbl) {
                    if (log_radius) {
                        faunus_logger->warn("radius mismatches in PQR file {0}: using atomlist values", filename);
                        log_radius = false;
                    }
                    faunus_logger->info("radius mismatch on atom {0} {1}: using {2} (atomlist) over {3} (PQR)",
                                        atom_data.name, atom_index, atom_data.sigma / 2, radius);
                }
                destination.push_back(particle);
            } else // is it a CRYST1 record? (unit cell information)
                readCrystalRecord(record, box_length);
        }
    }
    return box_length;
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

/**
 * @param particles Particle vector
 * @param boxlength Unit cell dimensions (optional)
 * @param n Number of atoms in each residue (default: 1e9)
 */
std::string FormatPQR::toString(const Tpvec &particles, Point boxlength, int n) {
    int nres = 1, natom = 1;
    std::ostringstream o;
    if (boxlength.norm() > 1e-6)
        o << writeCryst1(boxlength);
    for (auto &i : particles) {
        AtomData &d = atoms.at(i.id);        // atom properties as defined in topology
        Point pos = i.pos + 0.5 * boxlength; // move origin to corner of box (faunus used the middle)
        o << fmt::format("ATOM  {:5d} {:4} {:4}{:5d}    {:8.3f} {:8.3f} {:8.3f} {:.3f} {:.3f}", natom++, d.name, d.name,
                         nres, pos.x(), pos.y(), pos.z(), i.charge, 0.5 * d.sigma);
        if (d.name == "CTR")
            nres++;
        else if (natom % n == 0)
            nres++;
    }
    o << "END\n";
    return o.str();
}

/**
 * @param filename Output PQR filename
 * @param particles Particle vector
 * @param boxlength Unit cell dimensions (default: [0,0,0], not printed)
 * @param n Number of atoms in each residue (default: 1e9)
 */
bool FormatPQR::save(const std::string &filename, const Tpvec &particles, Point boxlength, int n) {
    return IO::writeFile(filename, toString(particles, boxlength, n));
}

bool FormatXYZ::save(const std::string &file, const FormatXYZ::Tpvec &particles, const Point &box) {
    std::ostringstream o;
    o << particles.size() << "\n" << box.transpose() << "\n";
    for (auto &i : particles)
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
            faunus_logger->error("mxyz load error: number of particles in xyz file {0} does not match input file ({1})",
                                 n, p.size());
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

std::vector<Particle> fastaToParticles(const std::string &fasta_sequence, double bond_length, const Point &origin) {
    std::vector<Particle> particles; // particle vector
    particles.reserve(fasta_sequence.size());
    Particle p;                                      // single particle
    p.pos = origin;                                  // first atom placed at `origin`
    for (auto id : fastaToAtomIds(fasta_sequence)) { // loop over atom ids
        p.id = id;
        p.charge = Faunus::atoms[id].charge;
        if (not particles.empty())
            p.pos = particles.back().pos + ranunit(random) * bond_length;
        particles.push_back(p);
    }
    return particles;
}

bool loadStructure(const std::string &file, std::vector<Particle> &dst, bool append, bool keep_charges) {
    if (append == false)
        dst.clear();
    std::string suffix = file.substr(file.find_last_of(".") + 1);
    if (suffix == "aam")
        FormatAAM::load(file, dst, keep_charges);
    if (suffix == "pqr")
        FormatPQR::load(file, dst, keep_charges);
    if (suffix == "xyz")
        FormatXYZ::load(file, dst);
    return !dst.empty();
}

} // namespace Faunus
