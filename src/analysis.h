#pragma once

#include "space.h"
#include "io.h"
#include "scatter.h"
#include "reactioncoordinate.h"
#include "auxiliary.h"
#include <set>

namespace Faunus {

namespace Energy {
class Hamiltonian;
class Energybase;
}

namespace Analysis {

/**
 * @brief Base class for all analysis functions
 *
 * Analysis classes are used to perform system analysis
 * over the course or a simulations or for post-analysis
 * of an existing trajectory. The base class adds basic
 * functionality such as timing, number of steps, json IO
 * and enforce a common interface.
 */
class Analysisbase {
  private:
    virtual void _to_json(json &) const;
    virtual void _from_json(const json &);
    virtual void _sample() = 0;
    virtual void _to_disk(); //!< save data to disk
    int stepcnt = 0;
    int totstepcnt = 0;
    TimeRelativeOfTotal<std::chrono::microseconds> timer;

  protected:
    int steps = 0; //!< Sample interval (do not modify)
    int nskip = 0; //!< MC steps to skip before sampling
    int cnt = 0;   //!< number of samples

  public:
    std::string name; //!< descriptive name
    std::string cite; //!< url, doi etc. describing the analysis

    void to_json(json &) const;    //!< JSON report w. statistics, output etc.
    void from_json(const json &);  //!< configure from json object
    void to_disk();                //!< Save data to disk (if defined)
    virtual void sample();
    virtual ~Analysisbase() = default;
};

void to_json(json &j, const Analysisbase &base);

/*
 * @brief Sample and save reaction coordinates to a file
 */
class FileReactionCoordinate : public Analysisbase {
  private:
    Average<double> avg;
    std::string type, filename;
    std::ofstream file;
    std::shared_ptr<ReactionCoordinate::ReactionCoordinateBase> rc = nullptr;

    void _to_json(json &j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    FileReactionCoordinate(const json &j, Space &spc);
};

/**
 * @brief Excess chemical potential of molecules
 */
class WidomInsertion : public Analysisbase {
    Space &spc;
    Energy::Hamiltonian *pot;
    RandomInserter rins;
    std::string molname; // molecule name
    int ninsert;
    int molid; // molecule id
    bool absolute_z = false;
    Average<double> expu;
    Change change;

    void _sample() override;
    void _to_json(json &j) const override;
    void _from_json(const json &j) override;

  public:
    WidomInsertion(const json &j, Space &spc, Energy::Hamiltonian &pot);
};

/**
 * @brief Density of atom along axis
 *
 * Calculates the summed density of `atoms` in spherical, cylindrical or planar shells around
 * `origo` which by default is the center of the simulation box
 */
class AtomProfile : public Analysisbase {
    Space &spc;
    Equidistant2DTable<double, double> tbl;
    std::vector<std::string> names; // atom names to analyse
    std::set<int> ids;              // atom ids to analyse
    std::string file;               // output filename
    Point ref = {0, 0, 0};
    Eigen::Vector3i dir = {1, 1, 1};
    double dr; // radial resolution
    bool count_charge = false;
    std::string atom_com;
    int id_com = -1; // center at COM of id_com atoms?

    void _from_json(const json &j) override;
    void _to_json(json &j) const override;
    void _to_disk() override;
    void _sample() override;

  public:
    AtomProfile(const json &j, Space &spc);
};

/**
 * @brief Measures the density of atoms along z axis
 */
class SlicedDensity : public Analysisbase {
    Space &spc;
    Table2D<double, unsigned int> N; // N(z)
    std::vector<std::string> names;
    std::vector<int> ids;
    std::string file;
    double dz;
    std::string atom_com;
    int id_com = -1; // center at COM of id_com atoms?

    void _from_json(const json &j) override;
    void _to_json(json &j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    SlicedDensity(const json &j, Space &spc);
};

/**
 * @brief Analysis of particle densities
 */
class Density : public Analysisbase {
    Space &spc;
    typedef Equidistant2DTable<unsigned int, double> Ttable; // why double?

    std::map<int, Ttable> swpdhist; // Probability density of swapping atoms
    std::map<int, Ttable> atmdhist; // Probability density of atomic molecules
    std::map<int, Ttable> moldhist; // Probability density of polyatomic molecules
    std::map<int, Average<double>> rho_mol, rho_atom;
    std::map<int, int> Nmol, Natom;
    Average<double> Lavg, Vavg, invVavg;

    // int capacity_limit = 10; // issue warning if capacity get lower than this

    void _sample() override;
    void _to_json(json &) const override;
    void _to_disk() override;

  public:
    Density(const json &, Space &);
};

class ChargeFluctuations : public Analysisbase {
  private:
    Space &spc;
    typename decltype(Faunus::molecules)::const_iterator mol_iter; // selected molecule type

    std::vector<std::map<int, int>> idcnt; // populations of types of atomic indexes
    std::vector<Average<double>> charge;   // average charges of atomic indexes
    std::string file;                      // name of PQR file with average charges
    bool verbose;                          // set to true for more output

    void _sample() override;

    void _to_json(json &) const override;

    /**
     * @brief Saves average PQR file to disk if `pqrfile` input is given
     */
    void _to_disk() override;

  public:
    ChargeFluctuations(const json &j, Space &spc);
}; // Fluctuations of atomic charges

class Multipole : public Analysisbase {
    const Space &spc;
    struct data {
        Average<double> Z, Z2, mu, mu2;
    };
    std::map<int, data> _map; //!< Molecular moments and their fluctuations

    void _sample() override;
    void _to_json(json &) const override;

  public:
    Multipole(const json &, const Space &);
}; // Molecular multipoles and their fluctuations

class SystemEnergy : public Analysisbase {
    std::string file, sep = " ";
    std::ofstream f;
    std::function<std::vector<double>()> energyFunc;
    Average<double> uavg, u2avg; //!< mean energy and mean squared energy
    std::vector<std::string> names;
    Table2D<double, double> ehist; // Density histograms
    double uinit;

    void normalize();
    void _sample() override;
    void _to_json(json &) const override;
    void _from_json(const json &) override;
    void _to_disk() override;

  public:
    SystemEnergy(const json &, Energy::Hamiltonian &);
}; //!< Save system energy to disk. Keywords: `nstep`, `file`.

/**
 * @brief Checks if system is sane. If not, abort program.
 */
class SanityCheck : public Analysisbase {
  private:
    Space &spc;
    void _sample() override;

  public:
    SanityCheck(const json &, Space &);
};

class SaveState : public Analysisbase {
  private:
    std::function<void(std::string)> writeFunc = nullptr;
    std::string file;
    bool saverandom;
    void _to_json(json &) const override;
    void _to_disk() override;
    void _sample() override;

  public:
    SaveState(const json &, Space &);
};

/**
 * @brief Base class for distribution functions etc.
 */
class PairFunctionBase : public Analysisbase {
  protected:
    int dim = 3;            // dimentions to use when normalizing
    int id1 = -1, id2 = -1; // particle id (mol or atom)
    double dr = 0;          // distance resolution
    Eigen::Vector3i slicedir = {0, 0, 0};
    double thickness = 0;
    Equidistant2DTable<double, double> hist;
    std::string name1, name2; // atom/molecule names
    std::string file;         // output filename
    double Rhypersphere = -1; // Radius of 2D hypersphere
    Average<double> V;        // average volume (angstrom^3)

  private:
    void _from_json(const json &) override;
    void _to_json(json &) const override;
    void _to_disk() override;

  public:
    PairFunctionBase(const json &);
};

class PairAngleFunctionBase : public PairFunctionBase {
  protected:
    Equidistant2DTable<double, Average<double>> hist2;

  private:
    void _from_json(const json &) override;
    void _to_disk() override;

  public:
    PairAngleFunctionBase(const json &);
};

/** @brief Atomic radial distribution function, g(r) */
class AtomRDF : public PairFunctionBase {
    Space &spc;

    void _sample() override;

  public:
    AtomRDF(const json &, Space &);
};

/** @brief Same as `AtomRDF` but for molecules. Identical input. */
class MoleculeRDF : public PairFunctionBase {
    Space &spc;
    void _sample() override;
  public:
    MoleculeRDF(const json &, Space &);
};

/** @brief Dipole-dipole correlation function, <\boldsymbol{\mu}(0)\cdot\boldsymbol{\mu}(r)> */
class AtomDipDipCorr : public PairAngleFunctionBase {
    Space &spc;
    void _sample() override;
  public:
    AtomDipDipCorr(const json &, Space &);
};

/** @brief Write XTC trajectory file */
class XTCtraj : public Analysisbase {
    std::vector<int> molids;        // molecule ids to save to disk
    std::vector<std::string> names; // molecule names of above
    std::function<bool(Particle &)> filter; // function to filter molecule ids

    void _to_json(json &) const override;
    void _from_json(const json &) override;

    FormatXTC xtc;
    Space &spc;
    std::string file;

    void _sample() override;

  public:
    XTCtraj(const json &j, Space &s);
};

/**
 * @brief Excess pressure using virtual volume move
 */
class VirtualVolume : public Analysisbase {
    std::string file; // output filename
    std::ofstream output_file; // output filestream
    double dV; // volume perturbation
    Change c;
    Energy::Energybase &pot;
    std::function<double()> getVolume;
    std::function<void(double)> scaleVolume;
    Average<double> duexp; // < exp(-du/kT) >

    void _sample() override;
    void _from_json(const json &j) override;
    void _to_json(json &j) const override;
    void _to_disk() override;

  public:
    VirtualVolume(const json &j, Space &spc, Energy::Energybase &pot);
};

/**
 * @brief Multipolar decomposition between groups as a function of separation
 * @date Malmo 2014
 * @todo Add option to use charge center instead of mass center
 */
class MultipoleDistribution : public Analysisbase {
    typedef typename Space::Tgroup Tgroup;

    struct data {
        Average<double> tot, ii, id, iq, dd, mucorr;
    };

    std::vector<std::string> names; //!< Molecule names (len=2)
    std::vector<int> ids;           //!< Molecule ids (len=2)
    std::string filename;           //!< output file name
    // int id1, id2;                   //!< pair of molecular id's to analyse
    double dr;             //!< distance resolution
    std::map<int, data> m; //!< Energy distributions
    Space &spc;

    double g2g(const Tgroup &g1, const Tgroup &g2); //<! exact ion-ion energy between particles
    void save() const;                              //!< save to disk
    void _sample() override;
    void _to_json(json &j) const override;
    void _to_disk() override;

  public:
    MultipoleDistribution(const json &j, Space &spc);

}; // end of multipole distribution

/**
 * @brief Sample scattering intensity
 */
class ScatteringFunction : public Analysisbase {
  private:
    enum Schemes { DEBYE, EXPLICIT }; // two different schemes
    Schemes scheme = DEBYE;
    Space &spc;
    bool usecom;                    // scatter from mass center, only?
    std::string filename;           // output file name
    std::vector<Point> p;           // vector of scattering points
    std::vector<int> ids;           // Molecule ids
    std::vector<std::string> names; // Molecule names
    typedef Scatter::FormFactorUnity<double> Tformfactor;

    std::shared_ptr<Scatter::DebyeFormula<Tformfactor>> debye;
    std::shared_ptr<Scatter::StructureFactor<double>> explicit_average;

    void _sample() override;
    void _to_disk() override;
    void _to_json(json &j) const override;

  public:
    ScatteringFunction(const json &j, Space &spc);
};

/*
 * @brief Sample and save gyration eigenvalues of all particles having the same id
 */
class AtomInertia : public Analysisbase {
  private:
    Space &spc;
    std::string filename;
    int index; // atom id
    std::ofstream file;

    Point compute();
    void _to_json(json &j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    AtomInertia(const json &j, Space &spc);
};

/*
 * @brief Sample and save the eigenvalues of the inertia tensor for a range of indexes within a molecule
 */
class InertiaTensor : public Analysisbase {
  private:
    Space &spc;
    std::string filename;
    std::vector<size_t> indexes; // range of indexes within the group
    int index; // group indes
    std::ofstream file;
    struct Data {
        Point eivals, eivec; // eigenvalues and principal axis
    };

    Data compute();
    void _to_json(json &j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    InertiaTensor(const json &j, Space &spc);
};

/*
 * @brief Sample and save charge, dipole and quadrupole moments for a range of indexes within a molecule
 */
class MultipoleMoments : public Analysisbase {
  private:
    Space &spc;
    std::string filename;
    std::vector<size_t> indexes; // range of indexes within the group
    size_t index; // group index
    bool mol_cm;
    std::ofstream file;
    struct Data {
        int q = 0; // total charge
        Point mu {0,0,0}; // dipole vector
        Point eivals, eivec, center; // quadrupole eigenvalues and major axis
    };

    Data compute();
    void _to_json(json &j) const override;
    void _sample() override;
    void _to_disk() override;

  public:
    MultipoleMoments(const json &j, Space &spc);
};

/**
 * @brief Analysis of polymer shape - radius of gyration, shape factor etc.
 * @date November, 2011
 *
 * This will analyse polymer groups and calculate Rg, Re and the shape factor. If
 * sample() is called with different groups these will be distinguished by their
 * *name* and sampled individually.
 */
class PolymerShape : public Analysisbase {
    Space &spc;
    std::map<int, Average<double>> Rg2, Rg, Re2, Re, Rs, Rs2, Rg2x, Rg2y, Rg2z;
    std::vector<int> ids; // molecule id's to analyse

    void _to_json(json &j) const override;
    Point vectorgyrationRadiusSquared(typename Space::Tgroup &g) const;
    void _sample() override;

  public:
    PolymerShape(const json &j, Space &spc);
};

/**
 * @brief "Trajectory" with charge and radius, only, for all (active, inactive) particles
 *
 * For use w. VMD to visualize charge fluctuations and grand canonical ensembles
 */
class QRtraj : public Analysisbase {
  private:
    std::string file;
    std::ofstream f;
    std::function<void()> write_to_file;
    void _sample() override;
    void _to_json(json &j) const override;
    void _to_disk() override;

  public:
    QRtraj(const json &j, Space &spc);
};

struct CombinedAnalysis : public BasePointerVector<Analysisbase> {
    CombinedAnalysis(const json &j, Space &spc, Energy::Hamiltonian &pot);
    void sample();
    void to_disk(); // prompt all analysis to safe to disk if appropriate
}; //!< Aggregates analysis

/** @brief Example analysis */
template <class T, class Enable = void> struct _analyse {
    void sample(T &) { std::cout << "not a dipole!" << std::endl; } //!< Sample
};                                                                  // primary template

/** @brief Example analysis */
template <class T> struct _analyse<T, typename std::enable_if<std::is_base_of<Dipole, T>::value>::type> {
    void sample(T &) { std::cout << "dipole!" << std::endl; } //!< Sample
};                                                            // specialized template

} // namespace Analysis

} // namespace Faunus
