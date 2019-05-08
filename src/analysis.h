#pragma once

#include "space.h"
#include "io.h"
#include "mpi.h"
#include "scatter.h"
#include "reactioncoordinate.h"
#include "auxiliary.h"

namespace Faunus {

namespace Energy {
class Hamiltonian;
class Energybase;
}

namespace Analysis {

class Analysisbase {
  private:
    virtual void _to_json(json &j) const;
    virtual void _from_json(const json &j);
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
    std::string cite; //!< reference, url, doi etc. describing the analysis

    void to_json(json &j) const;   //!< JSON report w. statistics, output etc.
    void from_json(const json &j); //!< configure from json object
    void to_disk();                //!< Save data to disk (if defined)
    virtual void sample();
    virtual ~Analysisbase();
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

  public:
    FileReactionCoordinate(const json &j, Tspace &spc);
};

/**
 * @brief Excess chemical potential of molecules
 */
class WidomInsertion : public Analysisbase {
    typedef typename Tspace::Tpvec Tpvec;

    Tspace &spc;
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
    WidomInsertion(const json &j, Tspace &spc, Energy::Hamiltonian &pot);
};

class AtomProfile : public Analysisbase {
    Tspace &spc;
    Equidistant2DTable<double, double> tbl;
    std::vector<std::string> names; // atom names to analyse
    std::set<int> ids;              // atom ids to analyse
    std::string file;               // output filename
    Point ref = {0, 0, 0};
    double dr; // radial resolution
    bool count_charge = false;
    // bool Vnormalise = true;

    void _from_json(const json &j) override;
    void _to_json(json &j) const override;
    void _sample() override;

  public:
    AtomProfile(const json &j, Tspace &spc);
    ~AtomProfile();
};

/**
 * @brief Measures the density of atoms along z axis
 */
class SlicedDensity : public Analysisbase {
    Tspace &spc;
    Table2D<double, unsigned int> N; // N(z)
    std::vector<std::string> names;
    std::vector<int> ids;
    std::string file;
    double dz;
    std::string atomCOM;
    int idCOM = -1; // center at COM of idCOM atoms?

    void _from_json(const json &j) override;
    void _to_json(json &j) const override;
    void _sample() override;

  public:
    SlicedDensity(const json &j, Tspace &spc);
    ~SlicedDensity();
};

/**
 * @brief Analysis of particle densities
 */
class Density : public Analysisbase {
    Tspace &spc;
    typedef typename Tspace::Tpvec Tpvec;
    typedef Equidistant2DTable<unsigned int, double> Ttable; // why double?

    std::map<int, Ttable> swpdhist; // Probability density of swapping atoms
    std::map<int, Ttable> atmdhist; // Probability density of atomic molecules
    std::map<int, Ttable> moldhist; // Probability density of polyatomic molecules
    std::map<int, Average<double>> rho_mol, rho_atom;
    std::map<int, int> Nmol, Natom;
    Average<double> Lavg, Vavg, invVavg;

    // int capacity_limit = 10; // issue warning if capacity get lower than this

    void _sample() override;
    void _to_json(json &j) const override;

  public:
    Density(const json &j, Tspace &spc);
    virtual ~Density();
};

class ChargeFluctuations : public Analysisbase {
  private:
    Tspace &spc;
    typedef typename Tspace::Tpvec Tpvec;
    typename decltype(Faunus::molecules)::const_iterator mol_iter; // selected molecule type

    std::vector<std::map<int, int>> idcnt; // populations of types of atomic indexes
    std::vector<Average<double>> charge;   // average charges of atomic indexes
    std::string file;                      // name of PQR file with average charges
    bool verbose;                          // set to true for more output

    void _sample() override;

    void _to_json(json &j) const override;

    /**
     * @brief Saves average PQR file to disk if `pqrfile` input is given
     */
    void _to_disk() override;

  public:
    ChargeFluctuations(const json &j, Tspace &spc);
}; // Fluctuations of atomic charges

class Multipole : public Analysisbase {
    const Tspace &spc;
    struct data {
        Average<double> Z, Z2, mu, mu2;
    };
    std::map<int, data> _map; //!< Molecular moments and their fluctuations

    void _sample() override;
    void _to_json(json &j) const override;

  public:
    Multipole(const json &j, const Tspace &spc);
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
    void _to_json(json &j) const override;
    void _from_json(const json &j) override;

  public:
    SystemEnergy(const json &j, Energy::Hamiltonian &pot);
}; //!< Save system energy to disk. Keywords: `nstep`, `file`.

/**
 * @brief Checks if system is sane. If not, abort program.
 */
class SanityCheck : public Analysisbase {
  private:
    Tspace &spc;
    void _sample() override;

  public:
    SanityCheck(const json &j, Tspace &spc);
};

class SaveState : public Analysisbase {
  private:
    std::function<void(std::string)> writeFunc = nullptr;
    std::string file;
    bool saverandom;
    void _to_json(json &j) const override;
    void _sample() override;

  public:
    SaveState(const json &j, Tspace &spc);
    ~SaveState();
};

/**
 * @brief Base class for distribution functions etc.
 */
class PairFunctionBase : public Analysisbase {
  protected:
    int dim = 3;
    int id1 = -1, id2 = -1; // particle id (mol or atom)
    double dr = 0;
    Eigen::Vector3i slicedir = {0, 0, 0};
    double thickness = 0;
    Equidistant2DTable<double, double> hist;
    std::string name1, name2, file;
    double Rhypersphere = -1; // Radius of 2D hypersphere
    Average<double> V;        // average volume (angstrom^3)

  private:
    void _from_json(const json &j) override;
    void _to_json(json &j) const override;

  public:
    PairFunctionBase(const json &j);
    virtual ~PairFunctionBase();
};

/** @brief Atomic radial distribution function, g(r) */
class AtomRDF : public PairFunctionBase {
    Tspace &spc;

    void _sample() override;

  public:
    AtomRDF(const json &j, Tspace &spc);
};

/** @brief Same as `AtomRDF` but for molecules. Identical input. */
class MoleculeRDF : public PairFunctionBase {
    typedef typename Tspace::Tpvec Tpvec;
    Tspace &spc;

    void _sample() override;

  public:
    MoleculeRDF(const json &j, Tspace &spc);
};

/** @brief Write XTC trajectory file */
class XTCtraj : public Analysisbase {
    typedef typename Tspace::Tparticle Tparticle;
    std::vector<int> molids;        // molecule ids to save to disk
    std::vector<std::string> names; // molecule names of above
    std::function<bool(Tparticle &)> filter = [](Tparticle &) { return true; };

    void _to_json(json &j) const override;
    void _from_json(const json &j) override;

    FormatXTC xtc;
    Tspace &spc;
    std::string file;

    void _sample() override;

  public:
    XTCtraj(const json &j, Tspace &s);
};

class VirtualVolume : public Analysisbase {
    double dV;
    Change c;
    Energy::Energybase &pot;
    std::function<double()> getVolume;
    std::function<void(double)> scaleVolume;
    Average<double> duexp; // < exp(-du/kT) >

    void _sample() override;
    void _from_json(const json &j) override;
    void _to_json(json &j) const override;

  public:
    VirtualVolume(const json &j, Tspace &spc, Energy::Energybase &pot);
}; //!< Excess pressure using virtual volume move

/**
 * @brief Multipolar decomposition between groups as a function of separation
 * @date Malmo 2014
 * @todo Add option to use charge center instead of mass center
 */
class MultipoleDistribution : public Analysisbase {
    typedef typename Tspace::Tgroup Tgroup;
    typedef typename Tspace::Tparticle Tparticle;

    struct data {
        Average<double> tot, ii, id, iq, dd, mucorr;
    };

    std::vector<std::string> names; //!< Molecule names (len=2)
    std::vector<int> ids;           //!< Molecule ids (len=2)
    std::string filename;           //!< output file name
    // int id1, id2;                   //!< pair of molecular id's to analyse
    double dr;             //!< distance resolution
    std::map<int, data> m; //!< Energy distributions
    Tspace &spc;

    double g2g(const Tgroup &g1, const Tgroup &g2); //<! exact ion-ion energy between particles
    void save() const;                              //!< save to disk
    void _sample() override;
    void _to_json(json &j) const override;

  public:
    MultipoleDistribution(const json &j, Tspace &spc);
    ~MultipoleDistribution();

}; // end of multipole distribution

/** @brief Sample scattering intensity */
class ScatteringFunction : public Analysisbase {
    Tspace &spc;
    bool usecom;                    // scatter from mass center, only?
    std::string filename;           // output file name
    std::vector<Point> p;           // vector of scattering points
    std::vector<int> ids;           // Molecule ids
    std::vector<std::string> names; // Molecule names
    typedef Scatter::FormFactorUnity<double> Tformfactor;
    Scatter::DebyeFormula<Tformfactor> debye;

    void _sample() override;
    void _to_json(json &j) const override;

  public:
    ScatteringFunction(const json &j, Tspace &spc);
    ~ScatteringFunction();
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
    Tspace &spc;
    std::map<int, Average<double>> Rg2, Rg, Re2, Re, Rs, Rs2, Rg2x, Rg2y, Rg2z;
    std::vector<int> ids; // molecule id's to analyse

    void _to_json(json &j) const override;
    Point vectorgyrationRadiusSquared(typename Tspace::Tgroup &g) const;
    void _sample() override;

  public:
    PolymerShape(const json &j, Tspace &spc);
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

  public:
    QRtraj(const json &j, Tspace &spc);
};

struct CombinedAnalysis : public BasePointerVector<Analysisbase> {
    CombinedAnalysis(const json &j, Tspace &spc, Energy::Hamiltonian &pot);
    void sample();
    ~CombinedAnalysis();
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
