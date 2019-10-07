#pragma once

#include "group.h"
#include "auxiliary.h"
#include <set>

template<typename T> class ExprFunction;

namespace Faunus {

struct Change;
struct Space;

namespace Energy {

    /**
 * All energies inherit from this class
 */
class Energybase {
  public:
    enum keys { OLD, NEW, NONE };
    keys key = NONE;
    std::string name;                                     //!< Meaningful name
    std::string cite;                                     //!< Possible reference. May be left empty
    TimeRelativeOfTotal<std::chrono::microseconds> timer; //!< Timer for measure speed of each term
    virtual double energy(Change &) = 0;                  //!< energy due to change
    virtual void to_json(json &j) const;                  //!< json output
    virtual void sync(Energybase *, Change &);
    virtual void init();                               //!< reset and initialize
    virtual inline void force(std::vector<Point> &){}; // update forces on all particles
    inline virtual ~Energybase(){};
};

void to_json(json &j, const Energybase &base); //!< Converts any energy class to json object

/**
 * @brief Base class for external potentials
 *
 * This will apply an external energy to a defined
 * list of molecules, either acting on individual
 * atoms or the mass-center. The specific energy
 * function, `func` is injected in derived classes.
 */
class ExternalPotential : public Energybase {
  protected:
    typedef typename Faunus::ParticleVector Tpvec;
    bool COM = false; // apply on center-of-mass
    Space &spc;
    std::set<int> molids;                                   // molecules to act upon
    std::function<double(const Particle &)> func = nullptr; // energy of single particle
    std::vector<std::string> _names;

    double _energy(const Group<Particle> &) const; //!< External potential on a single particle
  public:
    ExternalPotential(const json &, Space &);

    /*
     * @todo The `dN` check is very inefficient
     * as it calculates the external potential on *all*
     * particles.
     */
    double energy(Change &) override;
    void to_json(json &) const override;
}; //!< Base class for external potentials, acting on particles

/**
 * @brief Returns a functor for a Gouy-Chapman electric potential
 * @details Gouy-Chapman equations:
 * @f[ \rho = \sqrt{\frac{2 c_0}{\pi l_B} } \sinh ( \beta \phi_0 e / 2 ) @f]
 * @f[ \beta e \phi_0 = 2\mbox{asinh} \left ( \rho \sqrt{\frac{\pi \lambda_B} {2 c_0}} \right ) @f]
 * @f[ \Gamma_0=\tanh{ \beta \phi_0 z e / 4 } @f]
 * where `lB` is the Bjerrum length, `kappa` the inverse Debye length, and `c_0` the
 * bulk salt concentration.
 *
 * @warning: under construction
 */
std::function<double(const Particle &)> createGouyChapmanPotential(const json &);

/**
 * @brief Custom external potential on molecules
 */
class CustomExternal : public ExternalPotential {
  private:
    std::unique_ptr<ExprFunction<double>> expr;
    struct Data { // variables
        double q = 0, x = 0, y = 0, z = 0;
    };
    Data d;
    json jin; // initial json input

  public:
    CustomExternal(const json &, Space &);
    void to_json(json &) const override;
};

/**
 * @brief Mean field electric potential from outside rectangular simulation box.
 * @date Asljunga, December 2010.
 *
 * Calculates the average potential outside a simulation box due to ion
 * densities inside the box, extended to infinity.
 * The update() function calculates the average charge densities in slits in the xy-plane
 * and with a spacing dz. This is used to evaluate the electric potential along
 * the z-axis by using the method by Torbjorn and CO (Mol. Phys. 1996, 87:407). To avoid
 * energy drifts, update() returns the energy change brought about by updating the charge profile.
 */
class ExternalAkesson : public ExternalPotential {
  private:
    std::string filename; //!< File name for average charge
    bool fixed;
    unsigned int nstep = 0;       //!< Internal between samples
    unsigned int nphi = 0;        //!< Distance between phi updating
    unsigned int updatecnt = 0;   //!< Number of time rho has been updated
    double epsr;                  //!< Relative dielectric constant
    double dz;                    //!< z spacing between slits (A)
    double lB;                    //!< Bjerrum length (A)
    double halfz;                 //!< Half box length in z direction
    Equidistant2DTable<double> Q; //!< instantaneous net charge

  public:
    unsigned int cnt = 0;                            //!< Number of charge density updates
    Equidistant2DTable<double, Average<double>> rho; //!< Charge density at z (unit A^-2)
    Equidistant2DTable<double> phi;                  //!< External potential at z (unit: beta*e)

  private:
    void to_json(json &) const override;
    void save();
    void load();

    /*
     * This is Eq. 15 of the mol. phys. 1996 paper by Greberg et al.
     * (sign typo in manuscript: phi^infty(z) should be "-2*pi*z" on page 413, middle)
     */
    double phi_ext(double, double) const;
    void sync(Energybase *, Change &) override;

    // update average charge density
    void update_rho();

    // update average external potential
    void update_phi();

  public:
    ExternalAkesson(const json &, Space &);
    double energy(Change &) override;
    ~ExternalAkesson();
};

/**
 * @brief Confines molecules inside geometric shapes
 */
class Confine : public ExternalPotential {
  public:
    enum Variant { sphere, cylinder, cuboid, none };
    Variant type = none;

  private:
    Point origo = {0, 0, 0}, dir = {1, 1, 1};
    Point low, high;
    double radius, k;
    bool scale = false;
    std::map<std::string, Variant> m = {{"sphere", sphere}, {"cylinder", cylinder}, {"cuboid", cuboid}};

  public:
    Confine(const json &, Space &);
    void to_json(json &) const override;
}; //!< Confine particles to a sub-region of the simulation container

/**
 * @brief Used to add self energies to atoms from i.e. electrostatic schemes
 *
 * The energy class "Nonbonded" detects if a pair-potential carries a self-energy
 * and will automatically inject an instance of ParticleSelfEnergy into the
 * Hamiltonian.
 *
 * The constructor requires a functor that operators on each particle
 * and returns the resulting energy (kT)
 */
class ParticleSelfEnergy : public ExternalPotential {
  public:
    ParticleSelfEnergy(Space &, std::function<double(const Particle &)>);
};

} // namespace Energy
} // namespace Faunus
