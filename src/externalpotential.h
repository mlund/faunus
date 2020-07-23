#pragma once

#include "group.h"
#include "aux/timers.h"
#include "aux/equidistant_table.h"
#include <set>

template<typename T> class ExprFunction;

namespace Faunus {

struct Change;
class Space;

namespace Energy {

/**
 * All energies inherit from this class
 */
class Energybase {
  public:
    enum keys { ACCEPTED_MONTE_CARLO_STATE, TRIAL_MONTE_CARLO_STATE, NONE };
    keys key = NONE;
    std::string name;                                     //!< Meaningful name
    std::string citation_information;                     //!< Possible reference. May be left empty
    TimeRelativeOfTotal<std::chrono::microseconds> timer; //!< Timer for measure speed of each term
    virtual double energy(Change &) = 0;                  //!< energy due to change
    virtual void to_json(json &) const;                   //!< json output
    virtual void sync(Energybase *, Change &);
    virtual void init();                                  //!< reset and initialize
    virtual inline void force(PointVector &){};           //!< update forces on all particles
    inline virtual ~Energybase() = default;
};

void to_json(json &j, const Energybase &base); //!< Converts any energy class to json object

/**
 * @brief Base class for external potentials
 *
 * Apply an external energy to a defined list of molecules, either acting on individual
 * atoms or the mass-center. The specific energy function, `externalPotentialFunc`
 * is defined in derived classes.
 *
 * @todo The `dN` check is inefficient as it calculates the external potential on *all* particles.
 */
class ExternalPotential : public Energybase {
  private:
    bool act_on_mass_center = false;                   //!< apply only on center-of-mass
    std::set<int> molecule_ids;                        //!< ids of molecules to act on
    std::vector<std::string> molecule_names;           //!< corresponding names of molecules to act on
    double groupEnergy(const Group<Particle> &) const; //!< external potential on a single group
  protected:
    Space &space;                                                            //!< reference to simulation space
    std::function<double(const Particle &)> externalPotentialFunc = nullptr; //!< energy of single particle
  public:
    ExternalPotential(const json &, Space &);
    double energy(Change &) override;
    void to_json(json &) const override;
};

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
std::function<double(const Particle &)> createGouyChapmanPotential(const json &j, const Geometry::Chameleon &);

/**
 * @brief Custom external potential on molecules
 */
class CustomExternal : public ExternalPotential {
  private:
    std::unique_ptr<ExprFunction<double>> expr;
    struct ParticleData { // storage for particle properties
        double charge = 0, x = 0, y = 0, z = 0;
    };
    ParticleData particle_data;
    json json_input_backup; // initial json input

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
    std::string filename;                            //!< input/output filename to charge density profile
    bool fixed_potential = false;                    //!< If true, the potential is never updated
    unsigned int nstep = 0;                          //!< Internal between samples
    unsigned int phi_update_interval = 0;            //!< Distance between phi updating
    unsigned int num_rho_updates = 0;                //!< Number of time rho has been updated
    unsigned int num_density_updates = 0;            //!< Number of charge density updates
    double dielectric_constant;                      //!< Relative dielectric constant
    double dz;                                       //!< z spacing between slits (A)
    double bjerrum_length;                           //!< Bjerrum length (A)
    double half_box_length_z;                        //!< Half box length in z direction
    Equidistant2DTable<double> charge_profile;       //!< instantaneous charge as func. of z
    Equidistant2DTable<double, Average<double>> rho; //!< Charge density at z (unit A^-2)
    Equidistant2DTable<double> phi;                  //!< External potential at z (unit: beta*e)

    double phi_ext(double, double) const; //!< Calculate external potential
    void update_rho();                    //!< update average charge density
    void update_phi();                    //!< update average external potential
    void save_rho();                      //!< save charge density profile to disk
    void load_rho();                      //!< load charge density profile from disk
    void to_json(json &) const override;
    void sync(Energybase *, Change &) override;

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
 * The constructor requires a functor that operates on each particle
 * and returns the resulting energy (kT)
 */
class ParticleSelfEnergy : public ExternalPotential {
  public:
    ParticleSelfEnergy(Space &, std::function<double(const Particle &)>);
};

} // namespace Energy
} // namespace Faunus
