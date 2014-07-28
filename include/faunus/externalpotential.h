#ifndef FAUNUS_EXTERNALPOTENTIAL_H
#define FAUNUS_EXTERNALPOTENTIAL_H

#include <faunus/potentials.h>
#include <faunus/analysis.h>

namespace Faunus {

  namespace Potential {

    /**
     * @brief Base class for external potentials acting on particles
     *
     * If the external potential depends on the particle position
     * one may specify a functor that transforms the particle
     * position into another coordinate. The default signature of this
     * function object is
     *
     * - `double(const Point&)`
     *
     * and is set with `setCoordinateFunc`.
     * The signature can be changed by the
     * template argument, `TcoordFunc`.
     *
     * @tparam TcoordFunc Type of optional functor to transform particle
     *                    position into internal coordinate.
     */
    template<class TcoordFunc=std::function<double(const Point&)> >
      class ExternalPotentialBase {
        private:
          virtual string _info()=0;
        protected:
          TcoordFunc p2c; // function to convert point to internal coordinate
        public:
          ExternalPotentialBase() : p2c(nullptr) {}

          virtual ~ExternalPotentialBase() {}

          template<typename T> void setCoordinateFunc(T f) { p2c=f; }

          string info() {
            return _info();
          }

          template<class Tparticle>
            Point field(const Tparticle &) { return Point(0,0,0); }

          std::string name;
      };

    /**
     * @brief Excess chemical potential for ions using Debye-Huckel theory
     */
    template<class T=double>
      class ExcessDH : public ExternalPotentialBase<> {
        private:
          std::string _info() { return _dh.info(20); }
          Potential::DebyeHuckel _dh;
        public:
          ExcessDH(Faunus::InputMap &in) : _dh(in) {
            name="Debye-Huckel Single Ion Excess";
          }
          template<typename Tparticle> T operator()(const Tparticle &p) {
            return _dh.excessChemPot(p.charge, p.radius);
          }
      };

    /**
     * @brief Gouy-Chapman potential
     *
     * This is an external potential due to a charged Gouy-Chapman surface
     *
     * During construction, the `InputMap` is searched for the following keywords:
     *
     * Keyword                 | Description
     * ----------------------- | ------------------------------------------------
     * `dh_ionicstrength`      | Ionic strength [mol/l]
     * `gouychapman_phi0`      | Surface potential [unitless, i.e. phi_0*e/kT]
     * `gouychapman_qarea`     | Surface charge density (if phi0 not defined)
     * `gouychapman_rho`       | Surface charge density [1/A^2] (if qarea not defined)
     * `gouychapman_offset`    | Shift GC surface [A] (default: 0)
     *
     * Code example:
     * 
     *     typedef Potential::GouyChapman<> Txp;
     *     InputMap in("input");
     *     Cuboidslit geo(in);
     *     Energy::ExternalPotential<Txp> pot(in);
     *     pot.expot.setSurfPositionZ( &geo.len_half.z() );
     *
     * @note Salt is assumed monovalent
     * @date Lund/Asljunga, 2011-2012
     */
    template<class T=double, bool linearize=false>
      class GouyChapman : public ExternalPotentialBase<> {
        private:
          Potential::DebyeHuckel dh;
          T c0;          //!< Ion concentration (1/A^3)
          T rho;         //!< Surface charge density (e/A^2)
          T phi0;        //!< Unitless surface potential @f\frac{\phi0 e}{kT}@f
          T gamma0;      //!< Gouy-chapman coefficient ()
          T lB;          //!< Bjerrum length (A)
          T k;           //!< Inv. debye len (1/A)
          T offset;      //!< Distance off set for hiding GC surface behind the box surface
          std::string _info();
        public:
          GouyChapman(InputMap&);          //!< Constructor
          void setSurfPositionZ(T*);       //!< Set surface position on z-axis
          T surfDist(const Point&);        //!< Point<->GC surface distance
          template<typename Tparticle>
            T operator()(const Tparticle&);//!< Particle<->GC interaction energy
      };

    /**
     * @details Gouy-Chapman equations:
     * @f[ \rho = \sqrt{\frac{2 c_0}{\pi l_B} } \sinh ( \beta \phi_0 e / 2 ) @f]
     * @f[ \beta e \phi_0 = 2\mbox{asinh} \left ( \rho \sqrt{\frac{\pi \lambda_B} {2 c_0}} \right ) @f]
     * @f[ \Gamma_0=\tanh{ \beta \phi_0 z e / 4 } @f]
     * where `lB` is the Bjerrum length, `kappa` the inverse Debye length, and `c_0` the
     * bulk salt concentration.
     */
    template<class T, bool linearize>
      GouyChapman<T,linearize>::GouyChapman(InputMap &in) : dh(in) {
        string prefix = "gouychapman_";
        name = "Gouy-Chapman";
        c0 = dh.ionicStrength() * pc::Nav / 1e27; // assuming 1:1 salt, so c0=I
        lB = dh.bjerrumLength();
        k  = 1/dh.debyeLength();
        //linearize = in(prefix+"linearize",false);
        phi0 = in(prefix+"phi0",0.); // Unitless potential = beta*e*phi0
        if ( std::abs(phi0)>1e-6 )
          rho = sqrt(2*c0/(pc::pi*lB))*sinh(.5*phi0); //Evans&Wennerstrom,Colloidal Domain p 138-140
        else {
          rho = 1/in(prefix+"qarea",0.);
          if (rho>1e9)
            rho = in(prefix+"rho",0.);
          phi0 = 2.*asinh(rho * sqrt(.5*lB*pc::pi/c0 ));//[Evans..]
        }
        gamma0 = tanh(phi0/4); // assuming z=1  [Evans..]
        offset = in(prefix+"offset",0.);
      }

    /**
     * Before using this function make sure to set a surface calculations
     * method either with `setSurfPositionZ` or `setSurfPositionZ`.
     */
    template<class T, bool linearize>
      T GouyChapman<T,linearize>::surfDist(const Point &p) {
        assert(p2c!=nullptr && "Surface distance function not set!");
        return p2c(p);
      }

    template<class T, bool linearize>
      void GouyChapman<T,linearize>::setSurfPositionZ(T* z) {
        setCoordinateFunc
          (
           [=](const Point &p) { return std::abs(*z-p.z()) + offset; }
          ); // c++11 lambda
      }

    /**
     * @details Interaction of charged particle with GC potential:
     * @f[ \beta e \Phi(r_i) = 2\ln{\frac{1+\Gamma_0 \exp{-\kappa r_i}}{1-\Gamma_0 \exp{-\kappa r_i}}}@f]
     * @f[ \beta u = z_i \cdot \beta e \Phi(r_i) @f]
     * where `z_i` is the charge number and `r_i` the distance from the surface.
     */
    template<class T, bool linearize>
      template<typename Tparticle>
      T GouyChapman<T,linearize>::operator()(const Tparticle &p) {
        if (p.charge!=0) {
#ifdef FAU_APPROXMATH
          T x=exp_cawley(-k*surfDist(p));
#else
          T x=exp(-k*surfDist(p));
#endif
          if (linearize)
            return p.charge*phi0*x;
          else {
            x=gamma0*x;
            return 2*p.charge*log( (1+x)/(1-x) );
          }
        }
        return 0;
      }

    template<class T, bool linearize>
      std::string GouyChapman<T,linearize>::_info() {
        using namespace textio;
        char w=30;
        std::ostringstream o;
        o << dh.info(w)
          << pad(SUB,w,"Surface potential") << phi0 << kT+"/e = "
          << phi0*pc::kB*pc::T()/pc::e << " V=J/C" << endl
          << pad(SUB,w,"Surface charge density") << rho*pc::e*1e20 << " C/m"+squared << endl
          << pad(SUB,w,"Area per charge") << 1/rho << _angstrom+squared << endl
          << pad(SUB,w,"GC-coefficient "+Gamma+"o") << " " << gamma0 << endl
          << pad(SUB,w,"Linearize") << ((linearize) ? "yes" : "no") << endl
          << pad(SUB,w,"Surface position") << surfDist(Point(0,0,0)) << _angstrom << endl;
        return o.str();
      }

    /**
     * @brief Mean field correction
     *
     * This external potential will add a mean field correction
     * for electrostatics outside a cylindrical cut-off.
     * 
     * @author Anil Kurut
     * @warning untested!
     * @todo test
     */
    template<typename T=double>
      class CylindricalCorrectionDH : public ExternalPotentialBase<> {
        private:
          typedef Analysis::Table2D<T,Average<T> > Ttable;
          T threshold;       //!< Threshold for mean field; must be equal to radius of cylinderical container
          T bin;             //!< Resolution for the container slices 
          T prefactor;       //!< exp(-kappa*threshold)
          Ttable qdensity;   //!< Tabulated charge desity for each slice of the container.
          string filename;   //!< Filename of charge density file
          bool loadfromdisk; //!< `True` = Load from disk, `False` = Sample charge density
          string _info();
        public:
          CylindricalCorrectionDH(InputMap&, std::string="mfc_");  //!< Constructor
          template<class Tparticle> T operator()(const Tparticle&);//!< External potential on particle
          template<class Tpvec> void sample(const Tpvec&, T, T);   //!< Sample charge density
      };

    template<class T>
      template<class Tparticle>
      T CylindricalCorrectionDH<T>::operator()(const Tparticle &p) {
        return prefactor*p.charge*qdensity( p.z() );
      }

    /**
     * Sampling is done only when `loadfromdisk==false`.
     * After each sampling event, the charge
     * density table is saved to disk.
     *
     * @param p Particle vector
     * @param zmin Minimum z distance to sample
     * @param zmax Maximum z distance to sample
     */
    template<class T>
      template<class Tpvec>
      void CylindricalCorrectionDH<T>::sample(const Tpvec& p, T zmin, T zmax){
        if (!loadfromdisk) {
          typedef Analysis::Table2D<T,T> Ttable;
          Ttable qsum(bin,Ttable::XYDATA);   //total bin charge
          T dV=pc::pi*pow(threshold, 2)*bin; // bin volume
          for (auto &i : p)
            qsum(i.z())+=i.charge;
          for (T z=zmin; z<=zmax; z+=bin)
            qdensity(z)+=qsum(z)/dV;
          qdensity.save(filename);
        }
      }

    /**
     * @brief Construct from InputMap
     *
     * The addtion to `Potential::DebyeHuckel` keywords, we look for:
     *
     * Keyword       | Description
     * ------------- | ----------------------------------------
     * mfc_load      | `True` load from disk (default: `false`)
     * mfc_filename  | Charge density to load/save
     * mfc_radius    | Radius of cut-off
     */
    template<class T>
      CylindricalCorrectionDH<T>::CylindricalCorrectionDH(InputMap& in, string pfx) :
        bin(in.get<T>("CylindricalCorrectionDH_binsize", 2)),
        qdensity(bin,Ttable::XYDATA) {
          name = "Cylindrical DH Correction";  
          Potential::DebyeHuckel dh(in);
          threshold=in(pfx+"radius", pc::infty);
          loadfromdisk=in(pfx+"load", false);
          filename=textio::prefix+pfx+"qdensity";
          prefactor=std::exp( -threshold/dh.debyeLength() )
            * dh.bjerrumLength()*pc::pi*2*bin*dh.debyeLength();
          if (loadfromdisk)
            qdensity.load(filename);
        }

    template<class T>
      string CylindricalCorrectionDH<T>::_info() {
        using namespace textio;
        char w=30;
        std::ostringstream o;
        o << pad(SUB,w,"Mean Field hole radius") << threshold << _angstrom << endl
          << pad(SUB,w,"Mean Field bin width") << bin << _angstrom << endl
          << pad(SUB,w,"Prefactor") << prefactor << _angstrom+cubed << endl;
        return o.str();
      }

    /**
     * @brief External, uniform electric field
     *
     * This is a "potential" for `ExternalPotential` to
     * mimic a uniform electric field. The template argument must
     * be a pair potential that implements a `fieldEnergy` function.
     * The `InputMap` is searched for `field_x`, `field_y`, and
     * `field_z` which indicate the electric field vector. Default
     * value is `(0,0,0)`.
     *
     * Example:
     *
     * ~~~~
     * typedef Potential::ElectricField<Potential::DipoleDipole> Texpot;
     * ExternalPotential<Tspace,Texpot> pot(in);
     * ~~~~
     *
     * @todo unifinished and lacking doc (units?)
     */
    template<class Tpairpot>
      class ElectricField : ExternalPotentialBase<> {
        private:
          typedef ExternalPotentialBase<> base;
          Tpairpot pairpot;
          Point E; // electric field vector

          string _info() { return pairpot.info(); }

        public:
          ElectricField(InputMap &in) : pairpot(in) {
            base::name = "Uniform electric field";
            E.x() = in.get<double>("field_x", 0, "Electric field (x)");
            E.y() = in.get<double>("field_y", 0, "Electric field (y)");
            E.z() = in.get<double>("field_z", 0, "Electric field (z)");
          }

          /** @brief Interaction of particle `p` with electric field */
          template<class Tparticle>
            double operator()(const Tparticle &p) {
              return pairpot.fieldEnergy(p,E);
            }

          /** @brief Electric field on particle `p` */
          template<class Tparticle>
            Point field(const Tparticle &p) { return E; }
      };

    /**
     * @brief Hydrophobic wall potential
     * @author Joao Henriques
     * @date Lund, 2014
     *
     * This (external) potential class is used to simulate hydrophobic interactions between 
     * particle(s) and a surface, using a simple square well (the default) or a Lennard-Jones
     * potential. Surface position must be specified in the program even if one has already
     * done it for the Gouy-Chapman potential (both classes inherit from ExternalPotentialBase<>
     * but are for the most part independent).
     *
     * See doi:10.1021/la300892p for more details on the square well potential.
     *
     * The Lennard-Jones potential has the form:
     * @f$
     * \beta u=\epsilon
     * \left ((\sigma_{i}/r_{i,s})^{12}-2(\sigma_{i}/r_{i,s})^6\right )
     * @f$
     * where
     * \f$\sigma_{i}\f$ is the residue/particle radius. The potential reaches its minimum when
     * \f$r_{i,s} = \sigma_{i}\f$, ie. the residue/particle is in close contact with the wall.
     *   
     * The InputMap parameters are:
     *
     * Key                  | Description
     * :------------------- | :---------------------------
     * `hydrwl_type`        | Type of potential, ie. square well ("sqwl", default) or Lennard-Jones ("lj") 
     * `hydrwl_depth`       | Depth, \f$\epsilon\f$ [kT] (positive number)
     * `hydrwl_threshold    | Threshold, [angstrom] (particle center-to-wall distance) - for "sqwl" type only!
     *
     */
    template<class T=double>
      class HydrophobicWall : public ExternalPotentialBase<> {
        protected:
          T _depth;
          T _threshold;
          std::string _type;
          std::string _info();
        public:
          HydrophobicWall(InputMap&);
          void setSurfPositionZ(T*);        // sets position of surface
          template<typename Tparticle>
            T operator()(const Tparticle &p); // returns energy
      };

    template<class T>
      HydrophobicWall<T>::HydrophobicWall(InputMap &in) {
        string prefix = "hydrwl_";
        name          = "Hydrophobic Wall";
        _type         = in.get<string>(prefix + "type", "sqwl");
        _depth        = in.get<double>(prefix + "depth"    , 0);
        _threshold    = in.get<double>(prefix + "threshold", 0);
      }

    template<class T>
      void HydrophobicWall<T>::setSurfPositionZ(T* z) {
        this->setCoordinateFunc
          (
           [=](const Point &p) { return std::abs(*z-p.z()); }
          ); // c++11 lambda
      }

    template<class T>
      template<typename Tparticle>
      T HydrophobicWall<T>::operator()(const Tparticle &p) {
        if (p.hydrophobic) {
          if (_type == "sqwl")
            if (this->p2c(p) < _threshold)
              return -_depth;
          if (_type == "lj") {
            double r2  = (p.radius * p.radius) / (this->p2c(p) * this->p2c(p));
            double r6  = r2 * r2 * r2;
            double val = _depth * ((r6 * r6) - (2 * r6));
            return val;
          } 
        }
        return 0;
      }

    template<class T>
      std::string HydrophobicWall<T>::_info() {
        std::ostringstream o;
        if (_type == "sqwl")
          o << pad(textio::SUB, 50, ">>> USING: square well potential <<<") << endl
            << pad(textio::SUB, 26, "Depth, " + textio::epsilon + "(SQWL)") << _depth << textio::kT + " = " << pc::kT2kJ(_depth) << " kJ/mol" << endl
            << pad(textio::SUB, 25, "Threshold") << _threshold << textio::_angstrom << " (particle - wall distance)" << endl;
        if (_type == "lj")
          o << pad(textio::SUB, 50, ">>> USING: Lennard-Jones potential <<<") << endl
            << pad(textio::SUB, 26, "Depth, " + textio::epsilon + "(LJ)") << _depth << textio::kT + " = " << pc::kT2kJ(_depth) << " kJ/mol"<< endl;
        return o.str();
      }

    /**
     * @brief Hydrophobic wall potential w. linear slope
     * @date Lund, 2014
     *
     * As `HydrophobicWall` but the potential is varying linearly with
     * distance. It the potential is zero at `threshold` and `depth` when
     * the particle-surface separation is zero.
     */
    template<class T=double>
      struct HydrophobicWallLinear : public HydrophobicWall<T> {
        HydrophobicWallLinear(InputMap &in) : HydrophobicWall<T>::HydrophobicWall(in) {
          this->name+=" Linear";
        }
        template<typename Tparticle>
          T operator()(const Tparticle &p) {
            if (p.hydrophobic) {
              double d = this->p2c(p);
              assert(d>0 && "Particle-surface distance must be positive");
              if (d<this->_threshold)
                return -(this->_depth) * (1-d/this->_threshold);
            }
            return 0;
          }
      };

  } //namespace
} //namespace
#endif

