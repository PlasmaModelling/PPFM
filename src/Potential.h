#ifndef POTENTIAL_H
#define POTENTIAL_H 

#include "Species.h"

/// @brief Abstract base class for pair interaction potentials.
/// @details Implementations must return the potential energy V(r) given the
///          inter-particle separation r.
/// @note Input: r in Å (unless a specific class states otherwise).
///       Output: V(r) in eV.
class Potential {

    public:
    
    /// @brief Evaluate the potential at distance r.
    /// @param r distance [Å]
    /// @return potential energy [eV]
    virtual double Pot(double r) = 0;
};

/// @brief Debye–Hückel (screened Coulomb) potential.
/// @details V(r) ~ (z1 z2 e^2 / (4π ε0 r)) * exp(-r/λ_D).
/// @note Input r in Å; output in eV (assuming internal constants handle units).
class ShieldedCoulombPotential : public Potential {
    double z1;      ///< charge number of particle 1 [-]
    double z2;      ///< charge number of particle 2 [-]
    double lambdaD; ///< Debye length [Å] (or [m] if handled consistently inside Pot)


    public:
    
    /// @brief Construct a screened Coulomb potential.
    /// @param Z1 charge number of particle 1 [-]
    /// @param Z2 charge number of particle 2 [-]
    /// @param LambdaD Debye length [Å] (or consistent unit with Pot)
    ShieldedCoulombPotential(double Z1, double Z2, double LambdaD)
        : z1(Z1), z2(Z2), lambdaD(LambdaD) {}

    /// @copydoc Potential::Pot
    double Pot(double r) override;
};

/// @brief Lennard–Jones (12-6) potential.
/// @details V(r) = 4 ε [ (σ/r)^12 - (σ/r)^6 ].
/// @note Input r in Å; output in eV.
class LennardJones : public Potential {

    protected:
    double eps; ///< ε [eV]
    double sig; ///< σ [Å]


    public:
    
    /// @brief Construct with explicit parameters.
    /// @param eps epsilon [eV]
    /// @param sig sigma [Å]
    LennardJones(double eps, double sig) : eps(eps), sig(sig) {}

    /// @copydoc Potential::Pot
    double Pot(double r) override;

    // Hard-coded constructors for specific pairs (implemented in HCPotentialDB.cpp)
    /// @brief Predefined parameters for Ar–Ar.
    LennardJones(Argon* ar, Argon* arr);
    /// @brief Predefined parameters for He–He.
    LennardJones(Helium* t1, Helium* t2);
    /// @brief Predefined mixed rule for Ar–He (e.g., Lorentz–Berthelot).
    LennardJones(Argon* t1, Helium* t2);
};

/// @brief Capitelli potential.
/// @details V(r) uses a generalized inverse-power repulsion with attractive tail;
///          typical for neutral/charged interactions.
/// @note Input r in Å; output in eV.
class Capitelli : public Potential {

    protected:
    
    Species* t1; ///< species 1
    Species* t2; ///< species 2

    double beta; ///< shape/exponent parameter [-]
    double De;   ///< well depth [eV]
    double Re;   ///< equilibrium distance [Å]


    public:
    
    /// @brief Construct with explicit parameters.
    /// @param sp1 species 1
    /// @param sp2 species 2
    /// @param beta shape parameter [-]
    /// @param De well depth [eV]
    /// @param Re equilibrium distance [Å]
    Capitelli(Species* sp1, Species* sp2, double beta, double De, double Re)
        : beta(beta), De(De), Re(Re) {}

    /// @copydoc Potential::Pot
    double Pot(double r) override;

    // Hard-coded predefined pairs (implemented in HCPotentialDB.cpp)
    Capitelli(Argon* t1, Argon* t2);
    Capitelli(Argon* t1, ArgonI* t2);
    Capitelli(Argon* t1, ArgonII* t2);
    Capitelli(Argon* t1, Nitrogen* t2);
    Capitelli(Argon* t1, NitrogenI* t2);
    Capitelli(Argon* t1, MolecularNitrogen* t2);
    Capitelli(ArgonI* t1, Nitrogen* t2);
    Capitelli(ArgonI* t1, MolecularNitrogen* t2);
    Capitelli(ArgonII* t1, Nitrogen* t2);
    Capitelli(ArgonII* t1, MolecularNitrogen* t2);
};

/// @brief Reduced Hulburt–Hirschfelder potential.
/// @details Classic diatomic potential form with polynomial correction on top of Morse-like terms.
/// @note Input r in Å; output in eV.
class HulburtHirschfelder : public Potential {

    protected:
    
    double De; ///< well depth [eV]
    double Re; ///< equilibrium distance [Å]

    double a; ///< reduced parameter a [-] (often related to curvature at Re)
    double b; ///< reduced cubic correction coefficient [-]
    double c; ///< reduced quartic correction coefficient [-]


    public:
    
    /// @brief Construct with explicit parameters.
    /// @param De well depth [eV]
    /// @param Re equilibrium distance [Å]
    /// @param a reduced parameter a [-]
    /// @param b reduced parameter b [-]
    /// @param c reduced parameter c [-]
    HulburtHirschfelder(double De, double Re, double a, double b, double c)
        : De(De), Re(Re), a(a), b(b), c(c) {}

    /// @copydoc Potential::Pot
    double Pot(double r) override;
};

/// @brief Unreduced Hulburt–Hirschfelder potential.
/// @details Parameterization using spectroscopic constants.
/// @note Input r in Å (or meters—must be consistent with constructor); output in eV.
class HulburtHirschfelderUnreduced : public Potential {

    protected:
    
    double eps;    ///< well depth [eV]
    double re;     ///< equilibrium distance [Å] (or [m], must match usage)
    double we;     ///< vibrational constant [re]^-1
    double weXe;   ///< anharmonicity constant [re]^-1
    double Be;     ///< rotational constant [re]^-1
    double Alphae; ///< rovibrational coupling [re]^-1

    public:
    
    /// @brief Construct with spectroscopic constants.
    /// @param eps well depth [eV]
    /// @param re equilibrium distance [Å] (or [m])
    /// @param we vibrational constant [re]^-1
    /// @param weXe anharmonicity constant [re]^-1
    /// @param Be rotational constant [re]^-1
    /// @param Alphae rovibrational coupling [re]^-1
    HulburtHirschfelderUnreduced(double eps, double re, double we,
                                 double weXe, double Be, double Alphae)
        : eps(eps), re(re), we(we), weXe(weXe), Be(Be), Alphae(Alphae) {}

    /// @copydoc Potential::Pot
    double Pot(double r) override;
};

/// @brief Base class for Morse-type potentials.
class Morse : public Potential {};

/// @brief Power-law repulsive potential V(r) = v0 * r^-n.
/// @note Input r in Å; output in eV.
class PowerPot : public Morse {
    
    double v0; ///< scale factor [eV Å^n]
    double n;  ///< exponent n [-]

    public:
    
    /// @brief Construct a power-law potential.
    /// @param v0 scale factor [eV Å^n]
    /// @param n exponent [-]
    PowerPot(double v0, double n) : v0(v0), n(n) {}

    /// @copydoc Potential::Pot
    double Pot(double r) override;
};

/// @brief Two-parameter Morse exponential form V(r) = V0 * exp(-B r).
/// @note Input r in Å; output in eV. References also refer to this as the classic exponential potential.
class Morse2Param : public Morse {

    protected:
    double V0; ///< V0 [eV]
    double B;  ///< B [Å^-1]


    public:
    
    /// @brief Construct with V0 and B.
    /// @param V0 [eV]
    /// @param B [Å^-1]
    Morse2Param(double V0, double B) : V0(V0), B(B) {}

    /// @copydoc Potential::Pot
    double Pot(double r) override;
};

/// @brief Standard three-parameter Morse potential.
/// @details V(r) = De [e^{-2β(r-Re)} - 2 e^{-β(r-Re)}].
/// @note Input r in Å; output in eV.
class Morse3Param : public Morse {

    protected:

    double De;   ///< well depth [eV]
    double Beta; ///< inverse range [Å^-1]
    double Re;   ///< equilibrium distance [Å]


    public:

    /// @brief Construct classical Morse(De, Beta, Re).
    /// @param De [eV]
    /// @param Beta [Å^-1]
    /// @param Re [Å]
    Morse3Param(double De, double Beta, double Re) : De(De), Re(Re), Beta(Beta) {}

    /// @copydoc Potential::Pot
    double Pot(double r) override;
};

/// @brief Five-parameter generalized Morse potential.
/// @note Input r in Å; output in eV.
class Morse5Param : public Morse {

    double re;     ///< equilibrium distance [Å]
    double D;      ///< well depth [eV]
    double b0;     ///< base slope parameter [- or Å^-1 depending on formulation]
    double gamma;  ///< linear correction coefficient [- or Å^-1 *depending on model*]
    double lambda; ///< quadratic correction coefficient [- or Å^-2 *depending on model*]

    public:

    /// @brief Construct a 5-parameter Morse-like potential.
    /// @param RE re [Å]
    /// @param DD D [eV]
    /// @param B0 b0 [-]
    /// @param GAMMA gamma [-]
    /// @param LAMBDA lambda [-]
    Morse5Param(double RE, double DD, double B0, double GAMMA, double LAMBDA)
        : re(RE), D(DD), b0(B0), gamma(GAMMA), lambda(LAMBDA) {}

    /// @copydoc Potential::Pot
    double Pot(double r) override;
};

/// @brief Charge–induced-dipole polarization potential.
/// @details One species must be charged; the other neutral with polarizability α.
/// @note Input r in Å; output in eV.
class Polarization : public Potential {
    
    std::tuple<Species*, double> inducedDipole; ///< (neutral species, polarizability α [Å^3])
    ChargedSpecies* ion;                         ///< charged species (ion)

    public:
   
    /// @brief Construct polarization interaction between (ion, neutral).
    /// @param sp1 species 1 (ion or neutral)
    /// @param sp2 species 2 (neutral or ion)
    /// @param Polarizability neutral polarizability α [Å^3]
    /// @throws std::invalid_argument if both or neither species are charged
    Polarization(Species* sp1, Species* sp2, double Polarizability);

    /// @copydoc Potential::Pot
    double Pot(double r) override;
};

/// @brief Ad-hoc HFDTCS2 potential for Ar–Ar used to validate Zhang’s study.
/// @details As in Aziz & Slaman, J. Chem. Phys. 92, 1030 (1990). See paper for
///          the piecewise definitions and damping terms.
/// @note Input r in Å; output in eV.
class HFDTCS2_ArAr : public Potential {

    public:
    
    /// @brief Construct the predefined HFDTCS2 Ar–Ar potential.
    HFDTCS2_ArAr();

    /// @copydoc Potential::Pot
    double Pot(double r) override;
};

/// @brief Hartree–Fock Dispersion type-B (HFD-B) potential 
/// Input: r in Å → Output: V(r) in eV.
/// Dimensional form:
/// V(r) = A * exp( -alpha * r + beta * r^2 )
///        - F(x) * ( C6/r^6 + C8/r^8 + C10/r^10 ),  x = r/Rm
///
/// F(x) = exp( -((D/x) - 1)^2 )  for x < D, otherwise 1
class HFD_B : public Potential {
    
    public:
    
    /// @brief Main constructor: all parameters must be given in dimensional units (eV, Å).
    /// @param A [eV]
    /// @param alpha [Å^-1]
    /// @param beta [Å^-2]
    /// @param C6 [eV Å^6]
    /// @param C8 [eV Å^8]
    /// @param C10 [eV Å^10]
    /// @param Rm [Å]
    /// @param D [-]
    HFD_B(double A, double alpha, double beta,
          double C6, double C8, double C10,
          double Rm, double D);

    /// @brief Overloaded constructor with final flag:
    /// If atomic_units == true, all parameters must be provided in atomic units:
    ///   A [Ha], alpha [a0^-1], beta [a0^-2],
    ///   Cn [Ha a0^n], Rm [a0], D dimensionless.
    /// Parameters will be converted internally into eV/Å.
    /// If atomic_units == false, forwards to the main constructor.
    HFD_B(double A, double alpha, double beta,
          double C6, double C8, double C10,
          double Rm, double D, bool atomic_units);

    /// @brief Evaluate potential at distance r.
    /// @param r distance [Å]
    /// @return potential energy [eV]
    double Pot(double r) override;

    private:
    
    void check_dimensional_() const;

    // Stored always in dimensional units (eV / Å)
    double A;     // [eV]
    double alpha; // [Å^-1]
    double beta;  // [Å^-2]
    double C6;    // [eV Å^6]
    double C8;    // [eV Å^8]
    double C10;   // [eV Å^10]
    double Rm;    // [Å]
    double D;     // [-]
};


#endif
