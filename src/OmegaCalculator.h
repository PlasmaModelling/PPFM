 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef OMEGA_CALCULATOR_H
#define OMEGA_CALCULATOR_H 

/**
 * @file OmegaCalculator.h
 * @brief This file contains the OmegaCalculator class hierarchy for computing collision integrals.
 * Collision integrals are essential for calculating transport properties in gas mixtures, 
 * and they can be computed using different methods depending on the interaction type.
 * The OmegaCalculator class defines a common interface for these computations, 
 * and derived classes implement specific algorithms for different interactions.
 * @cite devoto1976c
*/

/* #include"TransportCrossSection.h" */
#include"DataLoader.h"

// forward declarations 
class Species ; 
class ChargedSpecies ;
class InteractionInterface ;
class CsCalculator ; 
class CsHolder ; 
class MultiCs ; 
class Potential ;
class TcsInterface ;  
class Composition ;

/**
 * @class OmegaCalculator
 * @brief Interface class for algorithms to compute on Collision Integrals.
 * A pointer to this base class is possessed by a CollisionIntegral object which 
 * initialize the method polimorphically 
*/
class OmegaCalculator {

    protected:

    /// @brief "SpeciesFormula_SpeciesFormula"
    std::string interactionName ; 

    /// @brief Pointer to specie 1 of the interaction.
    Species* sp1 ; 

    /// @brief Pointer to specie 2 of the interaction.
    Species* sp2 ; 

    /// @brief Base constructor of CsCalculator classes assign species and interactionName.
    OmegaCalculator(Species*sp1,Species*sp2) ;
    
    /// @brief Base constructor of CsCalculator classes assign interactionName.
    OmegaCalculator(InteractionInterface* i) ;

    public:

    /** @brief Interface function to compute the collision integral
     * @param l-th momentum of the TCS
     * @param s-th order of the Sonine polynomials 
     * @param T temperature of collision [K°] 
     * @param lambda debye length [m]
     * @param TcS Pointer to the TransportCrossSection interface to access TCS methods
     * without occurring in circular dependance 
     * @return double value of the collision integral funtion of l,s,T and lambda  */
    virtual double Compute( int l, int s, double T, double lambda, CsCalculator* TcS ) = 0;
    
    /// @brief info on the algorithm for the specific interaction.
    virtual void info() = 0 ;
    
    virtual ~OmegaCalculator() = default;

};

/** 
 * @class LoaderOmega
 * @brief class for the Loading of the collision integrals through raw .txt files.
 * further developments will involve more flexible formats to load on to */
class LoaderOmega : public OmegaCalculator, public DataLoader {

    /// @brief Termperature grid [K] for the loaded collision integral data.
    std::vector<double> temperatures;

    /// @brief Matrix of collision integral values for different (l,s) l=1..4, s = 1..7, as a function of temperature.
    std::vector<std::vector<double>> qColumns;

    /// @brief Initialize the loader. Realization of DataLoader method. 
    void Init() override { if(!loaded) LoadData("Collision_Integrals",interactionName); }

    ~LoaderOmega(){};

    public:

    /** @brief construct a Loader object searching for the raw file named as the interaction
     * (convention for this program for raw computed 4-th order 
     * collision integrals .txt datafiles) and storing data on this class members */
    LoaderOmega( InteractionInterface* i ) ;

    /** @brief construct a Loader object with custom prefix searching for the raw file named as the interaction
     * (convention for this program for raw computed 4-th order 
     * collision integrals .txt datafiles) and storing data on this class members */
    LoaderOmega( const std::string& prefix , InteractionInterface* i ) ;

    /** @brief compute Omega^(l,s)(T) interpolating with a positive cubic Spline the loaded data 
     * @see OmegaCalculator for info on the parameters */
    double Compute(int l, int s, double temperatura, double lambda, CsCalculator* TcS) override ;

    /// @brief info on the object.
    void info() override { std::cout << "Raw loader"; }

    /// @brief Realization of DataLoader method to build the file name for loading.
    std::string BuildFileName(const std::string& name) override {
        return "Q4th_"+customPrefix+name+".txt";
    }
    
    /// @brief Realization of DataLoader method to parse the loaded file and store data in class members.
    void ParseFile(std::ifstream& file) override;
    
};


/** 
 * @class CoulombOmega
 * @brief class for the calculation of ChargedSpecie-ChargedSpecie (Coulomb)
 * collision integral. 
 * @details Coulomb collision integrals only depend on the charge state and the mass 
 * of the interacting species through the shielded Coulomb potential. 
 * As so, they've been computed several times and are interpolated in PPFM 
 * from Devoto developments. */
class CoulombOmega : public OmegaCalculator {

    /** @brief Matrix that stores computed Coulomb collision integrals on adimensional 
     * temperatures (the first row) for attractive (even rows) and repulsive (odd rows) 
     * shielded coulomb potential, (l,s) combinations necessary for the 4-th order approximation
     * to transport properties. @see CoulombOmega::CoulombOmega for more. */
    std::vector<std::vector<double>> qcTxt ;

    public:

    /// @brief Constructor to CoulombOmega class, initialize Qctxt hard-coded data.
    CoulombOmega(ChargedSpecies* s1, ChargedSpecies* s2) ;

    /// @brief Constructor to CoulombOmega class, initialize Qctxt hard-coded data. 
    CoulombOmega(InteractionInterface* i ) ; 

    /// @brief Print info when called.    
    virtual void info() override {std::cout<<"Coulomb collision"; } ;

    /** @brief compute Omega^(l,s)(T) for Coulomb collision, 
     * @see OmegaCalculator for info on the parameters */
    double Compute( int l, int s, double T, double Lam, CsCalculator* TcS ) override ;

};

/** 
 * @class NonCoulombOmega
 * @brief class for the calculation of Non-Coulomb collision integrals.
 * @details Non-Coulomb collision integrals depend on the 
 * potential between the interacting species.
 * Also, this class's a reference for others, 
 * with the compute being called whenever integration of TCS is required.   
*/
class NonCoulombOmega : public OmegaCalculator {

    protected:

    /// @brief 32th order Laguerre polynomial nodes to integrate Transport Cross Section. 
    std::vector<double> x;

    /// @brief 32th order Laguerre polynomial weights.
    std::vector<double> w;

    /// @brief Compute Omega^(l,s) from a CsHolder object.
    double ModuleCompute(int l, int s, double T, double Lam, CsHolder* TcS);

    /// @brief Compute Omega^(l,s) from a standard CsCalculator via interpolation.
    double IntegrateOmega(int l, int s, double T, double Lam, CsCalculator* TcS);

    /// @brief Compute degeneracy-weighted Omega^(l,s) for MultiCs calculators.
    double MultiCompute(int l, int s, double T, double Lam, MultiCs* p);

    public:

    /// @brief Initialize x and w (Laguerre nodes and weights).
    NonCoulombOmega(InteractionInterface* i);

    /// @brief Initialize x and w, and assign ChiIntegrator based on Potential type.
    NonCoulombOmega(InteractionInterface* i, TcsInterface* t, Potential* pot);

    /// @brief Print info when called.
    virtual void info() override { std::cout << "Non-coulomb coll."; }

    /// @brief Compute Omega^(l,s)(T) based on the type of calculator.
    double Compute(int l, int s, double T, double Lam, CsCalculator* TcS) override;

};

/**
 * @class ChargeExchangeOmega
 * @brief Charge Exchange collision integrals class called for simplified setting
 * of the charge-exchange cross section. It copies the elastic interaction already set
 * and initialize the charge exchange cross section.
 * @cite GhediniColomboSanibondi
*/
class ChargeExchangeOmega : public OmegaCalculator {

    /** @brief Inner elastic collision calculator copied from the Collision Integral
     * object during contruction. */
    OmegaCalculator* omegaEl ; 
  
    /// @brief Coefficients A and B for the analytical formula of the charge exchange collision integral.
    double A,B;

    public: 
    
    /// @brief Constructor for ChargeExchangeOmega, initializes the elastic part and sets A and B coefficients. 
    ChargeExchangeOmega( InteractionInterface* i , OmegaCalculator* OMel, double A , double B ) ;

    /** @brief Compute the elastic part with the previous calculator and the inelastic
     * part with the analytical formula 12 of Devoto R S 1967 Phys. Fluids 10 354. 
     * for the charge exchange collision integral. 
     * @return the square modulus of Qel and Qe for odd l values . 
     * @see OmegaCalculator for info on the parameters */
    double Compute ( int l, int s, double temperatura, double lambda, CsCalculator* TcS ) override ;

    /// @brief Print info when called 
    void info() override { omegaEl->info(); std::cout << " & Charge Exchange " ; }

};

#endif