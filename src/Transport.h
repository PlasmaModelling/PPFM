 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef TRANSPORT_H
#define TRANSPORT_H 

/**
 * @file Transport.h
 * @brief This file contains the base class for transport properties calculations.
 * The Transport class defines the common interface and data members 
 * for computing transport properties such as viscosity, thermal conductivity, 
 * and diffusion coefficients. 
 * It also includes a pointer to a CiBox object for accessing collision integrals 
 * needed in the calculations. 
 * Derived classes will implement specific theories (e.g., Chapman-Enskog, Devoto) 
 * for computing these properties based on the gas mixture composition and state.
 * Due to the inner complexity of the formulations involved in the Chapman-Enskog expansion, 
 * some "hiding behind encapsulation" has been needed, 
 * and, also, one of the design choices that demanded the need for PPFM implementation.
*/

#include <vector>
#include <string>
#include <numbers>

// forward declarations
class GasMixture ; 
class CiBox ; 

#define kB 1.380650524e-8 
#define PI std::numbers::pi  
#define q_e 1.60217646e-19

/** @class Transport
 ** @brief Base class common to all Transport Properties theories
 ** @details This class contains members to collect transport properties
 ** and Diffusion Coefficients in SI units. \n
 ** It also serves as a bridge to Collision Integrals computation module \n
 ** through a pointer to a class CiBox, that's used to compute and store \n
 ** Collision Integrals in the matrix Qt through the method QtCalc(). */
class Transport {

    friend class ZhangTpCsv ;
    friend class DevotoTpCsv ;

    protected :

    /// @brief Vector to store computed Transport Properties values
    std::vector<double> Tp ; 

    /// @brief Ordinary Diffusion coefficients
    std::vector<std::vector<double>> D ; 

    /// @brief Thermal Diffusion coefficients
    std::vector<double> DT ; 

    /// @brief Pointer to CollisionIntegral collection @see class CiBox for further details.
    CiBox* Ci ; 

    /** @brief Calculated Collision Integrals stored at every temperature step as a double matrix.
     **  @details the rows span from 0 to 15, the number of (m,p) approximations in the Chapman-Enskog method, \n
     **  the columns span from 0 to Number of Species in the Mixture . \n
     **  default units are [ micron^2 ] */
    std::vector<std::vector<double>> Qt ; 

    /// @brief Constructor for Transport class with assigned CiBox object.
    Transport ( CiBox* cbx ) : Ci(cbx)  {} 

    /// @brief Constructor for Transport class with assigned GasMixture object.
    Transport ( GasMixture* mix ) ;

    /** @brief Compute and store Collision Integrals values 
     ** in Qt via the pointer to CiBox object. \n 
     ** @attention ! Ensure a call to this method at every 
     ** computeTransport implementation ! \n 
     ** @details the default method compute and stores Collision Integrals 
     ** in the adimensional form
     ** @see <a href="../../articles/TransportFormulas.pdf"> TransportFormulas.pdf </a>
     ** eq. 17 */
    void QtCalc ( GasMixture* gasmix ) ;

    public : 

    /** @brief Common Interface to compute Transport properties 
     ** with different theories, \n  
     ** please do refer to implemented modules and place 
     ** computed Transport Coefficients in SI units when implementing. */
    virtual void computeTransport ( GasMixture* gasmix ) = 0 ; 
    
} ;

/**
 ** @class Properties
 ** @brief Class to incapsulate computation methods
 ** @details methods will be unaccesible to user, \n
 ** the user can access to Transport Properties via the methods of the base class Transport \n 
 ** after calling for the computation of them all via the method computeTransport(GasMixture*). \n
 ** Methods incapsulated here will then serve to implement derived classes of Transport Properties computational theories. */
class Properties : public Transport {

    protected :

    /** @brief Vector to store the order of Chapman-Enskog approximations 
     ** for each transport property. To be filled by derived classes constructors 
     ** with default values. Users can use the setOrders method to specify them as a whole. */
    std::vector<int> orders ; 

    /// @brief Constructor for Properties class with assigned CiBox object.
    Properties ( CiBox* cbx  ) : Transport ( cbx ) {}

    /// @brief Constructor for Properties class with assigned GasMixture object.
    Properties ( GasMixture* mix ) : Transport ( mix ) {} ;

    /** @brief Interface for the calculation of the Thermal Conductivity of electrons
     **  @param gasmix GasMixture object 
     **  @param order Desired order of approximation */
    virtual double ThermalCondEl ( GasMixture* gasmix, int order ) = 0 ;
    
    /** @brief Interface for the calculation of the Thermal Conductivity of heavier species  
     **  @param gasmix GasMixture object 
     **  @param order Desired order of approximation */
    virtual double ThermalCondHeavy ( GasMixture* gasmix, int order ) = 0 ;
    
    /** @brief Interface for the calculation of the Viscosity   
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation */
    virtual double Viscosity ( GasMixture* gasmix, int order ) = 0 ;
    
    /** @brief Interface for the calculation of Electrical conductivity 
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation */
    virtual double ElCond ( GasMixture* gasmix, int order ) = 0 ;
    
    /** @brief Interface for the calculation of Heat Exchange between electrons and heavier species   
     ** @param gasmix GasMixture object */
    virtual double Qeh ( GasMixture* gasmix ) = 0 ; 
    
    /** @brief Interface for the calculation of Diffusion Coefficients
     ** @param gasmix GasMixture object
     ** @param order Desired order of approximation
     ** @param i-th specie 
     ** @param j-th specie */
    virtual double Dij ( GasMixture* gasmix , int order , int i , int j ) = 0 ;
    
    /** @brief Interface for the calculation of Thermal Diffusion Coefficients
     ** @param gasmix GasMixture object
     ** @param order Desired order of approximation
     ** @param i-th specie */
    virtual double DiT ( GasMixture* gasmix , int order , int i ) = 0 ;

    public :

    /** @brief Sets the Chapman-Enskog orders for the transport properties. \n
     ** Orders have to be provided for all transport properties as a vector of integers. \n
     ** For class DevotoTP components orders are: \n
     ** { D, DT, λₑ, λₕ, η, σ }. \n
     ** For class ZhangMurphyTP components orders are: \n
     ** { D, DT, λₑ, λₕ, η, σ, λₑⁿᵉ, λₕⁿᵉ, Dⁿᵉ, DTⁿᵉ }. \n
     ** @see class DevotoTP and class ZhangMurphyTP constructors for further details. */
    virtual void setOrders ( const std::vector<int>& ord ) ;

};

/** @class Appendix
 ** @brief class to incapsulate methods that compute bracket expressions needed \n 
 ** for the Chapman-Enskog approximation to Transport Properties. \n
 ** Implemented methods are those needed for the Devoto and Bonneoi theory,
 ** please override it to implement differences in other theories transport modules.
 ** @see 1) RS Devoto "Transport Properties of Ionized Monoatomic Gases",
 **  The Physics of fluid, 9,6,June(1966) 
 ** @see 2) RS Devoto "Simplified Expressions for the Transport Properties 
 ** of Ionized Monoatomic Gases", \n  The Physics of fluid, 10,10,October(1967) */
class Appendix : public Properties {

    protected : 

    /// @brief Constructor for Appendix class with assigned CiBox object.
    Appendix ( CiBox* cbx  ) : Properties ( cbx ) {}

    /// @brief Constructor for Appendix class with assigned GasMixture object.
    Appendix ( GasMixture* mix ) : Properties(mix) {} 
    
    /** @brief Some Bracket Expression Coefficients for 
     ** thermal conductivity and diffusion Transport properties \n 
     ** see this class reference 1.
     ** @param gasmix object of class GasMixture
     ** @param m Chapman-Enskog order of approximation
     ** @param p Chapman-Enskog order of approximation
     ** @param i-th specie
     ** @param j-th specie */
    virtual double qmpij ( GasMixture* gasmix , int m , int p , int i , int j ) ; 

    /** @brief some simplified Bracket Expression for electron Transport properties \n 
     ** see this class reference 2. 
     ** @param gasmix object of class GasMixture
     ** @param m Chapman-Enskog order of approximation
     ** @param p Chapman-Enskog order of approximation */
    virtual double qsimpmpij ( GasMixture* gasmix , int m , int p ) ;

    /** @brief Some Bracket Expression Coefficients for viscosity 
     ** see this class reference 1.
     ** @param gasmix object of class GasMixture
     ** @param m Chapman-Enskog order of approximation
     ** @param p Chapman-Enskog order of approximation
     ** @param i-th specie
     ** @param j-th specie */
    virtual double qcapmpij ( GasMixture* gasmix, int m, int p, int i, int j ) ;

    /** @brief Function to select the right collision integral collected in Qt
     ** @param l order of momentum
     ** @param s order of Sonine polynomials
     ** @param i-th specie 
     ** @param j-th specie 
     ** @param N_specs Number of species */
    virtual double Qmpil ( GasMixture* gasmix , int l , int s , int i , int j ) ; 

    /** @brief Delta function needed for 1st order approcimations 
     ** @param gasmix object of class GasMixture
     ** @param i-th specie
     ** @param j-th specie */
    double DeltaIJ1 (GasMixture* gasmix , int i , int j ) ;

    /** @brief Delta function needed for 2nd order approcimations 
     ** @param gasmix object of class GasMixture
     ** @param i-th specie
     ** @param j-th specie */
    double DeltaIJ2 (GasMixture* gasmix , int i , int j ) ;

    /** @brief Alpha function needed for 1st order approcimations 
     ** @param mass vector of species masses
     ** @param i-th specie
     ** @param j-th specie */
    double alphaIJ(std::vector<double> mass, int i, int j ) ;

    /** @brief Function for 1st order approximation of thermal conductivity.
     ** @details BEWARE: It doesn't separate electrons and heavy particle contributions. 
     ** Use only in LTE.
     ** @param gasmix object of class GasMixture
     ** @return double 1st order approximation of thermal conductivity */
    double Ktr1(GasMixture* gasmix);

    /** @brief Function for 1st order approximation of viscosity.
     ** @param gasmix object of class GasMixture
     ** @return double 1st order approximation of viscosity */
    double Viscosity1(GasMixture* gasmix);

    /** @brief Function for 1st order approximation of electrical conductivity.
     ** @param gasmix object of class GasMixture
     ** @return double 1st order approximation of electrical conductivity */
    double Sigma1(GasMixture* gasmix);

} ;

#endif
