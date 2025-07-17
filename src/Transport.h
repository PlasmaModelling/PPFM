#ifndef TRANSPORT_H
#define TRANSPORT_H 

#include <vector>
#include <string>
#include <numbers>

// forward declarations
class GasMixture ; 
class CiBox ; 

#define kB 1.380650524e-8 
#define PI std::numbers::pi  
#define q_e 1.60217646e-19
/**
 * @class Transport
 * @brief Base class common to all Transport Properties theories
 * @details This class contains members to collect transport properties
 * and Diffusion Coefficients in SI units. \n
 * It also serves as a bridge to Collision Integrals computation module \n
 * through a pointer to a class CiBox, that's used to compute and store \n
 * Collision Integrals in the matrix Qt through the method QtCalc().
*/
class Transport {

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
     *  @details the rows span from 0 to 15, the number of (m,p) approximations in the Chapman-Enskog method, \n
     *  the columns span from 0 to Number of Species in the Mixture . \n
     *  default units are [ micron^2 ]
    */
    std::vector<std::vector<double>> Qt ; 

    Transport ( CiBox* cbx ) : Ci(cbx)  {} 
    Transport ( GasMixture* mix ) ;

    /** 
     * @brief Compute and store Collision Integrals values 
     * in Qt via the pointer to CiBox object. \n 
     * @attention ! Ensure a call to this method at every 
     * computeTransport implementation ! \n 
     * @details the default method compute and stores Collision Integrals 
     * in the adimensional form
     * @see <a href="../../articles/TransportFormulas.pdf"> TransportFormulas.pdf </a>
     * eq. 17
    */
    void QtCalc ( GasMixture* gasmix ) ;

    public : 

    /** @brief Common Interface to compute Transport properties 
     * with different theories, \n  
     * please do refer to implemented modules and place 
     * computed Transport Coefficients in SI units when implementing.
    */
    virtual void computeTransport ( GasMixture* gasmix ) = 0 ; 
    
} ;

/**
 * @class Properties
 * @brief Class to incapsulate computation methods
 * @details methods will be unaccesible to user, \n
 * the user can access to Transport Properties via the methods of the base class Transport \n 
 * after calling for the computation of them all via the method computeTransport(GasMixture*). \n
 * Methods incapsulated here will then serve to implement derived classes of Transport Properties computational theories.
*/
class Properties : public Transport {

    protected :

    Properties ( CiBox* cbx  ) : Transport ( cbx ) {}
    Properties ( GasMixture* mix ) : Transport ( mix ) {} ;

    /** @brief Interface for the calculation of the Thermal Conductivity of electrons
     *  @param gasmix GasMixture object 
     *  @param order Desired order of approximation
    */
    virtual double ThermalCondEl ( GasMixture* gasmix, int order ) = 0 ;
    
    /** @brief Interface for the calculation of the Thermal Conductivity of heavier species  
     *  @param gasmix GasMixture object 
     *  @param order Desired order of approximation 
    */
    virtual double ThermalCondHeavy ( GasMixture* gasmix, int order ) = 0 ;
    
    /** @brief Interface for the calculation of the Viscosity   
     * @param gasmix GasMixture object 
     * @param order Desired order of approximation 
    */
    virtual double Viscosity ( GasMixture* gasmix, int order ) = 0 ;
    
    /** @brief Interface for the calculation of Electrical conductivity 
     * @param gasmix GasMixture object 
     * @param order Desired order of approximation 
    */
    virtual double ElCond ( GasMixture* gasmix, int order ) = 0 ;
    
    /** @brief Interface for the calculation of Heat Exchange between electrons and heavier species   
     * @param gasmix GasMixture object 
    */
    virtual double Qeh ( GasMixture* gasmix ) = 0 ; 
    /**
     * @brief Interface for the calculation of Diffusion Coefficients
     * @param gasmix GasMixture object
     * @param order Desired order of approximation
     * @param i-th specie 
     * @param j-th specie
    */
    virtual double Dij ( GasMixture* gasmix , int order , int i , int j ) = 0 ;
    /**
     * @brief Interface for the calculation of Thermal Diffusion Coefficients
     * @param gasmix GasMixture object
     * @param order Desired order of approximation
     * @param i-th specie 
    */
    virtual double DiT ( GasMixture* gasmix , int order , int i ) = 0 ;

};

/**
 * @class Appendix
 * @brief class to incapsulate methods that compute bracket expressions needed \n 
 * for the Chapman-Enskog approximation to Transport Properties. \n
 * Implemented methods are those needed for the Devoto and Bonneoi theory,
 * please override it to implement differences in other theories transport modules.
 * @see 1) RS Devoto "Transport Properties of Ionized Monoatomic Gases",
 *  The Physics of fluid, 9,6,June(1966) 
 * @see 2) RS Devoto "Simplified Expressions for the Transport Properties 
 * of Ionized Monoatomic Gases", \n  The Physics of fluid, 10,10,October(1967)
*/
class Appendix : public Properties {

    protected : 

    Appendix ( CiBox* cbx  ) : Properties ( cbx ) {}
    Appendix ( GasMixture* mix ) : Properties(mix) {} 
    
    /** @brief Some Bracket Expression Coefficients for 
     * thermal conductivity and diffusion Transport properties \n 
     * see this class reference 1.
     * @param gasmix object of class GasMixture
     * @param m Chapman-Enskog order of approximation
     * @param p Chapman-Enskog order of approximation
     * @param i-th specie
     * @param j-th specie
    */
    virtual double qmpij ( GasMixture* gasmix , int m , int p , int i , int j ) ; 

    /**
     * @brief some simplified Bracket Expression for electron Transport properties \n 
     * see this class reference 2. 
     * @param gasmix object of class GasMixture
     * @param m Chapman-Enskog order of approximation
     * @param p Chapman-Enskog order of approximation
    */
    virtual double qsimpmpij ( GasMixture* gasmix , int m , int p ) ;

    /** @brief Some Bracket Expression Coefficients for viscosity 
     *  see this class reference 1.
     *  @param gasmix object of class GasMixture
     *  @param m Chapman-Enskog order of approximation
     *  @param p Chapman-Enskog order of approximation
     *  @param i-th specie
     *  @param j-th specie
    */
    virtual double qcapmpij ( GasMixture* gasmix, int m, int p, int i, int j ) ;

    /**
     * @brief Function to select the right collision integral collected in Qt
     * @param l order of momentum
     * @param s order of Sonine polynomials
     * @param i-th specie 
     * @param j-th specie 
     * @param N_specs Number of species
    */
    virtual double Qmpil ( GasMixture* gasmix , int l , int s , int i , int j ) ; 

} ;

// ______________________________ Implementation ______________________________

#endif
