 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef DEVOTO_H
#define DEVOTO_H

#include "Transport.h"

/** @class DevotoTP
 ** @brief Implementation of the theory formalized by R.S Devoto in \n
 ** @see 1) R.S. Devoto "Transport Properties of Ionized Monoatomic Gases", \n 
 ** The Physics of fluid, 9,6,June(1966) \n
 ** and \n
 ** @see 2) R.S. Devoto "Simplified Expressions for the Transport Properties 
 ** of Ionized Monoatomic Gases", \n  The Physics of fluid, 10,10,October(1967) */
class DevotoTP : public Appendix {

    protected : 
    
    /**  @brief Function to compute Electrons Thermal Conductivity as in eq.20 of
     ** this class reference 2.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation */
    double ThermalCondEl ( GasMixture* gasmix, int order ) override ;

    /**  @brief Function to compute Electrons Thermal Conductivity as in eq.14 of
     ** this class reference 1.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation  */
    double ThermalCondHeavy ( GasMixture* gasmix, int order ) override ;

    /** @brief Function to compute Viscosity as in eq.21 of
     ** this class reference 1.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation  */
    double Viscosity ( GasMixture* gasmix, int order ) override ;

    /** @brief Function to compute Electrical Conductivity as in eq.16 of
     ** * this class reference 2.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation 
    */
    double ElCond ( GasMixture* gasmix, int order ) override ;
    
    /** @brief Function to compute Heat Exchange of electrons-heavy collisions
     ** @param gasmix GasMixture object  */
    double Qeh ( GasMixture* gasmix ) override ; 
    
    /** @brief Function to compute Diffusion Coefficients as in eq.8 of
     ** this class reference 1.
     ** @param gasmix GasMixture object
     ** @param order Desired order of approximation
     ** @param i-th specie 
     ** @param j-th specie */
    double Dij ( GasMixture* gasmix , int order , int i , int j ) override ; 

    /** @brief Function to compute Thermal Diffusion Coefficients as in eq.9 of
     ** this class reference 1.
     ** @param gasmix GasMixture object
     ** @param order Desired order of approximation
     ** @param i-th specie  */
    double DiT ( GasMixture* gasmix , int order , int i ) override ;

    public : 

    /** @brief Constructor with edited CiBox. 
     ** @details Default orders assigned for transport properties are  
     ** [D: 3 , DT:4 , Tp[0]:3 , Tp[1]:2 , Tp[2]:1 , Tp[3]:4 , Tp[4]: - ] \n
     ** Ordinary Diffusion coefficients,
     ** Thermal Diffusion coefficients,
     ** Electron Thermal Conductivity,
     ** Heavy species Thermal Conductivity,
     ** Viscosity,
     ** Electrical Conductivity, respectively. */
    DevotoTP (  CiBox* cbx  ) : Appendix ( cbx ) {
        orders = {3,3,3,2,2,3} ;
    } 

    /** @brief Constructor with default CiBox created from GasMixture. 
     ** @details Default orders assigned for transport properties.
     ** @see DevotoTP (  CiBox* cbx  ) constructor for details. */
    DevotoTP ( GasMixture* mix ) : Appendix ( mix ) {
        orders = {3,3,3,2,2,3} ;
    } 
    
    /** @brief Function that computes and store the Transport Coefficients in 
     ** the base class Transport members
     ** @param gasmix GasMixture object */
    void computeTransport ( GasMixture* gasmix ) override ;

    void lightComputeTransport ( GasMixture* gasmix ) ;

} ;

/**
 * @brief Class that extends DevotoTP adding the reactive contribution to heavy species thermal conductivity
 * in case of LTE assumption following Brokaw et al. and Colombo et. al. 
 * @see 3) "Thermodynamic and transport properties in non-equilibrium argon, oxygen and nitrogen thermal plasmas" */
class DevotoLteLambdaR : public DevotoTP {

    protected : 

    double Dbinij(int i,int j, GasMixture* gasmix) ;  

    public :

    /** @brief Constructor with edited CiBox and Thermodynamics object. 
     ** @details Default orders assigned for transport properties are  
     ** [D: 3 , DT:4 , Tp[0]:3 , Tp[1]:2 , Tp[2]:1 , Tp[3]:4 , Tp[4]: - ] \n
     ** Ordinary Diffusion coefficients,
     ** Thermal Diffusion coefficients,
     ** Electron Thermal Conductivity,
     ** Heavy species Thermal Conductivity,
     ** Viscosity,
     ** Electrical Conductivity, respectively. */
    DevotoLteLambdaR (  CiBox* cbx ) 
    : DevotoTP ( cbx ) {}

    /**
     * @brief Constructor with default CiBox created from GasMixture and Thermodynamics object.
     * @details Default orders assigned for transport properties.
     * @see DevotoLteLambdaR (  CiBox* cbx , Thermodynamics* thermo  ) constructor for details. */
    DevotoLteLambdaR ( GasMixture* mix ) 
    : DevotoTP ( mix ) {}

    /** @brief Computes as in Devoto and add the reactive contributions 
     * to heavy particle thermal conductivity following Brokaw et al. and 
     * especially with the formula described in Colombo et. al. 
     * "Thermodynamic and transport properties in non-equilibrium argon,
     * oxygen and nitrogen thermal plasmas"
     * @param gasmix GasMixture object */
    void computeTransport ( GasMixture* gasmix ) override ;

} ;

#endif