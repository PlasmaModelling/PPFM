 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef ZMTRANSPORT_H
#define ZMTRANSPORT_H

/** @file ZhangMurphyTP.h
 * @brief This file contains the ZhangMurphyTP class, 
 * which implements the transport properties calculation based on the theory of 
 * Zhang et. al.
*/

#include"ZMCoefficients.h"

/** @class ZhangMurphyTP
 ** @brief Implementation of the theory formalized by Zhang et. al. in \n
 ** @see 1) Xiao-Ning Zhang, He-Ping Li, Anthony B. Murphy, Wei-Dong Xia, \n 
 ** "A Numerical Model of Non-Equilibrium Thermal Plasmas. \n 
 ** I. Transport Properties," Phys. Plasmas, vol. 20, p. 033508, 2013.
 ** @see 2) Comparison of the transport properties of 
 ** two-temperature argon plasmas calculated using different methods. */
class ZhangMurphyTP : public ZMCoefficients {

    private:
 
    /// @brief Matrix of non equilibrium diffusion coefficients, see eq.22 of this class reference 1.
    std::vector<std::vector<double>> DiffTheta ;

    /// @brief Vector of non equilibrium thermal diffusion coefficients, see eq.36 of this class reference 1.
    std::vector<double> Dtheta ;

    /** @brief Function to compute Electrons Thermal Conductivity as in eq.34 of
     ** this class reference 1.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation
     ** @attention To confront electron thermal conductivity 
     ** with other theories it has to be divided by theta as stated in eq.14 
     ** of this class reference 2 */
    double ThermalCondEl    ( GasMixture* gasmix, int order ) override ;
    
    /** @brief Function to compute heavy particle Thermal Conductivity as in eq.34 of
     ** this class reference 1.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation */
    double ThermalCondHeavy ( GasMixture* gasmix, int order ) override ;
    
    /** @brief Function to compute Viscosity as in eq.A.18 of
     ** this class reference 2.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation */
    double Viscosity        ( GasMixture* gasmix, int order ) override ;

    /** @brief Function to compute Electrical Consuctivity as in eq.42 of
     ** this class reference 1.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation */
    double ElCond           ( GasMixture* gasmix, int order ) override ;

    double Qeh              ( GasMixture* gasmix )            override { return 0. ; } ; 

    /** @brief Function to compute Ordinary Diffusion Coefficients as in eq.20 of
     ** this class reference 1.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation */
    double Dij ( GasMixture* gasmix , int order , int i , int j ) override ;

    /** @brief Function to compute Thermal Diffusion Coefficients as in eq.21 of
     ** this class reference 2.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation */
    double DiT ( GasMixture* gasmix , int order , int i ) override ;

    /** @brief Function to compute n.e. Ordinary Diffusion Coefficients as in eq.22 of
     ** this class reference 1.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation */
    double DijTheta ( GasMixture* gasmix , int order , int i , int j ) ;
    
    /** @brief Function to compute n.e. thermal Diffusion Coefficients as in eq.36 of
     ** this class reference 1.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation */
    double DiTheta ( GasMixture* gasmix , int order , int i ) ;
    
    /** @brief Function to compute n.e. Electrons Thermal Conductivity as in eq.35 of
     ** this class reference 1.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation
     ** @attention To confront electron thermal conductivity 
     ** with other theories it has to be divided by theta as stated in eq.14 
     ** of this class reference 2 */
    double NeThermalCondEl    ( GasMixture* gasmix, int order ) ;
    
    /** @brief Function to compute n.e. Heavy particle Thermal Conductivity as in eq.35 of
     ** this class reference 1.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation */
    double NeThermalCondHeavy ( GasMixture* gasmix, int order ) ;
    
    public:

    /** @brief Constructor with edited CiBox. 
     ** @details Default orders assigned for transport properties.
     ** [D: 3 , DT:3 , Tp[0]:3 , Tp[1]:2 , Tp[2]:1 , Tp[3]:4 , Tp[4]: 3 , Tp[5] :2 , DiffTheta:3 , Dtheta:3 ]  \n
     ** Ordinary Diffusion coefficients,
     ** Thermal Diffusion coefficients,
     ** Electron Thermal Conductivity,
     ** Heavy species Thermal Conductivity,
     ** Viscosity,
     ** Electrical Conductivity,
     ** Electron Thermal Conductivity,
     ** n.e. Electron Thermal Conductivity,
     ** n.e. Heavy species Thermal Conductivity,
     ** n.e. Ordinary Diffusion coefficients,
     ** n.e. Thermal Diffusion coefficients, respectively. */
    ZhangMurphyTP ( CiBox* cbx  ) : ZMCoefficients ( cbx ) {
        orders = {3,3,3,2,1,4,3,2,3,3};
    } 

    /** @brief Constructor with default CiBox created from GasMixture. 
     ** @details Default orders assigned for transport properties.
     ** @see ZhangMurphyTP ( CiBox* cbx ) for details */
    ZhangMurphyTP ( GasMixture* mix ) : ZMCoefficients ( mix ) {
        orders = {3,3,3,2,1,4,3,2,3,3};
    } 

    /** @brief Function that computes and store the Transport Coefficients in 
     ** the base class Transport members
     ** @param gasmix GasMixture object */
    void computeTransport ( GasMixture* gasmix ) override ;

    
};

//______________________________ Implementation ______________________________

#endif