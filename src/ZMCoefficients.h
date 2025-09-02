 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef ZMCOEFFICIENTS_H
#define ZMCOEFFICIENTS_H

#include"ZMAppendix.h"
/**
 * @brief Coefficients of the solved Boltzmann equation as described in 
 * @see 1) Liheping Supplementary material for A numerical model of 
 * non-equilibrium thermal plasmas. \n I. Transport properties \n
 * and \n 
 * @see 2) Comparison of the transport properties of 
 * two-temperature argon plasmas calculated using different methods.
 * @see 3) Xiao-Ning Zhang, He-Ping Li, Anthony B. Murphy, Wei-Dong Xia, \n 
 * "A Numerical Model of Non-Equilibrium Thermal Plasmas. \n 
 * I. Transport Properties," Phys. Plasmas, vol. 20, p. 033508, 2013.
*/
class ZMCoefficients : public ZMAppendix {

    protected:

    std::vector<std::vector<double>> Qtilde ; 
    std::vector<std::vector<double>> Q1 ; 
    
    std::vector<double> n ; 
    std::vector<double> mass ; 
    int N ;
    double me,ne,T,theta,Te,rho,ntot ;


    protected : 

    ZMCoefficients ( CiBox* cbx  ) : ZMAppendix ( cbx ) {}
    ZMCoefficients ( GasMixture* mix  ) : ZMAppendix ( mix ) {}

    void init(GasMixture* gasmix ) ;

    /// @brief coefficient for viscosity see eq. A.20 of this class reference 2)
    double b10 ( GasMixture* gasmix ) ; 
    /// @brief coefficient for viscosity see eq. A.20 of this class reference 2)
    double b11 ( GasMixture* gasmix ) ; 

    /** @brief coefficient for thermal conductivity and diffusion, see eq. S-47
    * of this class reference 1 */
    double c1psk ( GasMixture* gasmix , int epsilon , int s , int k , int p ) ;
    
    /** @brief coefficient for thermal conductivity and diffusion, see eq. S-48
    * of this class reference 1 */
    double ci0   ( GasMixture* gasmix , int epsilon , int s , int k , int i ) ;

    /** @brief coefficient for thermal conductivity and diffusion, see eq. S-55
    * of this class reference 1 */
    double c11   ( GasMixture* gasmix , int epsilon , int s , int k ) ; 

    /** @brief coefficient for thermal conductivity and diffusion, see eq. S-56
    * of this class reference 1 */
    double ci1   ( GasMixture* gasmix , int epsilon , int s , int k , int i ) ;

    /** @brief coefficient for thermal conductivity and diffusion, see eq. S-52
    * of this class reference 1, 
    * also f1p, see TABLE III */
    double aip   ( GasMixture* gasmix , int epsilon , int i , int p ) ;

    /** @brief coefficient for thermal conductivity and diffusion, see eq. S-53
    * of this class reference 1 */
    double ai0   ( GasMixture* gasmix , int epsilon , int i ) ;

    /** @brief coefficient for thermal conductivity and diffusion, see eq. S-60
    * of this class reference 1 
    * also, eq. S-57 */
    double a11   ( GasMixture* gasmix , int epsilon ) ;
    
    /** @brief coefficient for thermal conductivity and diffusion, see eq. S-61
    * of this class reference 1 */
    double ai1   ( GasMixture* gasmix , int epsilon , int i ) ;

    double e1psk ( GasMixture* gasmix , int epsilon , int pp , int s , int k ) ;
    double ei0   ( GasMixture* gasmix , int epsilon , int ii , int s , int k ) ;    
    double e11sk ( GasMixture* gasmix , int epsilon , int s , int k ) ; 
    double ei1   ( GasMixture* gasmix , int epsilon , int ii , int s , int k ) ;    
    
    double f10   ( GasMixture* gasmix , int epsilon ) ; 
    double fi0   ( GasMixture* gasmix , int epsilon , int ii ) ; 
    double f11   ( GasMixture* gasmix , int epsilon ) ; 
    double fi1   ( GasMixture* gasmix , int epsilon , int ii ) ; 
    double wi    ( GasMixture* gasmix , int i ) ;

    /** @brief diffusion thermal conductivity, see eq. 32
    * of this class reference 3 */
    double lambdaijD    ( GasMixture* gasmix, int epsilon, int i , int j ) ;
    
    /** @brief first term of thermal conductivity, see eq. 30
    * of this class reference 3 */
    double lambdaiPrime ( GasMixture* gasmix, int epsilon, int i ) ;

    /** @brief first term of n.e. thermal conductivity, see eq. 31
    * of this class reference 3 */
    double lambdaiPrimeTheta ( GasMixture* gasmix, int epsilon, int i ) ;

    /** @brief first term of thermal conductivity, see eq. 23
    * of this class reference 3 */
    double DiThetaStar ( GasMixture* gasmix, int epsilon, int i ) ;

};

//______________________________ Implementation ______________________________

#endif