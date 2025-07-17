#ifndef ZMAPPENDIX_H
#define ZMAPPENDIX_H

#include"Transport.h"

/**
 * @brief class to incapsulate Appendix of bracket integrals expressions needed for the 
 * Transport Properties computation of the 
 * @see Xiao-Ning Zhang, He-Ping Li, Anthony B. Murphy, Wei-Dong Xia, \n 
 * "A Numerical Model of Non-Equilibrium Thermal Plasmas. \n 
 * I. Transport Properties," Phys. Plasmas, vol. 20, p. 033508, 2013.
*/
class ZMAppendix : public Appendix {

    protected : 

    ZMAppendix ( CiBox* cbx  ) : Appendix ( cbx ) {}
    ZMAppendix ( GasMixture* mix ) : Appendix ( mix ) {}
    
    /**
     * @brief Function to select the right collision integral collected in Qt
     * @details Differently from the theory of Devoto and Bonnefoi 
     * the Integrals involved in the computation of bracket expression 
     * are the dimensional ones. \n 
     * see this class reference eq. A1 to A5
     * @param gasmix GasMixture object
     * @param l order of momentum
     * @param s order of Sonine polynomials
     * @param i-th specie 
     * @param j-th specie 
    */
    double Qmpil ( GasMixture* gasmix , int l , int s , int i , int j )  override ; 

    /** @brief Some Bracket Expression Coefficients for viscosity \n
     *  see this class reference eq.14 and article supplementary material
     *  @param gasmix object of class GasMixture
     *  @param m Chapman-Enskog order of approximation
     *  @param p Chapman-Enskog order of approximation
     *  @param i-th specie
     *  @param j-th specie
    */
    double qcapmpij     ( GasMixture* gasmix , int m , int p , int i , int j ) override ;
    
    /** @brief Some Bracket Expression Coefficients for electron scattering in viscosity computation \n
     *  see this class reference eq.15 and article supplementary material
     *  @param gasmix object of class GasMixture
     *  @param m Chapman-Enskog order of approximation
     *  @param p Chapman-Enskog order of approximation
     *  @param i-th specie
     *  @param j-th specie
    */
    double qmpi1V       ( GasMixture* gasmix , int m , int p , int i ) ;
    
    /**
     * @brief Some Bracket Expression Coefficients for electron scattering in 
     * thermal conductivity and diffusion computation \n
     * see this class reference eq.15 and article supplementary material
     * @param gasmix object of class GasMixture
     * @param m Chapman-Enskog order of approximation
     * @param p Chapman-Enskog order of approximation
     * @param i-th specie
     * @param j-th specie
    */
    double qmpi1        ( GasMixture* gasmix , int m , int p , int i ) ;
    
    /** @brief Some Bracket Expression Coefficients for thermal conductivity and diffusion \n
     * @param gasmix object of class GasMixture
     * @param m Chapman-Enskog order of approximation
     * @param p Chapman-Enskog order of approximation
     * @param i-th specie
     * @param j-th specie
    */
    double qmpij        ( GasMixture* gasmix , int m , int p , int i , int j ) override ;  

    /**
     * @brief some simplified Bracket Expression for electron Transport properties \n
     * see this class reference eq.11 and article supplementary material
     * @param gasmix object of class GasMixture
     * @param m Chapman-Enskog order of approximation
     * @param p Chapman-Enskog order of approximation
     * @attention ! This equation preserve the recommendation stated in RS Devoto "Simplified Expressions for the Transport Properties 
     * of Ionized Monoatomic Gases", \n  The Physics of fluid, 10,10,October(1967) ! \n
     * i.e. neglecting all the omega integrals with l>1 in the first term of eq.11
    */
    double qsimpmpij    ( GasMixture* gasmix , int m , int p ) override ;
    
    /**
     * @brief Bracket expression \n
     * see this class reference eq.12 and article supplementary material
     * @param gasmix GasMixture object
     * @param m Chapman-Enskog order of approximation
     * @param p Chapman-Enskog order of approximation
     * @param i-th specie 
     * @param j-th specie 
    */
    double qmpijBar     ( GasMixture* gasmix, int m , int p , int i , int j ) ; 
    
    /**
     * @brief Bracket expression \n
     * see this class reference eq.13 and article supplementary material
     * @param gasmix GasMixture object
     * @param m Chapman-Enskog order of approximation
     * @param p Chapman-Enskog order of approximation
     * @param i-th specie 
     * @param j-th specie 
    */
    double qmpi1Bar     ( GasMixture* gasmix, int m , int p , int i ) ; 

};

// ______________________________ Implementation ______________________________

#endif