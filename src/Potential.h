#ifndef POTENTIAL_H
#define POTENTIAL_H 

#include "Species.h"

//lista dei parametri dedotta dal tipo implementata come specializzazione
class Potential {
    
    public : 

    virtual double Pot(double r) = 0 ;

};

class ShieldedCoulombPotential : public Potential {

    double z1,z2,lambdaD;
    
    public:
    
    ShieldedCoulombPotential(double Z1, double Z2, double LambdaD) : 
        z1(Z1), z2(Z2), lambdaD(LambdaD) {  }
    
    double Pot (double r) override ; 

};

// Potentiale di Lennard Jones
class LennardJones : public Potential {
    
    protected:
    // [eps] = eV 
    double eps  ;
    // [sig] = Ang
    double sig  ;

    public : 

    LennardJones(double eps, double sig) : eps(eps), sig(sig) {}

    virtual double Pot (double r) override ;

    // dichiarazioni potenziali per hard-code, realizzazioni overload in HCPotentialDB.cpp
    LennardJones(Argon* ar, Argon* arr)  ;
    LennardJones(Helium* t1, Helium* t2) ;
    LennardJones(Argon* t1, Helium* t2)  ;
};

// Potentiale di Capitelli
class Capitelli : public Potential {
    
    protected : 

    Species* t1 ; 
    Species* t2 ;

    // [beta] = # 
    double beta  ; 
    // [De] = eV 
    double De  ; 
    // [re] = Ang
    double Re  ;

    public :
    
    Capitelli(Species* sp1, Species* sp2, double beta, double De, double Re) : 
        beta(beta), De(De), Re(Re) {}

    virtual double Pot (double r) override ;

    // dichiarazioni potenziali per hard-code, realizzazioni in HCPotentialDB.cpp
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

// Potenziale Hilburt-Hirschfelder ridotto
class HulburtHirschfelder : public Potential {
    
    protected :
    // [De] = eV 
    double De  ;
    // [re] = Ang
    double Re  ;
    
    double a ,b ,c  ;

    public :

    HulburtHirschfelder(double De, double Re, double a , double b, double c) : 
        De(De), Re(Re), a(a), b(b), c(c) {}

    virtual double Pot (double r) override ;

};

// Potenziale Hilburt-Hirschfelder non ridotto
class HulburtHirschfelderUnreduced : public Potential {
    
    protected :
    
    double eps ; //[eV]
    double  re ; //[Ang] or [m]
    double we, weXe, Be, Alphae ; //[re]^(-1)

    public :

    HulburtHirschfelderUnreduced(double eps, double re, double we , double weXe, double Be, double Alphae) : 
        eps(eps), re(re), we(we), weXe(weXe), Be(Be), Alphae(Alphae) {}

    virtual double Pot (double r) override ;

};


class Morse : public Potential {} ;

class PowerPot : public Morse {

    double v0 , n ; 

    public :
    PowerPot( double v0, double n ) : v0(v0), n(n){}

    double Pot ( double r ) override ; 

} ; 
class Morse2Param : public Morse
{
    protected:
    // [V0] = eV
    double V0 ;
    // [B] = Ang^-1
    double B ;

    public:
    
    Morse2Param(double V0, double B) : V0(V0) , B(B) {}
    
    virtual double Pot ( double r ) override ; 

};

class Morse3Param : public Morse {
    
    protected : 
    // [De] = eV 
    double De  ;
    // [Beta] = Ang^-1 
    double Beta  ;
    // [Re] = Ang
    double Re  ;

    public :

    Morse3Param(double De, double Beta, double Re) : De(De) , Re(Re) , Beta(Beta) {}

    virtual double Pot (double r) override ;
};
class Morse5Param : public Morse {
    
    double re, D, b0, gamma, lambda ;

    public :
     
    Morse5Param( double RE, double DD , double B0, double GAMMA, double LAMBDA ) : 
        re(RE) , D(DD) , b0(B0), gamma(GAMMA), lambda(LAMBDA) {}
    
    virtual double Pot(double r) override ;

};

class Polarization : public Potential {

    std::tuple < Species*, double > inducedDipole ; 
    ChargedSpecies* ion ; 

    public : 

    Polarization ( Species* sp1, Species* sp2, double Polarizability ) ;

    virtual double Pot(double r) override ;

} ; 

/** @brief ad-hoc potential to validate Zhang study, 
 * @see "The repulsive wall of the Ar-Ar interatomic potential reexamined"
    Ronald A. Aziz and M. J. Slaman
    The Journal of Chemical Physics 92, 1030 (1990); doi: 10.1063/1.458165 */
class HFDTCS2_ArAr : public Potential {
public:

    HFDTCS2_ArAr () ;
    double Pot ( double r ) override ;
    Potential* clone () ;

} ;

#endif
