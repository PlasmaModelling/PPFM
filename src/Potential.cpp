
#ifndef HCPOTENTIALS_H
#define HCPOTENTIALS_H

#include <stdexcept>
#include <numbers>
#include "Potential.h"

// _____________________________ Potentials expressions _____________________________


Polarization::Polarization ( Species* sp1, Species* sp2, double Polarizability ) {
    
    auto p1 = dynamic_cast<ChargedSpecies*>(sp1) ; 
    auto p2 = dynamic_cast<ChargedSpecies*>(sp2) ; 

    if ( p1 != nullptr && p2 == nullptr ) {
        
        ion = p1 ; 
        inducedDipole = std::make_tuple(sp2,Polarizability) ;

    } else if ( p1 == nullptr && p2 != nullptr ) {

        ion = p2 ; 
        inducedDipole = std::make_tuple(sp1,Polarizability) ;

    } else {

        throw std::invalid_argument( "One species must be charged (ChargedSpecies*), "
            "and the other must be neutral (Species*)." ) ;

    }
    
}

// ShieldedCoulombPotential::Pot
double ShieldedCoulombPotential::Pot(double r) {
    
    return ((z1 * z2 * qe) / (4 * std::numbers::pi * eps0 * r * 1e-10)) * exp(-(r / lambdaD));

}

// LennardJones::Pot
double LennardJones::Pot(double r) {
    double sr = sig / r;
    return 4 * eps * ((pow(sr, 12)) - (pow(sr, 6)));
}

// Capitelli::Pot
double Capitelli::Pot(double r) {
    double m = (t1->getCharge() == 0 && t2->getCharge() == 0) ? 4 : 6;
    double x = r / Re;
    double n = beta + 4 * std::pow(x, 2);
    return De * ((m / (n - m)) * std::pow(1 / x, n) - (n / (n - m)) * std::pow(1 / x, m));
}

// HulburtHirschfelder::Pot
double HulburtHirschfelder::Pot(double r) {
    
    double rsd = r / Re;
    return De * (std::exp(-2 * a * (rsd - 1)) - 2 * std::exp(-a * (rsd - 1)) +
                 b * std::pow(rsd - 1, 3) * (1 + c * (rsd - 1)) * std::exp(-2 * a * (rsd - 1)));

}

double HulburtHirschfelderUnreduced::Pot(double r) {
    
    double a0,a1,a2 ;
    a0 = pow(we,2.)/(4.*Be) ; 
    a1 = -1.-((Alphae*we)/(6.*Be)) ; 
    a2 = (5./4.)*pow(a1,2.)-2.*weXe/(3.*Be) ;
   
    double b,c ; 
    c = 1+a1*sqrt(eps*a0) ; 
    b = 2.-((7./12.)-a2*eps/a0)/c ;
   
    double x ;  
    x = we/(2.*sqrt(Be*eps))*((r-re)/re) ; 

    return eps * ( pow(1.-std::exp(-x),2.) + c*pow(x,3.)*(1+b*x)*std::exp(-2*x) -1 ) ;

}

double PowerPot::Pot ( double r ) {
    
    return v0 * pow(r,-n) ; 

}

// Morse2Param::Pot
double Morse2Param::Pot ( double r ) {

    return V0 * std::exp ( -B * r ) ;

}

// Morse3Param::Pot
double Morse3Param::Pot ( double r ) {
    
    return De * ( std::exp ( -2. * Beta * (r - Re) ) - 2. * std::exp ( -Beta * (r - Re) ));

}

// Morse5Param::Pot
double Morse5Param::Pot ( double r ) {

    double b = b0 * ( 1 + gamma * (r - re) + lambda * pow( (r - re), 2 ) );
    return D * ( pow ( 1 + exp ( -b * (r - re)), 2 ) - 1.0 );

}

double Polarization::Pot ( double r ) {
    
    double rm = r ; // Ang 

    Species* neutral = std::get<0> ( inducedDipole ) ;
    double alpha = std::get<1> ( inducedDipole ) ; // Ang^3

    double q = ion->getCharge() * qe ; // C

    double V = - ( 1.e+10 * alpha * std::pow ( q, 2. ) ) / ( 8. * std::numbers::pi * eps0 * std::pow ( r, 4. ) ) ; // J

    V *= 1/ qe ;

    return V ; // eV

}

// Argon-Argon 
LennardJones::LennardJones(Argon* ar, Argon* arr) { eps = 0.010333 ; sig = 3.405 ; }
// Helium-Helium
LennardJones::LennardJones(Helium* t1, Helium* t2) { eps = 9.39289e-04 ; sig = 2.64 ; }
// Argon-Helium
LennardJones::LennardJones(Argon* t1, Helium* t2) { eps = std::sqrt(9.39289e-04 * 0.010333); sig = (2.64 + 3.405) / 2; }

// Argon-Argon
Capitelli::Capitelli(Argon* t1, Argon* t2) { beta = 8.12; De = 0.011603; Re = 3.794; }
// Argon-ArgonI
Capitelli::Capitelli(Argon* t1, ArgonI* t2) { beta = 7.6; De = 0.107609; Re = 3.213; }
// Argon-ArgonII
Capitelli::Capitelli(Argon* t1, ArgonII* t2) { beta = 7.48; De = 0.753338; Re = 2.624; }
// Argon-Nitrogen
Capitelli::Capitelli(Argon* t1, Nitrogen* t2) { beta = 6.94; De = 0.008505; Re = 3.695; }
// Argon-NitrogenI
Capitelli::Capitelli(Argon* t1, NitrogenI* t2) { beta = 7.37; De = 0.12184; Re = 3.054; }
// Argon-MolecularNitrogen
Capitelli::Capitelli(Argon* t1, MolecularNitrogen* t2) { beta = 9.1; De = 0.011955; Re = 3.906; }
// ArgonI-Nitrogen
Capitelli::Capitelli(ArgonI* t1, Nitrogen* t2) { beta = 6.82; De = 0.082229; Re = 3.121; }
// ArgonI-MolecularNitrogen
Capitelli::Capitelli(ArgonI* t1, MolecularNitrogen* t2) { beta = 7.59; De = 0.112256; Re = 3.229; }
// ArgonII-Nitrogen
Capitelli::Capitelli(ArgonII* t1, Nitrogen* t2) { beta = 6.79; De = 0.581253; Re = 2.537; }
// ArgonII-MolecularNitrogen
Capitelli::Capitelli(ArgonII* t1, MolecularNitrogen* t2) { beta = 7.47; De = 0.784957; Re = 2.639; }



 HFDTCS2_ArAr::HFDTCS2_ArAr() : Potential() {}

double HFDTCS2_ArAr::Pot(double r) {

    const double epsilon = 143.25 * KB / qe ; // Parametro di scala energetica (in eV)
    const double rm = 3.761;         
    const double A = 1.13211845e+5;  
    const double alpha = 9.00053441; 
    const double beta = -2.60270226; 
    const double c6 = 1.09971113;    
    const double c8 = 0.54511632;    
    const double c10 = 0.39278653;   

    // Parametri per la spline esponenziale cubica (regione piece-wise)
    const double r1 = 2.73;             // (in Å)
    const double r2 = 2.98;             // (in Å)
    const double a1 = 3.1790298;
    const double a2 = -18.991157;
    const double a3 = -80.983465;
    const double a4 = 1044.13619;
    const double B = 5780.0;            
    const double b = 3.6182;            

    // Variabile ridotta
    double x = r / rm;
    double x1 = r1 / rm;
    double x2 = r2 / rm;

    double Fx ;
    if ( x <= 1.04 )
        Fx = exp( -pow ( (1.04/x)-1 , 2. )) ; 
    else 
        Fx = 1 ; 
    
    double VStar ; 

    if (r < r1) {

        VStar =  B/epsilon * exp ( - b * rm * x ) ;

    } else if (r >= r1 && r <= r2) {
        
        VStar = exp ( a1 + ( x - x1 ) * ( a2 + ( x - x2 ) * ( a3 + ( x - x1 ) * a4 )));
    
    } else {
    
        VStar = A * exp ( - alpha * x + beta * x * x ) - 
            ( Fx * ( ( c6 / pow(x,6.) ) + ( c8 / pow(x,8.) ) + ( c10 / pow(x,10.) ) )) ;
    
    }

    return epsilon * VStar ;
}

Potential* HFDTCS2_ArAr::clone () { return new HFDTCS2_ArAr() ; }

#endif
