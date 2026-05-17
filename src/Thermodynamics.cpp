 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include<numbers>
#include"Thermodynamics.h"
#include"PfBox.h"
#include"GasMixture.h"
#include"PartitionFunction.h"

void Thermodynamics::computeThermodynamics ( GasMixture& gasmix ) {
    
    int N = gasmix.getN() ;
    std::vector<double> masses = gasmix.masses() ;
    std::vector<double> e_ion ;
    for (int i = 0; i < N; i++)
        e_ion.push_back( gasmix(i)->formationEnergy());
    
    double P = gasmix.getPressure() ;
    double T = gasmix.getTemperature() ;
    double theta = gasmix.theta->get() ;
    double Te = T * theta ;
    
    std::vector<double> n0 (N,1.);
    std::vector<double> nf (N,1.);
    std::vector<double> nb (N,1.);
    std::vector<double> Q (N,1.) ;
    std::vector<double> Qf (N,1.) ;
    std::vector<double> Qb (N,1.) ;
    double dRdT = 1. ;
    double dedT = 1. ;
    double dedT_old = 1. ;
    std::vector<double> dlogQdT_centered (N,1.) ;
    std::vector<double> dlogQdT_forward (N,1.) ;
    std::vector<double> dlogQdT_backward (N,1.) ;
    std::vector<double> dlogQdT_centered_old (N,1.) ;
    std::vector<double> dlogQdT_forward_old (N,1.) ;
    std::vector<double> dlogQdT_backward_old (N,1.) ;
    double dRdP = 1. ; 
    n0 = gasmix.getCompositionObj()->compositions();

    double dT = 50., Tf, Tef, Tb, Teb ;
    double dP = 500.; 

    double Tol = 1.e-3 ; 
    double alpha = 0.05 ; 

    double efTrial = 0., ebTrial = 0.;
    
    std::vector<double> hfTrial (N,1.);
    std::vector<double> hbTrial (N,1.);
    dhdT.resize(N,1.) ;
    std::vector<double> dhdT_old (N,1.) ;

    // Looping on dT : dRdT dlogQdT_centered, forward and backward finite differences.
    double err1 = 1., err2 = 1., err3 = 1., err4 = 1., err5 = 1., err6 = 1. ;
    double dT_old = dT, dRdT_old = 1.;
    while ( err1 > Tol || err2 > Tol || err3 > Tol || err4 > Tol || err5 > Tol || err6 > Tol )
    {    

        dT_old = dT ;
        
        dT = dT_old - alpha * dT ;

        dRdT_old = dRdT ;
        for (size_t i = 0; i < N; i++) {

            dlogQdT_centered_old[i] = dlogQdT_centered[i] ;
            dlogQdT_forward_old[i] = dlogQdT_forward[i] ;
            dlogQdT_backward_old[i] = dlogQdT_backward[i] ;
        
        }

        for (int i = 0; i < N; i++)
            Q[i] = gasmix.getCompositionObj()->getPfBox()->operator()(i) ;

        gasmix.setT(T+dT);
        nf= gasmix.getCompositionObj()->compositions();
        for (int i = 0; i < N; i++)
            Qf[i] = gasmix.getCompositionObj()->getPfBox()->operator()(i) ;

        gasmix.setT(T);
        gasmix.getCompositionObj()->setn0(n0) ;
        
        gasmix.setT(T-dT);
        nb = gasmix.getCompositionObj()->compositions();
        for (int i = 0; i < N; i++)
            Qb[i] = gasmix.getCompositionObj()->getPfBox()->operator()(i) ;

        gasmix.setT(T);
        gasmix.getCompositionObj()->setn0(n0) ;

        Tf = T+dT, Tb = T-dT ;
        Tef = Te + dT*theta, Teb = Te - dT*theta ;

        // R e dRdT
        double Rf = 0., Rb = 0.;
        Rf = masses[N-1] * nf[N-1] * Tef ;
        Rb = masses[N-1] * nb[N-1] * Teb ;    
        for (int i = 0; i < N-1; i++){
                Rf += masses[i]*nf[i]*Tf;
                Rb += masses[i]*nb[i]*Tb;
        }
        Rf = P/Rf ;
        Rb = P/Rb ;

        dRdT = (Rf-Rb)/(2.*dT) ;
        for (int i = 0; i < N; i++) {

            dlogQdT_centered[i] = ( std::log(Qf[i]) - std::log(Qb[i]) ) / ( 2.*dT ) ;
            dlogQdT_forward[i]  = ( std::log(Qf[i])-std::log(Q[i]) )/(dT) ;
            dlogQdT_backward[i] = ( std::log(Q[i])-std::log(Qb[i]) )/(dT) ;

        }

        err1 = std::abs( (dRdT - dRdT_old) / dRdT ) ;
        err2 = std::abs( (Norm(dlogQdT_centered) - Norm(dlogQdT_centered_old)) / Norm(dlogQdT_centered_old) ) ;
        err3 = std::abs( (Norm(dlogQdT_forward) - Norm(dlogQdT_forward_old)) / Norm(dlogQdT_forward_old) ) ;
        err4 = std::abs( (Norm(dlogQdT_backward) - Norm(dlogQdT_backward_old)) / Norm(dlogQdT_backward_old) ) ;

        efTrial = 0.; 
        ebTrial = 0.;
        for (int i = 0; i < N; i++){
            
            bool isElectron = (i == N-1) ;
            double Tif = isElectron ? Tef : Tf ;
            double Tib = isElectron ? Teb : Tb ;

            double Ei = e_ion[i] ;

            efTrial += KB * nf[i]* ( (1.5 *  Tif) + Ei/KB + (Tif * dlogQdT_forward[i]) ) ;
            ebTrial += KB * nb[i]* ( (1.5 *  Tib) + Ei/KB + (Tib * dlogQdT_backward[i]) ) ;

        }
        double rhof = 0., rhob = 0.;
        for (int i = 0; i < N; i++){
            rhof += nf[i]*masses[i] ;
            rhob += nb[i]*masses[i] ;
        }

        dedT_old = dedT ;
        dedT = (efTrial/rhof-ebTrial/rhob)/(2.*dT*theta) ; 

        err5 = std::abs( (dedT - dedT_old) / dedT ) ;

        dhdT_old = dhdT ;
        for (int i = 0; i < N; i++){
            bool isElectron = (i == N-1) ;
            double Tif = isElectron ? Tef : Tf ;
            double Tib = isElectron ? Teb : Tb ;
            double Ei = e_ion[i] ;

            hfTrial[i] = (2.5 * KB * Tif) + Ei + (KB * Tif * dlogQdT_forward[i]) ;
            hbTrial[i] = (2.5 * KB * Tib) + Ei + (KB * Tib * dlogQdT_backward[i]) ;
            dhdT[i] = (hfTrial[i]-hbTrial[i])/(2.*dT) ;
        }

        err6 = std::abs( (Norm(dhdT) - Norm(dhdT_old)) / Norm(dhdT_old) ) ;

    }

    // loop on dP : dRdP finite difference
    double dP_old = dP, dRdP_old = 1. ;
    std::vector<double> nPf (N,1.) ;
    std::vector<double> nPb (N,1.) ;
    err1 = 1. ;
    while ( err1 > Tol )
    {
        dP_old = dP ;
        dP = dP_old - alpha * dP ;

        dRdP_old = dRdP ;

        gasmix.setP(P+dP);
        nPf = gasmix.getCompositionObj()->compositions();
        gasmix.setP(P);
        gasmix.getCompositionObj()->setn0(n0) ;
        gasmix.setP(P-dP);
        nPb = gasmix.getCompositionObj()->compositions();
        gasmix.setP(P);
        gasmix.getCompositionObj()->setn0(n0) ;

        double rhopf = 0., rhopb = 0.;
        for (int i = 0; i < N; i++){
            rhopf += nPf[i]*masses[i] ;
            rhopb += nPb[i]*masses[i] ;
        }

        double Rf = 0., Rb = 0.;
        for (int i = 0; i < N-1; i++){
            Rf += masses[i]*nPf[i]*T;
            Rb += masses[i]*nPb[i]*T;
        }
        Rf += masses[N-1]*nPf[N-1]*Te;
        Rb += masses[N-1]*nPb[N-1]*Te;
        Rf = (P+dP)/Rf ;
        Rb = (P-dP)/Rb ;

        dRdP = (Rf-Rb)/(2.*dP) ;

        err1 = std::abs( (dRdP - dRdP_old) / dRdP ) ;
    }   

    // Iteration finished, now compute thermodynamic properties

    // std::cout << "dT = " << dT << " K, dP = " << dP << " Pa, dRdT = " << dRdT << " Pa/K, dRdP = " << dRdP << " 1/Pa" << std::endl ;

    double rho = 0.;
    for (int i = 0; i < N; i++)
        rho += n0[i]*masses[i] ;
    
    double R = 0.;
    for (int i = 0; i < N-1; i++)
        R += masses[i]*n0[i]*T;
    R += masses[N-1]*n0[N-1]*Te;
    R = P/R ;

    henthalpies.resize(N,0.);
    double hh = 0., he = 0., ee = 0., eh = 0. ;
    for (int i = 0; i < N; i++){
        
        bool isElectron = (i == N-1) ;
        double Ti = isElectron ? Te : T ;
        
        const double ni_over_rho = n0[i] / rho;

        double Ei = e_ion[i] ;

        // contributi (per particella) dentro la parentesi:
        const double e_part = (1.5 * KB * Ti) + Ei + (KB * Ti * Ti * dlogQdT_centered[i]) ;
        const double h_part = (2.5 * KB * Ti) + Ei + (KB * Ti * Ti * dlogQdT_centered[i]) ;

        henthalpies[i] = h_part ;

        if (isElectron)
        {
            ee += ni_over_rho * e_part ;
            he += ni_over_rho * h_part ;
        
        }else{
        
            eh += ni_over_rho * e_part ;
            hh += ni_over_rho * h_part ;
        
        }
        
    };
    
    double rhof = 0., rhob = 0.;
    for (int i = 0; i < N; i++){
        rhof += nf[i]*masses[i] ;
        rhob += nb[i]*masses[i] ;
    }

    double ef = 0., eb = 0.;
    for (int i = 0; i < N; i++){
        
        bool isElectron = (i == N-1) ;
        double Tif = isElectron ? Tef : Tf ;
        double Tib = isElectron ? Teb : Tb ;

        double Ei = e_ion[i] ;

        ef += KB * nf[i]* ( (1.5 *  Tif) + Ei/KB + (Tif * dlogQdT_forward[i]) ) ;

        eb += KB * nb[i]* ( (1.5 *  Tib) + Ei/KB + (Tib * dlogQdT_backward[i]) ) ;

    };


    // double dRdT = -(R/T) -(R/rho)*((rhof-rhob)/(2.*dT)) ;

    double cp = (ef/rhof-eb/rhob)/(2.*dT*theta) + R*(1+T/R*dRdT);

    // double dRdP = (R/P) - (R/rho)*((rhopf-rhopb)/(2.*dP)) ;

    double cv = cp - R*((1.+T/R*dRdT)*(1.+T/R*dRdT)/(1.-P/R*dRdP)) ;
    double gamma = cp/cv ;
    double a = std::sqrt( gamma*R*T/( 1.0 - P/R*dRdP ) ) ;
    Td[0] = rho ;
    Td[1] = R ;
    Td[2] = he ;
    Td[3] = hh ;
    Td[4] = ee ;
    Td[5] = eh ;
    Td[6] = cp ;
    Td[7] = cv ;
    Td[8] = gamma ;
    Td[9] = a ;
    dRdP_final = dRdP ; // save for DH corrections in ThermodynamicsDHcorrected

}

void ThermodynamicsDHcorrected::computeThermodynamics(GasMixture& gasmix) {

    Thermodynamics::computeThermodynamics(gasmix) ; 
    
    double T = gasmix.getTemperature();
    double P = gasmix.getPressure();
    double debyeL = gasmix.getCompositionObj()->getDebyeLength(T) ; 

    // Capitelli FACPP pag.104 eq 6.15
    double deltaUDH = -(1./8.)*(KB*T)*(1./( std::numbers::pi * pow(debyeL,3.)));

    // Capitelli FACPP pag. 106 formula 6.29, hentalphy
    Td[3] += (4./3.)*deltaUDH ; 

    // Internal energy
    Td[5] += deltaUDH ; 

    // Capitelli FACPP pag.107 formula 6.34, Cp
    Td[6] -= (4./3.)*deltaUDH/T ; 
    
    // Capitelli FACPP pag.107 formula 6.32, Cv
    Td[7] -= deltaUDH/(2.*T) ; 
    
    Td[8] = Td[6]/Td[7] ; 
    Td[9] = std::sqrt( Td[8] * Td[1] * T / ( 1.0 - P/Td[1] * dRdP_final ) );
    
}