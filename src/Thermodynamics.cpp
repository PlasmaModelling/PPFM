 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include"Thermodynamics.h"
#include"PfBox.h"
#include"GasMixture.h"
#include"PartitionFunction.h"

void Thermodynamics::computeThermodynamics ( GasMixture& gasmix ) {
    
    PfBox Qf = *(gasmix.getCompositionObj()->getPfBox()) ;
    PfBox Qb = *(gasmix.getCompositionObj()->getPfBox()) ;
    PfBox Q  = *(gasmix.getCompositionObj()->getPfBox()) ;

    double theta = gasmix.theta->get() ;
    double T = gasmix.getTemperature() ;
    double Te = T * theta;
    double P = gasmix.getPressure() ;
    int N =    gasmix.getN() ;
    henthalpies.resize(N);

    const double dT = 50 ; 
    const double dP = 1000. ;

    Qf.computePartitionFunctions(T+dT,P,gasmix.getCompositionObj()->getDebyeLength(T+dT));
    Qf[N-1]->computePartitionFunction(Te+dT*theta,P,gasmix.getCompositionObj()->getDebyeLength(Te+dT)) ; 
    Qf.updateCachedValues();


    Q.computePartitionFunctions(T,P,gasmix.getCompositionObj()->getDebyeLength(T));
    Q[N-1]->computePartitionFunction(Te,P,gasmix.getCompositionObj()->getDebyeLength(Te)) ; 
    Q.updateCachedValues();

    Qb.computePartitionFunctions(T-dT,P,gasmix.getCompositionObj()->getDebyeLength(T-dT)); 
    Qb[N-1]->computePartitionFunction(Te-dT*theta,P,gasmix.getCompositionObj()->getDebyeLength(Te-dT)) ; 
    Qb.updateCachedValues();

    double rho=0.0, R=0.0, he=0.0, hh=0.0, ee=0.0, eh=0.0;
    double cp=0.0, cv=0.0, gamma=0.0, vs=0.0;

    std::vector<double> n, mass, epsf ;
    n = gasmix.getCompositionObj()->compositions();;

    std::vector<double> n0 = gasmix.getCompositionObj()->compositions();;
    
    mass = gasmix.masses();
    rho = 0.;
    for (int i = 0; i < N; i++){
        rho += mass[i]*n[i];
        epsf.push_back(gasmix(i)->formationEnergy()) ;
    }
    
    R = 0.;
    R += mass[N-1]*n[N-1]*Te;
    for (int i = 0; i < N-1; i++)
        R += n[i]*mass[i]*T ;
    R = P/R ;
    
    std::vector<double> ei(N, 0.0);   // contributo specie i a e (J/kg miscela)
    std::vector<double> hi(N, 0.0);   // contributo specie i a h (J/kg miscela)
    for (int i = 0; i < N; i++) {

        bool isElectron = (dynamic_cast<Electron*>(gasmix(i)) != nullptr);
        double Ti = isElectron ? Te : T;

        // d/dT ln(Q_int) (tu stai usando forma (Qf-Qb)/(2 Q dT))
        // NB: meglio proteggere Q(i) da 0/negativi.
        const double Qi = std::max(Q(i), std::numeric_limits<double>::min());
        double dlogQdT = log(Qf(i) / Qb(i)) / (2.0 * dT);

        // fattore n_i / rho
        const double ni_over_rho = n[i] / rho;

        // Ei per particella: epsf[i]
        const double Ei = epsf[i];

        // contributi (per particella) dentro la parentesi:
        const double e_part = (1.5 * KB * Ti) + Ei + (KB * Ti * Ti * dlogQdT);
        const double h_part = (2.5 * KB * Ti) + Ei + (KB * Ti * Ti * dlogQdT);

        // contributi specie i (J/kg miscela)
        ei[i] = ni_over_rho * e_part;
        hi[i] = ni_over_rho * h_part;

        // se vuoi mantenere anche il tuo vettore henthalpies (per particella)
        henthalpies[i] = h_part;

        if (isElectron) {

            ee += ei[i];
            he += hi[i];
        
        } else {
        
            eh += ei[i];
            hh += hi[i];
        
        }
    }

    // derivatives dependent quantities 

    // Specific heat at constant pressure
    double Tf, rhof{0.}, Rf{0.}, hf{0.}, ef{0.};
    double Tb, rhob{0.}, Rb{0.}, hb{0.}, eb{0.};
    
    Tf = T + dT ;
    Tb = T - dT ;
    double Tef = Te + dT*theta ;
    double Teb = Te - dT*theta ;    

    gasmix.setT(Tf) ;
    std::vector<double> nf; 
    nf = gasmix.getCompositionObj()->compositions() ;
    Qf.computePartitionFunctions(Tf,P,gasmix.getCompositionObj()->getDebyeLength(Tf));
    Qf[N-1]->computePartitionFunction(Tef,P,gasmix.getCompositionObj()->getDebyeLength(Tef));
    Qf.updateCachedValues();
    
    gasmix.setT(T);
    gasmix.getCompositionObj()->setn0(n0) ;

    gasmix.setT(Tb) ;
    std::vector<double> nb;
    nb = gasmix.getCompositionObj()->compositions() ;
    Qb.computePartitionFunctions(Tb,P,gasmix.getCompositionObj()->getDebyeLength(Tb));
    Qb[N-1]->computePartitionFunction(Teb,P,gasmix.getCompositionObj()->getDebyeLength(Teb)) ;
    Qb.updateCachedValues();

    gasmix.setT(T);
    gasmix.getCompositionObj()->setn0(n0) ;
    
    // rho
    rhof = 0.;
    rhob = 0.;
    for (int i = 0; i < N; i++) {
        rhof += mass[i] * nf[i] ;
        rhob += mass[i] * nb[i]  ;
    }
    
    // e 
    ef = 0.;
    for (int i = 0; i < N; i++) {
        if (dynamic_cast<Electron*>(gasmix(i)) != nullptr ) {       
            ef += 1.5 * KB * Tef * nf[i] ;
        } else {
            ef += ( ( 1.5 * KB * Tf + epsf[i] ) * nf[i] ) ;
            ef += ( KB * Tf * Tf) * nf[i] * (log(Qf(i)/Q(i))) / (2. * dT) ;
        }
    }
    
    // e 
    eb = 0.;
    for (int i = 0; i < N; i++) {
        if (dynamic_cast<Electron*>(gasmix(i)) != nullptr ) {       
            eb += 1.5 * KB * Teb * nb[i] ;
        } else {
            eb += ( ( 1.5 * KB * Tb + epsf[i] ) * nb[i] ) ;
            eb += ( KB * Tb  * Tb) * nb[i] * (log(Q(i)/Qb(i))) / (2. * dT) ;
        }
    }

    // R e dRdT
    Rf = mass[N-1] * nf[N-1] * Tef ;
    Rb = mass[N-1] * nb[N-1] * Teb ;    
    for (int i = 0; i < N-1; i++){
            Rf += mass[i]*nf[i]*Tf;
            Rb += mass[i]*nb[i]*Tb;
    }
    Rf = P/Rf ;
    Rb = P/Rb ;
    double dRdT = (Rf-Rb)/(2. * dT) ;

    // GODIN eq.41
    cp = (( ef / rhof ) - ( eb / rhob )) / ( 2. * dT * theta ) + R * ( 1. + ( T/R ) * dRdT ) ;
    
    // calcolo calore specifico a volume costante
    
    Rf = 0.; Rb = 0.;
    double Pf = P + dP ;
    double Pb = P - dP ;
    
    gasmix.setP(Pf) ;
    nf = gasmix.getCompositionObj()->compositions() ;
    Qf.computePartitionFunctions(T,Pf,gasmix.getCompositionObj()->getDebyeLength(T));
    Qf[N-1]->computePartitionFunction(Te,Pf,gasmix.getCompositionObj()->getDebyeLength(Te)) ;
    Qf.updateCachedValues();

    gasmix.setP(P);
    gasmix.getCompositionObj()->setn0(n0) ;

    gasmix.setP(Pb) ;
    nb = gasmix.getCompositionObj()->compositions() ;
    Qb.computePartitionFunctions(T,Pb,gasmix.getCompositionObj()->getDebyeLength(T));
    Qb[N-1]->computePartitionFunction(Te,Pb,gasmix.getCompositionObj()->getDebyeLength(Te)) ;
    Qb.updateCachedValues();
    
    gasmix.setP(P);
    gasmix.getCompositionObj()->setn0(n0) ;

    Rf = mass[N-1] * nf[N-1] * Te ;
    Rb = mass[N-1] * nb[N-1] * Te ;    
    for (int i = 0; i < N-1; i++){
            Rf += mass[i]*nf[i]*T;
            Rb += mass[i]*nb[i]*T;
    }
    Rf = Pf / Rf ;
    Rb = Pb / Rb ;
    double dRdP = ( Rf - Rb ) / ( 2. * dP ) ;

    dRdP_final = dRdP ;

    // GODIN eq.42
    cv = cp - ( R * ( pow(( 1. + T/R * dRdT ),2.) / (1. - P/R *dRdP ) ));
    
    // eq.43 Cp/Cv Velocità del suono
    gamma = cp/cv ; 
    vs = sqrt( gamma * R * T / ( 1. - P/R * dRdP )) ;

    Td[0] = rho ;
    Td[1] = R ;
    Td[2] = he ;
    Td[3] = hh ;
    Td[4] = ee ;
    Td[5] = eh ;
    Td[6] = cp ;
    Td[7] = cv ;
    Td[8] = gamma ;
    Td[9] = vs ;
}

void ThermodynamicsDHcorrected::computeThermodynamics(GasMixture& gasmix) {

    Thermodynamics::computeThermodynamics(gasmix) ; 
    
    double T = gasmix.getTemperature();
    double P = gasmix.getPressure();
    double debyeL = gasmix.getCompositionObj()->getDebyeLength(T) ; 

    // Capitelli FACPP pag.104 eq 6.15
    double deltaUDH = -(1./8.)*(KB*T)*(1./(std::numbers::pi*pow(debyeL,3.)));

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