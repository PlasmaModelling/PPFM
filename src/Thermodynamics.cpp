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
    
    PfBox Qf (* gasmix.Comp->Qbox) ;
    PfBox Qb (* gasmix.Comp->Qbox) ;
    PfBox Q  (* gasmix.Comp->Qbox) ;

    double theta = gasmix.theta->get() ;
    double T = gasmix.getTemperature() ;
    double Te = T * theta;
    double P = gasmix.getPressure() ;
    int N =    gasmix.getN() ;

    const double dT = 0.01 ; 
    const double dP = 10. ;

    Qf.computePartitionFunctions(T+dT,P,gasmix.Comp->getDebyeLength(T+dT));
    Q.computePartitionFunctions(T,P,gasmix.Comp->getDebyeLength(T));
    Qb.computePartitionFunctions(T-dT,P,gasmix.Comp->getDebyeLength(T-dT)); 

    Qf[N-1]->computePartitionFunction(Te+dT*theta,P,gasmix.Comp->getDebyeLength(Te+dT)) ; 
    Q[N-1]->computePartitionFunction(Te,P,gasmix.Comp->getDebyeLength(Te)) ; 
    Qb[N-1]->computePartitionFunction(Te-dT*theta,P,gasmix.Comp->getDebyeLength(Te-dT)) ; 

    double rho , R , he , hh , ee , eh , cp , cv , gamma , vs ;

    std::vector<double> n, mass, epsf ;
    n = gasmix.Comp->compositions();

    std::vector<double> n0 = gasmix.Comp->compositions();
    
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
    
    he = 0.; ee = 0.; hh = 0.; eh = 0.;
    for (int i = 0; i < N; i++) {
        if (dynamic_cast<Electron*>(gasmix(i)) != nullptr ) {
            // 1o order
            he = ( 2.5 * KB * Te * n[i] ) / rho ; 
            ee = ( 1.5 * KB * Te * n[i] ) / rho ;
            // 2o order
            hh += ( KB * Te / rho ) * n[i] * ( log(Qf(i)) - log(Qb(i)) ) / (2. * dT) / rho ;
            eh += ( KB * Te / rho ) * n[i] * ( log(Qf(i)) - log(Qb(i)) ) / (2. * dT) / rho ;
        } else {
            // 1o order
            hh += ( ( 2.5 * KB * T + epsf[i] ) * n[i] ) / rho ; 
            eh += ( ( 1.5 * KB * T + epsf[i] ) * n[i] ) / rho ;
            // 2o order
            hh += ( KB * T / rho ) * n[i] * ( log(Qf(i)) - log(Qb(i)) ) / (2. * dT) / rho ;
            eh += ( KB * T / rho ) * n[i] * ( log(Qf(i)) - log(Qb(i)) ) / (2. * dT) / rho ;
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

    Gas gasf = ( gasmix ) , gasb = ( gasmix ) ;
    gasf.setT ( Tf ) ;
    gasb.setT ( Tb ) ;
     
    Composition nf( &gasmix, &gasmix, new PfBox(Q) ), nb( &gasmix, &gasmix, new PfBox(Q) ) ;

    nf.setn0(n0) ; 
    nb.setn0(n0) ;

    nf.compositionSolve ( &gasmix , &gasf ) ;
    nb.compositionSolve ( &gasmix , &gasb ) ;
    
    Qf.computePartitionFunctions(Tf + dT,P,gasmix.Comp->getDebyeLength(Tf + dT));
    Qb.computePartitionFunctions(Tb - dT,P,gasmix.Comp->getDebyeLength(Tb - dT));
    Qf[N-1]->computePartitionFunction(Tef,P,gasmix.Comp->getDebyeLength(Tef));
    Qb[N-1]->computePartitionFunction(Teb,P,gasmix.Comp->getDebyeLength(Teb));

    // rho
    for (int i = 0; i < N; i++) {
        rhof += mass[i] * nf(i) ;
        rhob += mass[i] * nb(i) ;
    }
    
    // e 
    Q.computePartitionFunctions(Tf,P,gasmix.Comp->getDebyeLength(Tf)) ; 
    Q[N-1]->computePartitionFunction(Tef,P,gasmix.Comp->getDebyeLength(Tef)) ;
    for (int i = 0; i < N; i++) {
        if (dynamic_cast<Electron*>(gasmix(i)) != nullptr ) {       
            ef += 1.5 * KB * Tef * nf(i) ;
        } else {
            ef += ( ( 1.5 * KB * Tf + epsf[i] ) * nf(i) ) ;
            ef += ( KB * Tf ) * nf(i) * (log(Qf(i)) - log(Q(i))) / (2. * dT) ;
        }
    }
    // e 
    Q.computePartitionFunctions(Tb,P,gasmix.Comp->getDebyeLength(Tb)) ; 
    for (int i = 0; i < N; i++) {
        if (dynamic_cast<Electron*>(gasmix(i)) != nullptr ) {       
            eb += 1.5 * KB * Teb * nb(i) ;
        } else {
            eb += ( ( 1.5 * KB * Tb + epsf[i] ) * nb(i) ) ;
            eb += ( KB * Tb ) * nb(i) * (log(Q(i)) - log(Qb(i))) / (2. * dT) ;
        }
    }

    // R e dRdT
    Rf = mass[N-1] * nf(N-1) * Tef ;
    Rb = mass[N-1] * nb(N-1) * Teb ;    
    for (int i = 0; i < N-1; i++){
            Rf += mass[i]*nf(i)*Tf;
            Rb += mass[i]*nb(i)*Tb;
    }
    Rf = P/Rf ;
    Rb = P/Rb ;
    double dRdT = (Rf-Rb)/(2. * dT) ;

    // GODIN eq.41
    cp = (( ef / rhof ) - ( eb / rhob )) / ( 2. * dT * theta ) + R * ( 1. + ( T/R ) * dRdT ) ;
    
    // calcolo calore specifico a volume costante
    gasf.setT(T) ;
    gasb.setT(T) ;

    Rf = 0.; Rb = 0.;
    double Pf = P + dP ;
    double Pb = P - dP ;
    gasf.setP ( Pf ) ; 
    nf.setn0(n0) ;
    nf.compositionSolve ( &gasmix , &gasf ) ; 
    gasb.setP ( Pb ) ; 
    nb.setn0(n0) ;
    nb.compositionSolve ( &gasmix , &gasb ) ;
    
    Rf = mass[N-1] * nf(N-1) * Te ;
    Rb = mass[N-1] * nb(N-1) * Te ;    
    for (int i = 0; i < N-1; i++){
            Rf += mass[i]*nf(i)*T;
            Rb += mass[i]*nb(i)*T;
    }
    Rf = Pf / Rf ;
    Rb = Pb / Rb ;
    double dRdP = ( Rf - Rb ) / ( 2. * dP ) ;

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

std::string ThermodynamicsCsv::BuildFileName(const std::string& filename ) const  {

    return "TH_" + filename + ".csv";  
    
}

void ThermodynamicsCsv::PrepareHeader() {
    
    header = "Te [K], ρ [kg/m³], Cₚ [J/(kg·K)], hₑ + hₕ [J/kg], γ [#], a [m/s]";

}
    
void ThermodynamicsCsv::PrintMessage(const std::string& filename)  { 
    std::cout << "Thermodynamic properties " << filename << " printed." << std::endl ;
}

ThermodynamicsCsv::ThermodynamicsCsv () : Thermodynamics() {}

ThermodynamicsCsv::ThermodynamicsCsv ( const std::string& folder ) : Thermodynamics() { customFolder = folder; }

void ThermodynamicsCsv::PrepareData ( const std::vector<double>& temperatureRange, GasMixture* gasmix ) {
        
    double T0 = gasmix->getTemperature() ;  

    double theta = gasmix->theta->get();

    data.resize (temperatureRange.size(), std::vector<double>(6) );

    for (int i = 0; i < temperatureRange.size() ; i++) {
        
        gasmix->setT(temperatureRange[i]) ;        
        computeThermodynamics(*gasmix) ;
    
        data[i][0] = temperatureRange[i]*theta ;
        data[i][1] = Td[0] ;
        data[i][2] = Td[6] ; 
        data[i][3] = Td[2] + Td[3] ; 
        data[i][4] = Td[8] ; 
        data[i][5] = Td[9] ; 
        
    }

    // SETBACK 
    gasmix->setT(T0); 
    gasmix->restartComposition() ; 
    computeThermodynamics(*gasmix) ;

} 
