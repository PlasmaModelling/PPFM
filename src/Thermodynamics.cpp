#include "Thermodynamics.h"
#include "PfBox.h"
#include "GasMixture.h"
#include "PartitionFunction.h"
#include "Composition.h"

void Thermodynamics::computeThermodynamics(GasMixture& gasmix) {
    
    PfBox Qf(*gasmix.getCompositionObj()->getPfBox());
    PfBox Qb(*gasmix.getCompositionObj()->getPfBox());
    PfBox Q (*gasmix.getCompositionObj()->getPfBox());

    double theta = gasmix.theta->get();
    double T  = gasmix.getTemperature();
    double Te = T * theta;
    double P  = gasmix.getPressure();
    int N     = gasmix.getN();

    const double dT = 0.01; 
    const double dP = 10.0;

    // Partition functions at shifted states
    Qf.computePartitionFunctions(T + dT, P, gasmix.getCompositionObj()->getDebyeLength(T + dT));
    Q .computePartitionFunctions(T,      P, gasmix.getCompositionObj()->getDebyeLength(T));
    Qb.computePartitionFunctions(T - dT, P, gasmix.getCompositionObj()->getDebyeLength(T - dT));

    Qf[N-1]->computePartitionFunction(Te + dT*theta, P, gasmix.getCompositionObj()->getDebyeLength(Te + dT));
    Q [N-1]->computePartitionFunction(Te,            P, gasmix.getCompositionObj()->getDebyeLength(Te));
    Qb[N-1]->computePartitionFunction(Te - dT*theta, P, gasmix.getCompositionObj()->getDebyeLength(Te - dT));

    double rho , R , he , hh , ee , eh , cp , cv , gamma , vs ;

    std::vector<double> n, mass, epsf;
    n  = gasmix.getCompositionObj()->compositions();
    std::vector<double> n0 = n;
    
    mass = gasmix.masses();
    rho = 0.;
    for (int i = 0; i < N; i++) {
        rho += mass[i] * n[i];
        epsf.push_back(gasmix(i)->formationEnergy());
    }

    R = mass[N-1] * n[N-1] * Te;
    for (int i = 0; i < N-1; i++)
        R += n[i] * mass[i] * T;
    R = P / R;
    
    // --- Energetics ---
    he = hh = ee = eh = 0.0;

    static bool ref_set = false;
    static double e_chem_ref = 0.0;

    double e_chem_tot = 0.0;

    for (int i = 0; i < N; i++) {
        bool isElectron = (dynamic_cast<Electron*>(gasmix(i)) != nullptr);
        double Ti = isElectron ? Te : T;

        double e_tr = 1.5 * KB * Ti * n[i] / rho;
        double e_ch = epsf[i] * n[i] / rho;
        e_chem_tot += e_ch;

        double dlogQdT = (std::log(Qf(i)) - std::log(Qb(i))) / (2.0 * dT);
        double e_in = KB * Ti * Ti * dlogQdT * n[i] / rho;

        double e_tot = e_tr + e_ch + e_in;
        double h_i   = (KB * Ti * n[i]) / rho + e_tot;

        if (isElectron) {
            ee += e_tot;
            he += h_i;
        } else {
            eh += e_tot;
            hh += h_i;
        }
    }

    if (!ref_set) {
        e_chem_ref = e_chem_tot;
        ref_set = true;
    }

    double delta_e_chem = e_chem_tot - e_chem_ref;
    hh += (delta_e_chem - e_chem_tot);
    eh += (delta_e_chem - e_chem_tot);

    // ----- Derivative-dependent quantities -----
    double Tf = T + dT;
    double Tb = T - dT;
    double Tef = Te + dT*theta;
    double Teb = Te - dT*theta;    

    Gas gasf = gasmix;
    Gas gasb = gasmix;
    gasf.setT(Tf);
    gasb.setT(Tb);

    auto* solver = gasmix.getCompositionObj();

    // Forward composition
    solver->setn0(n0);
    solver->CompositionSolve(&gasmix, &gasf);
    std::vector<double> nf = solver->compositions();

    // Backward composition
    solver->setn0(n0);
    solver->CompositionSolve(&gasmix, &gasb);
    std::vector<double> nb = solver->compositions();

    // rho forward/backward
    double rhof = 0.0, rhob = 0.0;
    for (int i = 0; i < N; i++) {
        rhof += mass[i] * nf[i];
        rhob += mass[i] * nb[i];
    }

    // Energy forward/backward for Cp
    double ef = 0.0, eb = 0.0;
    Q.computePartitionFunctions(Tf,P,gasmix.getCompositionObj()->getDebyeLength(Tf));
    Q[N-1]->computePartitionFunction(Tef,P,gasmix.getCompositionObj()->getDebyeLength(Tef));
    for (int i = 0; i < N; i++) {
        if (dynamic_cast<Electron*>(gasmix(i)) != nullptr) {
            ef += 1.5 * KB * Tef * nf[i];
        } else {
            ef += ((1.5 * KB * Tf + epsf[i]) * nf[i]);
            ef += (KB * Tf) * nf[i] * (log(Qf(i)) - log(Q(i))) / (2. * dT);
        }
    }
    Q.computePartitionFunctions(Tb,P,gasmix.getCompositionObj()->getDebyeLength(Tb));
    for (int i = 0; i < N; i++) {
        if (dynamic_cast<Electron*>(gasmix(i)) != nullptr) {
            eb += 1.5 * KB * Teb * nb[i];
        } else {
            eb += ((1.5 * KB * Tb + epsf[i]) * nb[i]);
            eb += (KB * Tb) * nb[i] * (log(Q(i)) - log(Qb(i))) / (2. * dT);
        }
    }

    // R forward/backward
    double Rf = mass[N-1] * nf[N-1] * Tef;
    double Rb = mass[N-1] * nb[N-1] * Teb;    
    for (int i = 0; i < N-1; i++){
        Rf += mass[i]*nf[i]*Tf;
        Rb += mass[i]*nb[i]*Tb;
    }
    Rf = P/Rf;
    Rb = P/Rb;
    double dRdT = (Rf - Rb)/(2. * dT);

    // GODIN eq.41
    cp = (( ef / rhof ) - ( eb / rhob )) / ( 2. * dT * theta ) + R * ( 1. + ( T/R ) * dRdT );

    // --- Cv from pressure perturbation ---
    gasf.setT(T);
    gasb.setT(T);

    double Pf = P + dP;
    double Pb = P - dP;

    gasf.setP(Pf);
    solver->setn0(n0);
    solver->CompositionSolve(&gasmix, &gasf);
    std::vector<double> nfP = solver->compositions();

    gasb.setP(Pb);
    solver->setn0(n0);
    solver->CompositionSolve(&gasmix, &gasb);
    std::vector<double> nbP = solver->compositions();

    Rf = mass[N-1] * nfP[N-1] * Te;
    Rb = mass[N-1] * nbP[N-1] * Te;    
    for (int i = 0; i < N-1; i++){
        Rf += mass[i]*nfP[i]*T;
        Rb += mass[i]*nbP[i]*T;
    }
    Rf = Pf / Rf;
    Rb = Pb / Rb;
    double dRdP = (Rf - Rb)/(2. * dP);

    // GODIN eq.42 & 43
    cv = cp - ( R * ( pow(( 1. + T/R * dRdT ),2.) / (1. - P/R * dRdP ) ));
    gamma = cp/cv;
    vs = sqrt( gamma * R * T / ( 1. - P/R * dRdP ));

    // Store results
    Td[0] = rho;
    Td[1] = R;
    Td[2] = he;
    Td[3] = hh;
    Td[4] = ee;
    Td[5] = eh;
    Td[6] = cp;
    Td[7] = cv;
    Td[8] = gamma;
    Td[9] = vs;
}

// ---------------- ThermodynamicsCsv ---------------- //

std::string ThermodynamicsCsv::BuildFileName(const std::string& filename) const {
    return "TH_" + filename + ".csv";  
}

void ThermodynamicsCsv::PrepareHeader() {
    header = "Te [K], ρ [kg/m³], Cₚ [J/(kg·K)], hₑ + hₕ [J/kg], γ [–], a [m/s]";
}
    
void ThermodynamicsCsv::PrintMessage(const std::string& filename) {
    std::cout << "Thermodynamic properties " << filename << " printed." << std::endl;
}

ThermodynamicsCsv::ThermodynamicsCsv(Thermodynamics* solver) : solver(solver) {}

ThermodynamicsCsv::ThermodynamicsCsv(Thermodynamics* solver, const std::string& folder) : solver(solver) {
    customFolder = folder;
}

void ThermodynamicsCsv::PrepareData(const std::vector<double>& temperatureRange, GasMixture* gasmix) {
        
    double T0 = gasmix->getTemperature();  
    double theta = gasmix->theta->get();

    data.resize(temperatureRange.size(), std::vector<double>(6));

    for (int i = 0; i < temperatureRange.size(); i++) {
        gasmix->setT(temperatureRange[i]);        
        solver->computeThermodynamics(*gasmix);
    
        data[i][0] = temperatureRange[i] * theta;
        data[i][1] = solver->rho();
        data[i][2] = solver->cp();
        data[i][3] = solver->he() + solver->hh();
        data[i][4] = solver->gamma();
        data[i][5] = solver->a();
    }

    gasmix->setT(T0); 
    gasmix->restartComposition(); 
    solver->computeThermodynamics(*gasmix);
} 
