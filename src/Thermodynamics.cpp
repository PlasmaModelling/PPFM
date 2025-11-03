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

    // initial values of dT
    double dT = 20; 
    double dP = 100;

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

    // --- parametri adattivi ---
    double dT0   = 2000.0;
    double dP0   = 10000.0;
    double tolCp = 1e-4;   // tolleranza relativa su Cp
    double tolRP = 1e-4;   // tolleranza relativa su dR/dP
    int    maxIt = 100;

    // --- helper: ricomputa cp, dRdT con un dato dT ---
    auto RecomputeWith_dT = [&](double dT_local,
                                double& cp_out,
                                double& dRdT_out,
                                double& R_out) {
        // ----- Derivative-dependent quantities -----
        double Tf  = T  + dT_local;
        double Tb  = T  - dT_local;
        double Tef = Te + dT_local*theta;
        double Teb = Te - dT_local*theta;

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

        // Partitions @ Tf / Tef
        Q.computePartitionFunctions(Tf, P, gasmix.getCompositionObj()->getDebyeLength(Tf));
        Q[N-1]->computePartitionFunction(Tef, P, gasmix.getCompositionObj()->getDebyeLength(Tef));
        for (int i = 0; i < N; i++) {
            if (dynamic_cast<Electron*>(gasmix(i)) != nullptr) {
                ef += 1.5 * KB * Tef * nf[i];
            } else {
                ef += ((1.5 * KB * Tf + epsf[i]) * nf[i]);
                // derivative of ln Q by central diff: (ln Qf - ln Q) / (2 dT)
                ef += (KB * Tf) * nf[i] * (std::log(Qf(i)) - std::log(Q(i))) / (2.0 * dT_local);
            }
        }

        // Partitions @ Tb / Teb
        Q.computePartitionFunctions(Tb, P, gasmix.getCompositionObj()->getDebyeLength(Tb));
        Q[N-1]->computePartitionFunction(Teb, P, gasmix.getCompositionObj()->getDebyeLength(Teb));
        for (int i = 0; i < N; i++) {
            if (dynamic_cast<Electron*>(gasmix(i)) != nullptr) {
                eb += 1.5 * KB * Teb * nb[i];
            } else {
                eb += ((1.5 * KB * Tb + epsf[i]) * nb[i]);
                // derivative of ln Q by central diff: (ln Q - ln Qb) / (2 dT)
                eb += (KB * Tb) * nb[i] * (std::log(Q(i)) - std::log(Qb(i))) / (2.0 * dT_local);
            }
        }

        // R forward/backward
        double Rf = mass[N-1] * nf[N-1] * Tef;
        double Rb = mass[N-1] * nb[N-1] * Teb;
        for (int i = 0; i < N-1; i++){
            Rf += mass[i]*nf[i]*Tf;
            Rb += mass[i]*nb[i]*Tb;
        }
        Rf = P / Rf;
        Rb = P / Rb;

        double dRdT = (Rf - Rb) / (2.0 * dT_local);

        // GODIN eq.41
        double cp_local = (( ef / rhof ) - ( eb / rhob )) / ( 2.0 * dT_local * theta ) + R * ( 1.0 + ( T / R ) * dRdT );

        // output
        cp_out   = cp_local;
        dRdT_out = dRdT;
        R_out    = R;  // R base (assumo già calcolato a T, P correnti)
    };

    // --- helper: ricomputa dRdP con un dato dP ---
    auto Recompute_dRdP = [&](double dP_local,
                            double& dRdP_out) {
        Gas gasf = gasmix;
        Gas gasb = gasmix;
        gasf.setT(T);
        gasb.setT(T);

        double Pf = P + dP_local;
        double Pb = P - dP_local;

        gasf.setP(Pf);
        auto* solver = gasmix.getCompositionObj();
        solver->setn0(n0);
        solver->CompositionSolve(&gasmix, &gasf);
        std::vector<double> nfP = solver->compositions();

        gasb.setP(Pb);
        solver->setn0(n0);
        solver->CompositionSolve(&gasmix, &gasb);
        std::vector<double> nbP = solver->compositions();

        double Rf = mass[N-1] * nfP[N-1] * Te;
        double Rb = mass[N-1] * nbP[N-1] * Te;
        for (int i = 0; i < N-1; i++){
            Rf += mass[i]*nfP[i]*T;
            Rb += mass[i]*nbP[i]*T;
        }
        Rf = Pf / Rf;
        Rb = Pb / Rb;

        dRdP_out = (Rf - Rb) / (2.0 * dP_local);
    };

    // ================== ADATTIVO SU dT ==================
    double dT_curr = dT0;
    double cp_curr = 0.0, cp_prev = std::numeric_limits<double>::quiet_NaN();
    double dRdT_curr = 0.0, dRdT_prev = 0.0;
    double R_base_dummy = 0.0;

    for (int it = 0; it < maxIt; ++it) {
        RecomputeWith_dT(dT_curr, cp_curr, dRdT_curr, R_base_dummy);

        if (it > 0) {
            double relCp = std::abs((cp_curr - cp_prev) / std::max(std::abs(cp_curr), 1e-30));
            if (relCp < tolCp) break;
        }
        cp_prev   = cp_curr;
        dRdT_prev = dRdT_curr;

        dT_curr *= 0.5; // raffina (top-down)
    }

    // valori finali da usare
    double cp_final   = cp_curr;
    double dRdT_final = dRdT_curr;

    // ================== ADATTIVO SU dP ==================
    double dP_curr = dP0;
    double dRdP_curr = 0.0, dRdP_prev = std::numeric_limits<double>::quiet_NaN();

    for (int it = 0; it < maxIt; ++it) {
        Recompute_dRdP(dP_curr, dRdP_curr);

        if (it > 0) {
            double relRP = std::abs((dRdP_curr - dRdP_prev) / std::max(std::abs(dRdP_curr), 1e-30));
            if (relRP < tolRP) break;
        }
        dRdP_prev = dRdP_curr;

        dP_curr *= 0.5; // raffina (top-down)
    }

    // ================== FINALI: cv, gamma, vs ==================
    cp   = cp_final;
    double dRdT_use = dRdT_final;
    double dRdP_use = dRdP_curr;

    dRdP_final = dRdP_use; 

    // GODIN eq.42 & 43
    cv    = cp - ( R * ( std::pow( ( 1.0 + T/R * dRdT_use ), 2.0 ) / (1.0 - P/R * dRdP_use ) ) );
    gamma = cp / cv;
    vs    = std::sqrt( gamma * R * T / ( 1.0 - P/R * dRdP_use ) );

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