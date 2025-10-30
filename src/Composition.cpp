// PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
// (University of Bologna, Italy)                   // 
// Licensed under CC BY 4.0.                        // 
// To view a copy of this license, visit:           // 
// https://creativecommons.org/licenses/by/4.0/     // 

#include <stdexcept>
#include <numbers>
#include <typeindex>
#include <cmath>
#include "Composition.h"
#include "PfBox.h"
#include "PartitionFunction.h"
#include "GasMixture.h"

// __________________________ Composition (abstract) __________________________ //

Composition::Composition(Mixture* mix, Gas* gas)
    : Composition(mix, gas, new PfBox(mix)) { initializeComposition(); }

Composition::Composition(Mixture* mix, Gas* gas, PfBox* qbox)
    : mixptr(mix), gasptr(gas), Qbox(qbox) {

        initializeComposition();
        
        epsf.resize(mix->getN(), 0.0);
        for (int i = 0; i < mix->getN(); i++) 
            epsf[i] = (*mix)(i)->formationEnergy();


}

std::vector<double> Composition::compositions(double conversion) {
    
    std::vector<double> ns;
    ns.reserve(ni.size());
    for (double val : ni) ns.push_back(val * conversion);
    return ns;

}

double Composition::ntot() {

    double sum = 0.0;
    for (double val : ni) sum += val;
    return sum;

}

void Composition::setDebyeModel(const std::string& modelName) {

    if (modelName == "Rat2002Th") 
        debyeChoice = DebyeModel::Rat2002Th;
    else if (modelName == "Rat2002Te") 
        debyeChoice = DebyeModel::Rat2002Te;
    else if (modelName == "Ghourui") 
        debyeChoice = DebyeModel::Ghourui;
    else {
    
        throw std::invalid_argument(
            "Invalid Debye model name: " + modelName +
            ". Please refer to the available options listed in "
            "Composition::DebyeModel in Composition.h"
        );
    
    }
}

double Composition::getDebyeLength(double T) {

    double lambdaD;

    switch (debyeChoice) {
        
        case DebyeModel::Rat2002Th:
            lambdaD = Debye_Rat2002(T);
            break;

        case DebyeModel::Rat2002Te:
            lambdaD = Debye_Rat2002(T * gasptr->theta->get());
            break;

        case DebyeModel::Ghourui:
            lambdaD = Debye_Ghourui(T);
            break;

        default:
            throw std::runtime_error("Unknown Debye model selected");
    }

    // --- ion-sphere cutoff ----------------------------------------------
    /* double nch = 0.0;
    
    int N = mixptr->getN();

    for (int i = 0; i < N; ++i) {
        Species* sp = (*mixptr)(i);
        const int z = sp->getCharge();
        if (z != 0)
            nch += ni[i];
    }

    if (nch > 0.0) {
        double a = std::pow(1.0 / (8.0 * std::numbers::pi * nch), 1.0 / 3.0);
        lambdaD = std::max(lambdaD, a);
    }
    */
    return lambdaD;
}


double Composition::Debye_Rat2002(double T) {

    return std::sqrt((eps0 * KB * T) / (qe * qe * ni.back()));

}

double Composition::Debye_Ghourui(double T) {

    double sum = 0.0;
    int M = mixptr->getM();
    
    for (size_t i = 0; i < ni.size(); i++)
        sum += ( std::pow((*mixptr)(i)->getCharge(), 2.0) * ni[i] ) ; 
    
    sum *= (qe * qe) / (eps0 * KB * T );

    return std::sqrt(1.0 / sum);
}

void Composition::initializeComposition() {

    int N = mixptr->getN();
    int M = mixptr->getM();
    double P = gasptr->getPressure();
    double T = gasptr->getTemperature();

    std::vector<Species*> species;
    for (int i = 0; i < N; i++)
        species.push_back((*mixptr)(i));

    ni.resize(N, 0.0);
    for (int i = 0; i < N; i++) {
        if (species[i]->getFormula() == "e-")
            ni[i] = 1.e15;
        else if (mixptr->isBase(species[i]))
            ni[i] = P / (KB * T);
        else
            ni[i] = 1.e10;
    }

}

std::vector<double> Composition::totalPartitionFunctions(double T, double P, double lambdaD) {
    
    int N = mixptr->getN();
    std::vector<double> Qtot(N, 0.0);

    Qbox->computePartitionFunctions(T, P, lambdaD);

    for (int i = 0; i < N - 1; i++) {

        if (dynamic_cast<BiatomicMolecule*>((*mixptr)(i))) {

            (*Qbox)[i]->computePartitionFunction(T, P, lambdaD);

            Qtot[i] = std::pow(((2.0 * std::numbers::pi * (*mixptr)(i)->getMass() * KB *
                T) / (hPlanck * hPlanck)), 1.5) * std::exp(-((epsf[i]) / (KB * T))) *
                (*Qbox)(i);

        } else {

            Qtot[i] = std::pow(((2.0 * std::numbers::pi * (*mixptr)(i)->getMass() * KB *
                T) / (hPlanck * hPlanck)), 1.5) * std::exp(-((epsf[i]) / (KB * T *
                    gasptr->theta->get()))) *
                (*Qbox)(i);

        }
    }

    Qtot[N - 1] = std::pow(((2. * std::numbers::pi * (*mixptr)(N - 1)->getMass() * KB * T *
        gasptr->theta->get()) /
        (hPlanck * hPlanck)), 1.5) * (*Qbox)(N - 1);

    return Qtot;

}

// __________________________ GodinTrepSahaSolver __________________________ //

GodinTrepSahaSolver::GodinTrepSahaSolver(Mixture* mix, Gas* gas)
    : Composition(mix, gas) {

    C = CompositionMatrix(*mix);

    Peff = gas->getPressure() ;

}

GodinTrepSahaSolver::GodinTrepSahaSolver(Mixture* mix, Gas* gas, PfBox* qbox)
    : Composition(mix, gas, qbox) {

    C = CompositionMatrix(*mix);

    Peff = gas->getPressure() ;

}

std::vector<double> GodinTrepSahaSolver::Crow(Species* specie,
    const std::map<std::type_index, int>& colmap) {

    std::vector<double> row(colmap.size() + 1, 0.0);

    Element* elem;
    if ((elem = dynamic_cast<Element*>(specie))) {

        if (colmap.find(tipo(elem)) != colmap.end())
            row[colmap.at(tipo(elem))] = 1.0;

    } else if (ChargedSpecies* chrgd = dynamic_cast<ChargedSpecies*>(specie)) {

        elem = chrgd->Constituent();
        if ((elem != nullptr) && (colmap.find(tipo(elem)) != colmap.end()))
            row[colmap.at(tipo(elem))] = 1.0;

    } else if (PolyAtomicMolecule* poly = dynamic_cast<PolyAtomicMolecule*>(specie)) {

        for (int i = 0; i < poly->numberOfCostituents(); i++) {

            elem = (*poly)[i];
            int part = (*poly)(i);
            if (colmap.find(tipo(elem)) != colmap.end())
                row[colmap.at(tipo(elem))] = (double)part;

        }

    } else {

        throw std::invalid_argument("Unrecognized species type");

    }

    row.back() = specie->getCharge(); // charge last
    return row;

}

std::vector<std::vector<double>> GodinTrepSahaSolver::CompositionMatrix(Mixture& mixx) {

    int N = mixx.getN();
    int M = mixx.getM();

    std::map<std::type_index, int> columns;

    int k = 0;
    for (int i = 0; i < N; i++)
        if (auto elem = dynamic_cast<Element*>(mixx(i)))
            columns[tipo(elem)] = k++;

    std::vector<std::vector<double>> CC(N, std::vector<double>(M));
    for (int i = 0; i < N; i++)
        CC[i] = Crow(mixx(i), columns);

    return CC;

}

std::vector<std::vector<double>> GodinTrepSahaSolver::ConservationMatrix(
    Mixture& mixx, Gas& gass, const std::vector<std::vector<double>>& C) {

    double theta = gass.theta->get();
    int N = mixx.getN();
    int M = mixx.getM();

    std::map<std::type_index, double> nuclei;
    for (int i = 0; i < N; i++)
        if (auto elem = dynamic_cast<Element*>(mixx(i)))
            nuclei[tipo(elem)] = 0.0;

    std::vector<Species*> mfSpecies;
    std::vector<double> molefractions;
    std::tie(mfSpecies, molefractions) = mixx.getMoleFractions();

    for (int i = 0; i < mfSpecies.size(); i++) {

        if (PolyAtomicMolecule* ptr = dynamic_cast<PolyAtomicMolecule*>(mfSpecies[i])) {
        
            for (size_t j = 0; j < ptr->numberOfCostituents(); j++) {
        
                auto elem = (*ptr)[j];
                double coeff = (*ptr)(j);
                nuclei[tipo(elem)] = molefractions[i] * coeff;
        
            }
        
        } else if (Element* ptr = dynamic_cast<Element*>(mfSpecies[i])) {

            nuclei[tipo(ptr)] = molefractions[i];

        } else {

            throw std::invalid_argument("Invalid Specie set as molefraction.");

        }
    }

    std::vector<double> stechiometry;
    for (const auto& pair : nuclei)
        stechiometry.push_back(pair.second);

    double T = gass.getTemperature();
    
    std::vector<std::vector<double>> AA(M, std::vector<double>(N));

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
    
            if (i < M - 2) {
    
                AA[i][j] = stechiometry[i + 1] * C[j][i] - stechiometry[i] * C[j][i + 1];
    
            } else if (i == M - 2) {
    
                AA[i][j] = C[j][M - 1]; // neutrality
    
            } else {
    
                AA[i][j] = KB * T;
                if (j == N - 1)
                    AA[i][j] = KB * T * theta;
    
            }
        }
    }
    
    return AA;

}

void GodinTrepSahaSolver::baseCalc(std::vector<int>& b, std::vector<int>& bs,
    const std::vector<std::vector<double>>& C) {

    int N = b.size();
    int M = b.size() - bs.size();
    int* perm = new int[N];
    double* n_ord = new double[N];

    for (int i = 0; i < N; i++) n_ord[i] = ni[i];
    sort(n_ord, perm, N);

    int tmp = 1;
    for (int i = 0; i < M; i++) {

        int j = 0;
        b[i] = perm[N - 1];
        if (C[perm[N - 1]][i] == 0) tmp = 0;

        for (int l = 0; l < i; l++)
            if (b[l] == b[i]) tmp = 0;

        while (tmp == 0) {

            j += 1;
            b[i] = perm[N - 1 - j];
            if (C[perm[N - 1 - j]][i] != 0) tmp = 1;

            for (int l = 0; l < i; l++)
                if (b[l] == b[i]) tmp = 0;

        }
    }

    int l = 0;
    for (int i = 0; i < N; i++) {

        tmp = 0;
        for (int j = 0; j < M; j++)
            if (i == b[j]) tmp = 1;

        if (tmp == 0) {

            bs[l] = i;
            l += 1;

        }
    }

    delete[] perm;
    delete[] n_ord;
}

void GodinTrepSahaSolver::CompositionSolve(Mixture* mix, Gas* gas) {
    
    double theta = gas->theta->get() ; 

    int N = mix->getN();
    int M = mix->getM();
    int L = N-M ;
    double T = gas->getTemperature() ; 
    
    std::vector<Species*> species ; 
    for (int i = 0; i < N; i++)
        species.push_back((*mix)(i)) ; 
        
    std::vector<std::vector<double>> B (M,std::vector<double>(M));
    std::vector<std::vector<double>> Bs (L,std::vector<double>(M));
    std::vector<std::vector<double>> Binv (M,std::vector<double>(M));
    std::vector<std::vector<double>> v (L,std::vector<double>(M));
    std::vector<int> b(N);
    std::vector<int> bs(L);
    std::vector<double> A0(M);


    baseCalc ( b, bs, C );
    for(int i = 0;i < M; i++) {
        for(int j = 0; j < M; j++) {
            B[i][j] = C[b[i]][j];
        }
    }
    for(int i=0;i<L;i++) {
        for(int j=0;j<M;j++) {
            Bs[i][j] = C[bs[i]][j];
        }
    }
    lu_inv(Binv,B,M);
    matrix_prod(v,Bs,Binv,L,M);

    std::vector<std::vector<double>> A = ConservationMatrix(*mix,*gas,C) ;
    
    for( int j = 0 ; j < M-1 ; j++) 
        A0[j] = 0.;
    A0[M-1] = Peff ;

    // algoritmo di calcolo
    
    double dnmax, ntot, prod;
    int iter;
    const int itermax = 10000;
    const double ERR = 1.e-15;
    std::vector<double> dn(M,0.0) ;
    std::vector<double> R(M,0.0) ;
    std::vector<double> nold(N,0.0) ;
    std::vector<std::vector<double>> J(M,std::vector<double>(M)) ;
    const double nmin = 1. ;
    dnmax = ERR + 1.;
    iter = 0;
    
    while(dnmax>ERR && iter<itermax){
        
        /* residuals as in - GODIN: eq 34 */
        residual(R,J,ni,A,A0,v,b,bs,N,M);
        for(int i=0; i<M; i++){
            R[i] = - R[i];
        }
        
        /* dn as in GODIN - eq 33 */
        lu_sistema(dn,J,R,M);

        /*Bdet = lu_det(J,M); */
        for(int i=0; i<N; i++) {
            nold[i] = ni[i];
        } 				
        
        /*update base densities GODIN eq.30*/
        for(int i = 0; i < M; i++ ) 
            ni[b[i]] += dn[i];
            
        for(int i = 0; i < M; i++ ) {
 
            if (ni[b[i]] < nmin) 
                ni[b[i]] = nmin;
            
        }
        
        // calcolo funzioni di partizione totali
        std::vector<double> Qtot = totalPartitionFunctions(T,Peff, getDebyeLength(T));
        
        /*update not-base densities: GODIN - eq 30*/
        for(int j = 0; j < (N-M); j++){
            
            prod = 1.;
            
            for(int i = 0; i < M; i++){

                prod *= pow((ni[b[i]]),v[j][i]) *
                    pow(1./Qtot[b[i]],v[j][i]);
            
            }
            
            ni[bs[j]] = prod * Qtot[bs[j]];
            if(ni[bs[j]] < nmin){ ni[bs[j]] = nmin;}
        
        }

        /* convergence of the method */
        ntot = 0.;	
        
        for(int i=0; i<N; i++){nold[i] = ni[i]-nold[i];}
        for(int i=0; i<N; i++){ntot += ni[i]; }
        
        dnmax = max_double(nold,N)/ntot ;
        iter += 1; 
    
    }
}

void GodinTrepSahaSolver::CompositionSolveLambdaFrozen(Mixture* mix, Gas* gas, double lambda) {
    
    double theta = gas->theta->get() ; 

    int N = mix->getN();
    int M = mix->getM();
    int L = N-M ;
    double T = gas->getTemperature() ; 
    
    std::vector<Species*> species ; 
    for (int i = 0; i < N; i++)
        species.push_back((*mix)(i)) ; 
        
    std::vector<std::vector<double>> B (M,std::vector<double>(M));
    std::vector<std::vector<double>> Bs (L,std::vector<double>(M));
    std::vector<std::vector<double>> Binv (M,std::vector<double>(M));
    std::vector<std::vector<double>> v (L,std::vector<double>(M));
    std::vector<int> b(N);
    std::vector<int> bs(L);
    std::vector<double> A0(M);


    baseCalc ( b, bs, C );
    for(int i = 0;i < M; i++) {
        for(int j = 0; j < M; j++) {
            B[i][j] = C[b[i]][j];
        }
    }
    for(int i=0;i<L;i++) {
        for(int j=0;j<M;j++) {
            Bs[i][j] = C[bs[i]][j];
        }
    }
    lu_inv(Binv,B,M);
    matrix_prod(v,Bs,Binv,L,M);

    std::vector<std::vector<double>> A = ConservationMatrix(*mix,*gas,C) ;
    
    for( int j = 0 ; j < M-1 ; j++) 
        A0[j] = 0.;
    A0[M-1] = Peff ;

    // algoritmo di calcolo
    
    double dnmax, ntot, prod;
    int iter;
    const int itermax = 10000;
    const double ERR = 1.e-15;
    std::vector<double> dn(M,0.0) ;
    std::vector<double> R(M,0.0) ;
    std::vector<double> nold(N,0.0) ;
    std::vector<std::vector<double>> J(M,std::vector<double>(M)) ;
    const double nmin = 1. ;
    dnmax = ERR + 1.;
    iter = 0;
    
    while(dnmax>ERR && iter<itermax){
        
        /* residuals as in - GODIN: eq 34 */
        residual(R,J,ni,A,A0,v,b,bs,N,M);
        for(int i=0; i<M; i++){
            R[i] = - R[i];
        }
        
        /* dn as in GODIN - eq 33 */
        lu_sistema(dn,J,R,M);

        /*Bdet = lu_det(J,M); */
        for(int i=0; i<N; i++) {
            nold[i] = ni[i];
        } 				
        
        /*update base densities GODIN eq.30*/
        for(int i = 0; i < M; i++ ) 
            ni[b[i]] += dn[i];
            
        for(int i = 0; i < M; i++ ) {
 
            if (ni[b[i]] < nmin) 
                ni[b[i]] = nmin;
            
        }
        
        // calcolo funzioni di partizione totali
        std::vector<double> Qtot = totalPartitionFunctions(T,Peff, lambda);
        
        /*update not-base densities: GODIN - eq 30*/
        for(int j = 0; j < (N-M); j++){
            
            prod = 1.;
            
            for(int i = 0; i < M; i++){

                prod *= pow((ni[b[i]]),v[j][i]) *
                    pow(1./Qtot[b[i]],v[j][i]);
            
            }
            
            ni[bs[j]] = prod * Qtot[bs[j]];
            if(ni[bs[j]] < nmin){ ni[bs[j]] = nmin;}
        
        }

        /* convergence of the method */
        ntot = 0.;	
        
        for(int i=0; i<N; i++){nold[i] = ni[i]-nold[i];}
        for(int i=0; i<N; i++){ntot += ni[i]; }
        
        dnmax = max_double(nold,N)/ntot ;
        iter += 1; 
    
    }
}

// __________________________ GTSahaDHcorrection __________________________ //

GTSahaDHcorrection::GTSahaDHcorrection(Mixture* mix, Gas* gas)
    : GodinTrepSahaSolver(mix, gas) {

    setDebyeModel("Ghourui");

}  

GTSahaDHcorrection::GTSahaDHcorrection(Mixture* mix, Gas* gas, PfBox* qbox)
    : GodinTrepSahaSolver(mix, gas, qbox) {

    setDebyeModel("Ghourui");

}   

double GTSahaDHcorrection::pressureCorrection(double T, double lambdaD) {
    
    // Capitelli Fundamental Aspects of Chemical Plasma Physics, 2016,
    // Eq. 6.23-6.24, pp. 105
    // P_DH = -k_B T / (24π λ_D^3)
    return - (KB * T) / (48.0 * std::numbers::pi * std::pow(lambdaD, 3));

}

std::vector<double> GTSahaDHcorrection::formationEnergyDHcorrected (double lambdaD) {

    // Capitelli Fundamental Aspects of Chemical Plasma Physics, 2016,
    // Eq. 6.45 pp. 109
    double KDH = (qe * qe) / (8. * std::numbers::pi * eps0 * lambdaD);
    
    // Eq. 6.44 pp. 109
    std::vector<double> delta_epsf;
    for (size_t i = 0; i < mixptr->getN()-1; i++)
        delta_epsf.push_back(KDH * std::pow((*mixptr)(i)->getCharge()+1, 2.));
        
    std::vector<double> corrected_epsf;
    for (size_t i = 0; i < mixptr->getN()-1; i++)
        corrected_epsf.push_back(epsf[i] - delta_epsf[i]);

    return corrected_epsf;

}
void GTSahaDHcorrection::CompositionSolve(Mixture* mix, Gas* gas) {

    double T = gas->getTemperature();
    double Pid = gas->getPressure();

    // Backup delle energie di formazione originali
    std::vector<double> epsf0(mix->getN(), 0.0);
    for (int i = 0; i < mix->getN(); i++) 
        epsf0[i] = (*mix)(i)->formationEnergy();

    // Inizializzazione variabili Debye-Hückel
    double lambdaOld = getDebyeLength(T);
    double lambdaNew = 0.0;
    double deltaLambda = 1.0;

    // Variabili per la convergenza sulla pressione efficace
    double PeffOld = gas->getPressure();
    double PeffNew = PeffOld;
    double deltaPeff = 1.0;

    int iter = 1;
    const double relax = 0.3; // coefficiente di rilassamento (0 < relax <= 1)

    while ((deltaLambda > tol && deltaPeff > tol) && iter <= maxIter) {

        // Ripristina le energie di formazione originali e la pressione ideale
        epsf = epsf0;
        Peff = Pid;

        // Correzione DH delle energie di formazione
        epsf = formationEnergyDHcorrected(lambdaOld);

        // Calcola correzione DH sulla pressione (negativa)
        double PDH = pressureCorrection(T, lambdaOld);
        PeffNew = Pid + PDH;  // Peff = P + PDH

        // Rilassamento sulla pressione per stabilità numerica
        double PeffMixed = (1.0 - relax) * PeffOld + relax * PeffNew;
        Peff = PeffMixed;

        // Risolve la composizione con le energie corrette
        GodinTrepSahaSolver::CompositionSolveLambdaFrozen(mix, gas, lambdaOld);

        // Aggiorna lunghezza di Debye
        lambdaNew = getDebyeLength(T);
        double lambdaMixed = (1.0 - relax) * lambdaOld + relax * lambdaNew;

        // Calcola variazioni relative
        deltaLambda = std::abs((lambdaMixed - lambdaOld) / lambdaOld);
        deltaPeff   = std::abs((PeffMixed - PeffOld) / Pid);

        // Aggiorna per iterazione successiva
        lambdaOld = lambdaMixed;
        PeffOld   = PeffMixed;
        iter += 1;

        // Debug log
        if (iter>=maxIter) {

            std::cout << "[Iter " << iter << "] T=" << T << " K "
            << "λ_D=" << lambdaMixed << " m, "
            << "Δλ/λ=" << deltaLambda << ", "
            << "Peff=" << Peff << " Pa, "
            << "ΔP/Pid=" << deltaPeff
            << std::endl;
        
        }
       
    }
}
