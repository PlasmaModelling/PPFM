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
    : mixptr(mix), gasptr(gas) { 

        initializeComposition(); 
        initializeFormationEnergies();

        // Just a default PfBox, assign it in the overload and wherever 
        Qbox = new PfBox(mix) ; 

}

Composition::Composition(Mixture* mix, Gas* gas, PfBox* qbox)
    : Composition(mix,gas) {

        Qbox = qbox ; 

}

void Composition::initializeFormationEnergies() {

    epsf.resize(mixptr->getN(), 0.0);
        for (int i = 0; i < mixptr->getN(); i++) 
            epsf[i] = (*mixptr)(i)->formationEnergy();

}

void Composition::initializeComposition() {

    ni.resize(mixptr->getN(), 0.0);
    for (int i = 0; i < mixptr->getN(); i++) {
        
        if ((*mixptr)(i)->getFormula() == "e-")
            ni[i] = 1.e+15;
        else if (mixptr->isBase((*mixptr)(i)))
            ni[i] = gasptr->getPressure() / (KB * gasptr->getTemperature());
        else
            ni[i] = 1.e+10;
    
    }
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

    /* no cut-off adobted, if the Debye Length goes too low  you're maybe 
    trying to compute outside the scope of the model itself. */

    return lambdaD;

}


double Composition::Debye_Rat2002(double T) {

    // ni.back() == ne , in Mixture e- is always the last specie
    return std::sqrt((eps0 * KB * T) / (qe * qe * ni.back()));

}

double Composition::Debye_Ghourui(double T) {

    double sum = 0.0;
    
    for (size_t i = 0; i < mixptr->getN(); i++)
        sum += ( std::pow((*mixptr)(i)->getCharge(), 2.0) * ni[i] ) ; 
    
    sum *= (qe * qe) / (eps0 * KB * T );

    return std::sqrt(1.0 / sum);
    
}

std::vector<double> Composition::totalPartitionFunctions(double T, double P, double lambdaD) {
    
    int N = mixptr->getN();
    std::vector<double> Qtot(N, 0.0);

    Qbox->computePartitionFunctions(T*gasptr->theta->get(), P, lambdaD);

    for (int i = 0; i < N - 1; i++) {

        if (dynamic_cast<BiatomicMolecule*>((*mixptr)(i))) {

            (*Qbox)[i]->computePartitionFunction(T, P, lambdaD);
            Qbox->updateCachedValues();
            
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

    N = mix->getN();
    M = mix->getM();
    L = N-M ;
    C = CompositionMatrix(*mix);
    B = std::vector<std::vector<double>> (M,std::vector<double>(M));
    Bs = std::vector<std::vector<double>> (L,std::vector<double>(M));
    Binv = std::vector<std::vector<double>> (M,std::vector<double>(M));
    v = std::vector<std::vector<double>> (L,std::vector<double>(M));
    b = std::vector<int> (N) ;
    bs = std::vector<int> (L) ;
    A0 = std::vector<double> (M) ;

    /* Being this the ideal solver Peff is always equal to gasptr->getPressure() 
    and can be initialized in the constructor, this is not the case for example for 
    DH corrections */
    Peff = gasptr->getPressure() ; 

}

GodinTrepSahaSolver::GodinTrepSahaSolver(Mixture* mix, Gas* gas, PfBox* qbox)
    : GodinTrepSahaSolver(mix,gas) {

    Qbox = qbox ; 

}

void GodinTrepSahaSolver::restart(){

    *this = GodinTrepSahaSolver(
    
        this->mixptr,
        this->gasptr,
        this->Qbox
    
    ) ;
}

std::vector<double> GodinTrepSahaSolver::Crow(Species* specie,
    const std::map<std::type_index, int>& colmap) {

    std::vector<double> row(colmap.size() + 1, 0.0);

    Element* elem;
    if ((elem = dynamic_cast<Element*>(specie))) {

        if (colmap.find(tipo(elem)) != colmap.end())
            row[colmap.at(tipo(elem))] = 1.0;

    }  else if (PolyAtomicMolecule* poly = dynamic_cast<PolyAtomicMolecule*>(specie)) {

        for (int i = 0; i < poly->numberOfCostituents(); i++) {

            elem = (*poly)[i];
            int part = (*poly)(i);
            if (colmap.find(tipo(elem)) != colmap.end())
                row[colmap.at(tipo(elem))] = (double)part;

        }

    }else if (ChargedSpecies* chrgd = dynamic_cast<ChargedSpecies*>(specie)) {

        elem = chrgd->Constituent();
        if ((elem != nullptr) && (colmap.find(tipo(elem)) != colmap.end()))
            row[colmap.at(tipo(elem))] = 1.0;

    } else {

        throw std::invalid_argument("Unrecognized species type");

    }

    row.back() = specie->getCharge(); // charge last
    return row;

}

std::vector<std::vector<double>> GodinTrepSahaSolver::CompositionMatrix(Mixture& mixx) {

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

    int Nn = b.size();
    int Mm = b.size() - bs.size();
    int* perm = new int[Nn];
    double* n_ord = new double[Nn];

    for (int i = 0; i < Nn; i++) n_ord[i] = ni[i];
    sort(n_ord, perm, Nn);

    int tmp = 1;
    for (int i = 0; i < Mm; i++) {

        int j = 0;
        b[i] = perm[Nn - 1];
        if (C[perm[Nn - 1]][i] == 0) tmp = 0;

        for (int l = 0; l < i; l++)
            if (b[l] == b[i]) tmp = 0;

        while (tmp == 0) {

            j += 1;
            b[i] = perm[Nn - 1 - j];
            if (C[perm[Nn - 1 - j]][i] != 0) tmp = 1;

            for (int l = 0; l < i; l++)
                if (b[l] == b[i]) tmp = 0;

        }
    }

    int l = 0;
    for (int i = 0; i < Nn; i++) {

        tmp = 0;
        for (int j = 0; j < Mm; j++)
            if (i == b[j]) tmp = 1;

        if (tmp == 0) {

            bs[l] = i;
            l += 1;

        }
    }

    delete[] perm;
    delete[] n_ord;
}

void GodinTrepSahaSolver::CompositionSolve( std::optional<double> lambda) {
    
    double theta = gasptr->theta->get() ; 
    double T = gasptr->getTemperature() ; 

    std::vector<Species*> species ; 
    for (int i = 0; i < N; i++)
        species.push_back((*mixptr)(i)) ; 

    baseCalc ( b, bs, C );

    for(int i = 0;i < M; i++) 
        for(int j = 0; j < M; j++) 
            B[i][j] = C[b[i]][j];
        
    for(int i=0;i<L;i++) 
        for(int j=0;j<M;j++) 
            Bs[i][j] = C[bs[i]][j];
        
    lu_inv(Binv,B,M);
    matrix_prod(v,Bs,Binv,L,M);

    std::vector<std::vector<double>> A = ConservationMatrix(*mixptr,*gasptr,C) ;
    
    for( int j = 0 ; j < M-1 ; j++) 
        A0[j] = 0.;
    A0[M-1] = Peff ;

    // Finally, the core of the algorithm described in the reference
    
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
        
        /* Newton base densities */

        /* residuals as in - GODIN: eq 34 */
        residual(R,J,ni,A,A0,v,b,bs,N,M);
        for(int i=0; i<M; i++)
            R[i] = - R[i];
        
        /* dn as in GODIN - eq 33 */
        lu_sistema(dn,J,R,M);

        /*Bdet = lu_det(J,M); */
        for(int i=0; i<N; i++) 
            nold[i] = ni[i];
        
        /*update base densities GODIN eq.30*/
        for(int i = 0; i < M; i++ ) 
            ni[b[i]] += dn[i];
            
        for(int i = 0; i < M; i++ ) {
 
            if (ni[b[i]] < nmin) 
                ni[b[i]] = nmin;
            
        }
        
        // In case of frozen lambda GodinTrepSahaSolver::CompositionSolve 
        // has to be called with lambda value assigned.
        double lambdaEff = lambda.has_value() ? *lambda : getDebyeLength(T);

        // Total partition functions calculation 
        std::vector<double> Qtot = totalPartitionFunctions(T, Peff, lambdaEff);
        
        /* Mass Action Law */
        /* update not-base densities: GODIN - eq 30 */
        for(int j = 0; j < (N-M); j++){
            
            prod = 1.;
            
            for(int i = 0; i < M; i++){

                prod *= pow((ni[b[i]]),v[j][i]) *
                    pow(1./Qtot[b[i]],v[j][i]);
            
            }
            
            ni[bs[j]] = prod * Qtot[bs[j]];
            if(ni[bs[j]] < nmin){ ni[bs[j]] = nmin;}
        
        }

        /* convergence */
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

    // Default Debye model for Debye-Huckle corrections
    setDebyeModel("Ghourui");

    /* Peff changes during computatio, has to be assigned in Solve */

}  

GTSahaDHcorrection::GTSahaDHcorrection(Mixture* mix, Gas* gas, PfBox* qbox)
    : GodinTrepSahaSolver(mix, gas) {
    
    Qbox = qbox ; 

}

void GTSahaDHcorrection::restart() {
    
    *this = GTSahaDHcorrection(
    
        this->mixptr,
        this->gasptr,
        this->Qbox
    
    ) ;

}

double GTSahaDHcorrection::pressureDHcorrected(double T, double lambdaD) {
    
    // Capitelli Fundamental Aspects of Chemical Plasma Physics, 2016,
    // Eq. 6.23-6.24, pp. 105
    // P_DH = -k_B T / (24π λ_D^3)
    return gasptr->getPressure() + ((KB * T) / (24.0 * std::numbers::pi * std::pow(lambdaD, 3.)));

}

std::vector<double> GTSahaDHcorrection::formationEnergyDHcorrected (double lambdaD) {

    /* Be sure corrections are applied only to original 
    formation energies from chemical species */
    initializeFormationEnergies() ; 

    // Capitelli Fundamental Aspects of Chemical Plasma Physics, 2016,
    // Eq. 6.45 pp. 109
    double KDH = (qe * qe) / (8. * std::numbers::pi * eps0 * lambdaD);
    
    // Eq. 6.44 pp. 109 
    std::vector<double> delta_epsf;
    for (size_t i = 0; i < mixptr->getN(); i++)
        delta_epsf.push_back( KDH * std::pow( (*mixptr)(i)->getCharge() + 1., 2.) );
        
    std::vector<double> corrected_epsf;
    for (size_t i = 0; i < mixptr->getN(); i++)
        corrected_epsf.push_back(epsf[i] - delta_epsf[i]);

    return corrected_epsf;

}

void GTSahaDHcorrection::CompositionSolve( std::optional<double> lambda ) {

    double T = gasptr->getTemperature();
    
    // Convergence on Debye Length
    double lambdaOld = getDebyeLength(T);
    double lambdaNew = 0.0;
    double deltaLambda = 1.0;

    // Convergence on Peff
    double PeffOld = gasptr->getPressure();
    double PeffNew = PeffOld;
    double deltaPeff = 1.0;

    int iter = 1;
    // relaxing coefficient (0 < relax <= 1)
    const double relax = 0.3; 

    while ((deltaLambda > tol && deltaPeff > tol) && iter <= maxIter) {

        // Formation Energy Corrections
        epsf = formationEnergyDHcorrected(lambdaOld);

        // Peff = P + PDH
        PeffNew = pressureDHcorrected(T, lambdaOld);  

        // Relaxing on pressure
        double PeffMixed = (1.0 - relax) * PeffOld + relax * PeffNew;
        Peff = PeffMixed;

        // Risolve la composizione con le energie corrette
        GodinTrepSahaSolver::CompositionSolve( lambdaOld );

        // Relaxing on Debye Length
        lambdaNew = getDebyeLength(T);
        double lambdaMixed = (1.0 - relax) * lambdaOld + relax * lambdaNew;

        // Errors
        deltaLambda = std::abs((lambdaMixed - lambdaOld) / lambdaOld);
        deltaPeff   = std::abs((PeffMixed - PeffOld) / PeffOld);

        // Update
        lambdaOld = lambdaMixed;
        PeffOld   = PeffMixed;
        iter += 1;

        // Debug log
        if (iter>=maxIter) {

            std::cout << "[Warning: IterMax reached] T=" << T << " K "
            << "λ_D=" << lambdaMixed << " m, "
            << "Δλ/λ=" << deltaLambda << ", "
            << "Peff=" << Peff << " Pa, "
            << "ΔP/Pid=" << deltaPeff
            << std::endl;
        
        }
       
    }
}
