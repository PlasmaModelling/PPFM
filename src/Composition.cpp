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

        Qbox = new PfBox(mix);
        ni.resize(mix->getN(), 0.0);

}

Composition::Composition(Mixture* mix, Gas* gas, PfBox* qbox)
    : mixptr(mix), gasptr(gas), Qbox(qbox) {

        ni.resize(mix->getN(), 0.0);

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

    switch (debyeChoice) {
        
        case DebyeModel::Rat2002Th:

            return Debye_Rat2002(T);
        
        case DebyeModel::Rat2002Te:
        
            return Debye_Rat2002( T * gasptr->theta->get());
        
        case DebyeModel::Ghourui:
        
            return Debye_Ghourui(T);
        
        default:
        
            throw std::runtime_error("Unknown Debye model selected");

    }
}

double Composition::Debye_Rat2002(double T) {

    return std::sqrt((eps0 * KB * T) / (qe * qe * ni.back()));

}

double Composition::Debye_Ghourui(double T) {

    double sum = 0.0;
    int M = mixptr->getM();
    
    for (size_t i = 0; i < ni.size(); i++)
        sum += ( std::pow((*mixptr)(i)->getCharge(), 2.0) * ni[i] ) / T; 
    
    sum *= (qe * qe) / (eps0 * KB);

    return std::sqrt(1.0 / sum);
}

// __________________________ GodinTrepSahaSolver __________________________ //

GodinTrepSahaSolver::GodinTrepSahaSolver(Mixture* mix, Gas* gas)
    : Composition(mix, gas) {

    int N = mix->getN();
    int M = mix->getM();
    double P = gas->getPressure();
    double T = gas->getTemperature();

    std::vector<Species*> species;
    for (int i = 0; i < N; i++)
        species.push_back((*mix)(i));

    for (int i = 0; i < N; i++) {
        if (species[i]->getFormula() == "e-")
            ni[i] = 1.e15;
        else if (mix->isBase(species[i]))
            ni[i] = P / (2.0 * KB * T);
        else
            ni[i] = 1.e10;
    }

    C = CompositionMatrix(*mix);
}

GodinTrepSahaSolver::GodinTrepSahaSolver(Mixture* mix, Gas* gas, PfBox* qbox)
    : Composition(mix, gas, qbox) {

    C = CompositionMatrix(*mix);

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
    double P = gas->getPressure() ;
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
    A0[M-1] = P ;

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
        for(int i = 0; i < M; i++ ) {
            ni[b[i]] += dn[i];
            }
        for(int i = 0; i < M; i++ ) {
            if (ni[b[i]] < nmin) { 
                ni[b[i]] = nmin;
                }
        }
        
        // calcolo funzioni di partizione totali
        std::vector<double> Qtot (N,0.0) ;
        
        Qbox->computePartitionFunctions( T*theta , P , getDebyeLength(T) ) ; 
        
        for (int i = 0; i < N-1; i++) {
            
            if ( dynamic_cast<BiatomicMolecule*>(species[i])) {
                
                (*Qbox)[i]->computePartitionFunction(T, P, getDebyeLength(T)) ;

                Qtot[i] = pow ((( 2. * std::numbers::pi * species[i]->getMass() * KB *
                T ) / ( hPlanck * hPlanck )) , 1.5) * std::exp ( - (( 
                    species[i]->formationEnergy() ) / (KB*T))) 
                        * (*Qbox)(i) ;
        
            } else {
        
                Qtot[i] = pow ((( 2. * std::numbers::pi * species[i]->getMass() * KB *
                    T ) / ( hPlanck * hPlanck )) , 1.5) * std::exp ( - (( 
                        species[i]->formationEnergy() ) / (KB*T*theta)) ) 
                            * (*Qbox)(i) ;
        
            }
        }
        
        Qtot[N-1] = pow ((( 2 * std::numbers::pi * species[N-1]->getMass() * KB * T * theta ) / 
                (hPlanck*hPlanck)) , 1.5 ) * (*Qbox)(N-1) ;
        
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
// __________________________ CompositionCsv __________________________ //

CompositionCsv::CompositionCsv(Composition* solver) : solver(solver) {}

CompositionCsv::CompositionCsv(Composition* solver, const std::string& folder)
    : solver(solver) {
        customFolder = folder;
}

std::string CompositionCsv::BuildFileName(const std::string& filename) const {
    return "Comp_" + filename + ".csv";
}

void CompositionCsv::PrepareHeader() {

    header = "Te [K],";
    for (int i = 0; i < solver->compositions().size(); i++) {
        header += "n_{" + (*solver->mixptr)(i)->getFormula() + "}[#/m^3],";
    }
    header += "n_{tot}[#/m^3]";
}

void CompositionCsv::PrepareData(const std::vector<double>& temperatureRange, GasMixture* mix) {

    double T0 = mix->getTemperature();
    double theta = mix->theta->get();

    std::vector<std::vector<double>> ns(temperatureRange.size(),
        std::vector<double>(mix->getN()));

    data.resize(temperatureRange.size(), std::vector<double>(mix->getN() + 2));

    for (int i = 0; i < temperatureRange.size(); i++) {
        
        mix->setT(temperatureRange[i]);
        solver->CompositionSolve(mix, mix);

        ns[i] = solver->compositions();

        data[i][0] = temperatureRange[i] * theta;
        for (int j = 1; j < ns[0].size() + 1; j++)
            data[i][j] = ns[i][j - 1];
        
        data[i][ns[0].size() + 1] = solver->ntot();
    }

    mix->setT(T0);
    mix->restartComposition();
}

void CompositionCsv::PrintMessage(const std::string& filename) {
    std::cout << "Composition " << filename << " printed." << std::endl;
}
