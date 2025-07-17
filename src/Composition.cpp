#include <stdexcept>
#include <numbers>
#include "Composition.h"
#include "PfBox.h"
#include "PartitionFunction.h"
#include <typeindex>
#include "GasMixture.h"

Composition::Composition ( Mixture* mix , Gas* gas ) {
    
    Qbox = new PfBox( mix ) ;

    int N = mix->getN();
    int M = mix->getM();
    int L = N-M ;
    double P = gas->getPressure() ;
    double T = gas->getTemperature() ; 
    std::vector<Species*> species ; 
    for (int i = 0; i < N; i++)
        species.push_back((*mix)(i)) ; 
    
    ni.resize(N,0.0);
    for (int i = 0; i < N; i++) {
        if (species[i]->getFormula()=="e-") 
            ni[i] = 1e+15 ;
        else if (mix->isBase(species[i]))
            ni[i] = P/(2.*KB*T);        
        else
            ni[i] = 1e+10 ;    
    }
    C = CompositionMatrix(*mix) ; 
} ;

Composition::Composition ( Mixture* mix , Gas* gas , PfBox* qbox) : Composition(mix,gas) { Qbox = qbox ; }

std::vector<double> Composition::Crow ( Species* specie, const std::map<std::type_index , int>& colmap )  {
    
    std::vector<double> row (colmap.size()+1 , 0.) ;

    // dynamic_cast logic to populate CompositionMatrix

    Element* elem ;
    if ( elem = dynamic_cast<Element*>(specie) ) {
    
        if ( colmap.find( tipo(elem) ) != colmap.end() )
            row[colmap.at(tipo(elem))] = 1. ;
    
    } else if ( ChargedSpecies* chrgd = dynamic_cast<ChargedSpecies*>(specie) ) {

        elem = chrgd->Constituent() ; 

        if ( ( elem != nullptr) && (colmap.find(tipo(elem)) != colmap.end() ) ) 
            row[colmap.at(tipo(elem))] = 1. ;
    
    } else if ( PolyAtomicMolecule* poly  = dynamic_cast<PolyAtomicMolecule*> ( specie ) ) {
        
        for (int i = 0; i < poly->numberOfCostituents() ; i++) {
            
            elem = (*poly)[i] ;
            int part = (*poly)(i) ;
            
            if (colmap.find(tipo(elem)) != colmap.end()) 
                row[colmap.at(tipo(elem))] = (double)part ; 
            
        }  
    
    } else {
        throw std::invalid_argument("Unrecognized species type");
    }
    /* Always place the charge last. */
    row.back() = specie->getCharge() ;
    
    return row ;
    
}

std::vector<std::vector<double>> Composition::CompositionMatrix ( Mixture& mixx ) {
    
    int N = mixx.getN() ; 
    int M = mixx.getM() ; 
    
    // build elements-in-the-mixture map for compositionMatrix
    std::map < std::type_index , int > columns ; 
    
    int k = 0 ;
    for (int i = 0; i < N; i++)
        if ( auto elem = dynamic_cast<Element*>(mixx(i)) )
            columns [tipo(elem)] = k++ ;    

    // populate composition matrix by rows
    std::vector<std::vector<double>> CC ( N, std::vector<double> (M) )  ; 

    for (int i = 0; i < N; i++)
        CC[i] = Crow( mixx(i) , columns ) ;
    
    return CC ; 
}
std::vector<std::vector<double>> Composition::ConservationMatrix ( Mixture& mixx, Gas& gass, const std::vector<std::vector<double>>& C ) {
    
    double theta = gass.theta->get() ; 

    int N = mixx.getN() ; 
    int M = mixx.getM() ;

    // build elements-in-the-mixture map for compositionMatrix
    std::map < std::type_index , double > nuclei ; 
    for (int i = 0; i < N; i++)
        if ( auto elem = dynamic_cast<Element*>(mixx(i)) )
            nuclei [tipo(elem)] = 0.0 ;    

    /* Mole fractions of Mixture constituents */
    std::vector<Species*> mfSpecies  ; 
    std::vector<double>   molefractions ; 
    std::tie ( mfSpecies, molefractions ) = mixx.molefractions ;

    /* compute stechiometry unppacking chemical structures */
    for (int i = 0; i < mfSpecies.size(); i++) {

        if (PolyAtomicMolecule* ptr = dynamic_cast<PolyAtomicMolecule*>(mfSpecies[i])) {

            for (size_t j = 0; j < ptr->numberOfCostituents(); j++) {

                auto elem = (*ptr)[j] ; 
                double coeff = (*ptr)(j) ;
                nuclei[tipo(elem)] = molefractions[i]*coeff ;

            }

        } else if (Element* ptr = dynamic_cast<Element*>( mfSpecies[i] ) ) {

            nuclei[tipo(ptr)] = molefractions[i] ;

        } else {

            throw std::invalid_argument("Invalid Specie set as molefraction. ") ; 

        }
    }

    std::vector<double> stechiometry ;
    for ( const auto& pair : nuclei )
        stechiometry.push_back(pair.second) ; 
    
    // calculation of conservation Matri2
    double pressure = gass.getPressure() ; 
    double temperature = gass.getTemperature() ; 
    
    std::vector<std::vector<double>> AA(M,std::vector<double>(N)) ; 
    
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (i < M-2) {
    
                AA[i][j] = stechiometry[i+1]*C[j][i] - stechiometry[i]*C[j][i+1] ; 
    
            } else if ( i == M-2 ) {

                //electrical neutrality
                AA[i][j] = C[j][M-1] ;

            } else {
                
                AA[i][j] = KB * temperature ; 
                if ( j == N-1 )
                    AA[i][j] = KB * temperature * theta ; //THETA QUI
            }
        }        
    }
    
    return AA ; 

}


void Composition::compositionSolve(Mixture* mix , Gas* gas) {
    
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
        
        /* 
            Va reso iterativo questo calcolo delle funzioni di partizione per tenere conto del 
            Lowering Ionization Potential, tuttavia, le iterazioni avverranno solo quando il 
            metodo scelto sarà l'ab-initio in quanto solo i criteri di cut-off per le configurazioni 
            elettroniche sarà dipendente dalla lunghezza di Debye.
        */
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

double Composition::getDebyeLength ( double temperatura ) {
    
    // [debye length] = m
   
    // formulazione Ghouroi
    /* double sum = 0. ; 
    for (size_t i = 0; i < ni.size(); i++)     
        sum += pow(C[i][1] , 2.) * ni[i] / temperatura ;
    
    sum *= (qe*qe) / (eps0*KB) ;
    return sqrt(1./sum);
    */

    // formulazione Rat 2002 se la temperatura passata è Th
    // source eq.45 UNIBO 2008 
    return sqrt((eps0*KB*temperatura)/(qe*qe*ni.back()));

}

double Composition::operator()(int i) {
    return ni[i];
} ;

std::vector<double> Composition::compositions() {
    
    std::vector<double> ns ;
    for (int i = 0; i < ni.size() ; i++) {
        ns.push_back( (*this)(i) ) ;
    }
    return ns ; 
}

std::vector<double> Composition::compositions(double conversion) {
    
    std::vector<double> ns ;
    for (int i = 0; i < ni.size() ; i++) {
        ns.push_back( (*this)(i) * conversion ) ;
    }
    return ns ; 
}

void Composition::baseCalc(std::vector<int>& b, std::vector<int>& bs, const std::vector<std::vector<double>>& C){

    int i, j, l, tmp;
    int N = b.size() ;
    int M = b.size() - bs.size() ;
    int *perm = new int[N];
    double *n_ord = new double[N];

    for (i = 0; i < N; i++)
        n_ord[i] = ni[i];
    sort(n_ord, perm, N);

    tmp = 1;
    for (i = 0; i < M; i++)
    {
        j = 0;
        b[i] = perm[N - 1];
        if (C[perm[N - 1]][i] == 0)
            tmp = 0;
        for (l = 0; l < i; l++)
        {
            if (b[l] == b[i])
                tmp = 0;
        }
        while (tmp == 0)
        {
            j += 1;
            b[i] = perm[N - 1 - j];
            if (C[perm[N - 1 - j]][i] != 0)
                tmp = 1;
            for (l = 0; l < i; l++)
            {
                if (b[l] == b[i])
                    tmp = 0;
            }
        }
    }

    l = 0;
    for (i = 0; i < N; i++)
    {
        tmp = 0;
        for (j = 0; j < M; j++)
            if (i == b[j])
                tmp = 1;
        if (tmp == 0)
        {
            bs[l] = i;
            l += 1;
        }
    }
}


std::string CompositionCsv::BuildFileName(const std::string& filename ) const {

    return "Comp_" + filename + ".csv";  

}

void CompositionCsv::PrepareHeader() {
        
    header = "Te [K],";

    for (int i = 0; i < ni.size(); i++) 
        header += "n_{" + (*mixptr)(i)->getFormula() + "}[#/m^3],";

    header += "n_{tot}[#/m^3]";

}

void CompositionCsv::PrepareData ( const std::vector<double>& temperatureRange, GasMixture* mix ) {
        
    double T0 = mix->getTemperature(); 
    double theta = mix->theta->get();

    std::vector<std::vector<double>> ns(temperatureRange.size(), std::vector<double>(mix->getN()));         
    data.resize(temperatureRange.size(), std::vector<double>(ns[0].size()+2)) ;

    for (int i = 0; i < temperatureRange.size(); i++) {
        
        mix->setT(temperatureRange[i]) ; 
        compositionSolve( mix, mix ) ; 
        
        ns[i] = this->compositions() ;
        
        data[i][0] = temperatureRange[i] * theta ; 
        
        size_t j = 1; 
        for ( j; j < ns[0].size()+1; j++)
            data[i][j] = ns[i][j-1] ;        
        
        data[i][j] = ntot() ;

    }
    
    // SETBACK  
    mix->setT(T0) ; 
    mix->restartComposition() ; 

}

void CompositionCsv::PrintMessage(const std::string& filename) {

    std::cout << "Composition " << filename << " printed." << std::endl ;

}

CompositionCsv::CompositionCsv ( Mixture* mix , Gas* gas ) : Composition(mix,gas) { 

    mixptr = mix; 

}
CompositionCsv::CompositionCsv ( Mixture* mix , Gas* gas , PfBox* qbox ) : Composition(mix,gas,qbox) { 

    mixptr = mix; 

}
CompositionCsv::CompositionCsv ( Mixture* mix , Gas* gas , const std::string& folder ) : CompositionCsv(mix,gas) { 

    customFolder = folder; 

}
CompositionCsv::CompositionCsv ( Mixture* mix , Gas* gas , PfBox* qbox, const std::string& folder ) : CompositionCsv(mix,gas,qbox) { 

    customFolder = folder; 

}