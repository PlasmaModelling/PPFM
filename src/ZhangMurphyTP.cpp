#include "ZhangMurphyTP.h"
#include "GasMixture.h"
#include <chrono>

void ZhangMurphyTP::computeTransport( GasMixture* gasmix ) {
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Computes collision Integrals
    Transport::QtCalc(gasmix) ;

    // init for ZMcoefficients (improves performances)
    ZMCoefficients::init(gasmix) ;

    D.resize(N, std::vector<double>(N));
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            D[i][j] = Dij(gasmix, 3, i, j);
        }
    }

    DT.resize(N, 0.);

    for (int i = 0; i < N; i++) {
        DT[i] = DiT(gasmix, 3, i);
    }

    Tp.resize(6, 0.);

    // Calculate properties in parallel or sequentially
    Tp[0] = ThermalCondEl(gasmix, 3);

    Tp[1] = ThermalCondHeavy(gasmix, 2);

    Tp[2] = Viscosity(gasmix, 1);

    Tp[3] = ElCond(gasmix, 4);

    Tp[4] = NeThermalCondEl(gasmix, 3);

    Tp[5] = NeThermalCondHeavy(gasmix, 2);

    DiffTheta.resize(N, std::vector<double>(N));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            DiffTheta[i][j] = DijTheta(gasmix, 3, i, j);
        }
    }

    Dtheta.resize(N, 0.);

    for (int i = 0; i < N; i++) {
        Dtheta[i] = DiTheta(gasmix, 3, i);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    // std::cout << "Zhang T = "<<gasmix->getTemperature()<<" elapsed in: "<< duration.count() << std::endl;

}

double ZhangMurphyTP::DijTheta ( GasMixture* gasmix , int order , int ii , int jj ) {
    
    double Ti = T ;
    if ( ii == N-1 )    
            Ti= Te ;

    double Ei0 = ei0(gasmix,order,jj,ii,ii) ;

    double D = ((n[ii]*rho*kB*Ti)/(ntot*mass[jj]))*sqrt((kB*Ti)/(2.*mass[ii]))*Ei0 ;    
    
    return D * 1.e-12 ;

}

double ZhangMurphyTP::DiTheta  ( GasMixture* gasmix , int order , int ii ) {
    
    double Ti = T ;
    if (ii == N-1)    
        Ti= Te ;
    
    double sum = 0. ; 
    for (int i = 0; i < N; i++)
        sum += mass[i] * DijTheta(gasmix,order,ii,i) * 1.e+12 * wi(gasmix,i) ; 
    
    sum *= (ntot/(n[ii]*rho*kB*Ti)) ;
    
    double Ditheta = DiThetaStar(gasmix,order,ii) / (n[ii]*mass[ii]) ;

    return ( Ditheta - sum ) * 1.e+3 ; 
    
}


double ZhangMurphyTP::NeThermalCondEl    ( GasMixture* gasmix, int order ) {
        
    double Ti = Te ;

    std::vector<std::vector<double>> Djk(N,std::vector<double>(N)) ;
    // #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)        
            Djk[i][j] = Dij(gasmix,order,i,j) * 1.e+12 * mass[j] ;

    std::vector<std::vector<double>> Ejk(N,std::vector<double>(N)) ;
    lu_inv(Ejk,Djk,N) ;

    std::vector <double> DT_theta (N,0.0) ;
    // #pragma omp parallel for
    for (int k = 0; k < N; k++) 
        DT_theta[k] = DiTheta(gasmix,order,k) ;
    
    std::vector <double> lambdaijDD ( N, 0.0 ) ;
    // #pragma omp parallel for
    for (int k = 0; k < N; k++) 
        lambdaijDD[k] = lambdaijD(gasmix,order,N-1,k) ;
    

    double sum1 = 0. ; double sum2 = 0. ; 
    // #pragma omp parallel for reduction(+:sum2)
    for (int j = 0; j < N; j++) {
        sum1 = 0. ;
        // #pragma omp parallel for reduction(+:sum1)
        for (int k = 0; k < N; k++) {
        
            double Tk = T ;
            if (k == N-1) 
                Tk = Te ;
            
            sum1 += 1.e-3 * DT_theta[k]*
                (lambdaijDD[j] * n[k] * Tk * Ejk[j][k]) / (n[j]*mass[j]*Ti) ;
        }
        sum2 += sum1 ; 
    }

    sum2 *= (rho*kB)/(ntot) ;
    double litheta = lambdaiPrimeTheta(gasmix,order,N-1) ;
    
    return (litheta + sum2) * 1.e-9 ;

}

double ZhangMurphyTP::NeThermalCondHeavy ( GasMixture* gasmix, int order ) {
        
    double Ti = T;
 
    std::vector<std::vector<double>> Djk(N, std::vector<double>(N));
    // #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Djk[i][j] = Dij(gasmix, order, i, j) * 1.e+12 * mass[j];
        }
    }

    std::vector<std::vector<double>> Ejk(N, std::vector<double>(N));
    lu_inv(Ejk, Djk, N);

    std::vector<double> lipp( N-1, 0.0 ) ;
    // #pragma omp parallel for
    for (int i = 0; i < N-1; i++) {
        lipp[i]= (lambdaiPrimeTheta(gasmix, order, i)) ;
    }

    // Precomputazione di lambdaijDD
    std::vector<std::vector<double>> lambdaijDD(N, std::vector<double>(N, 0.0));
    // #pragma omp parallel for collapse(2)
    for (int h = 0; h < N-1; h++) {
        for (int j = 0; j < N; j++) {
            lambdaijDD[h][j] = lambdaijD(gasmix, order, h, j);
        }
    }

    std::vector<double> DTtheta (N,0.0) ;
    // #pragma omp parallel for
    for (int k = 0; k < N; k++)
        DTtheta[k] = DiTheta(gasmix, order, k) ;
    
    

    double lambdaH = 0.0;
    
    // #pragma omp parallel for reduction(+:lambdaH)
    for (int h = 0; h < N-1; h++) {
        double sum2 = 0.0;
        // #pragma omp parallel for reduction(+:sum2)
        for (int j = 0; j < N; j++) {
            double sum1 = 0.0;
            // #pragma omp parallel for reduction(+:sum1)
            for (int k = 0; k < N; k++) {
                double Tk = T;
                if (k == N-1)
                    Tk = Te;

                sum1 += 1.e-3 * DTtheta[k] *
                    (lambdaijDD[h][j] * n[k] * Tk * Ejk[j][k]) / (n[j] * mass[j] * Ti);
            }

            sum2 += sum1;
        }

        sum2 *= (rho * kB) / ntot;
        lambdaH += lipp[h] + sum2;
    }

    return lambdaH * 1.e-9;
}

double ZhangMurphyTP::Viscosity ( GasMixture* gasmix, int order ) {
    
    double Th = T  ;
        
    int eps ;
    if (order<=2)    
        eps = order;
    else 
        eps = 1;

    double mue = 0.5*kB * ne * Te * b10 ( gasmix ) ; 

    std::vector<std::vector<double>> QQ (
        (eps*(N-1))+1 , std::vector<double>((eps*(N-1))+1)
    ) ;
    std::vector<std::vector<double>> Q (
        (eps*(N-1)) , std::vector<double>((eps*(N-1)))
    ) ;
    for (int m = 0; m < eps; m++) {
        for (int p = 0; p < eps; p++) {
            for (int i = 0; i < N-1; i++) {
                for (int j = 0; j < N-1; j++) {
                    Q[(m*(N-1))+i][(p*(N-1))+j] = 
                        qcapmpij(gasmix,m,p,i,j) ; 
                    QQ[(m*(N-1))+i][(p*(N-1))+j] = 
                        Q[(m*(N-1))+i][(p*(N-1))+j] ;
                }
            }
        }        
    }
    
    for (int p = 0; p < eps; p++) {
        for (int j = 0; j < N-1; j++) {
            if (p==0)            
                QQ[eps*(N-1)][(p*(N-1))+j] = n[j] ; 
            else
                QQ[eps*(N-1)][(p*(N-1))+j] = 0. ; 
        }
    }
    for (int m = 0; m < eps; m++) {
        for (int i = 0; i < N-1; i++) {
            QQ[(m*(N-1))+i][eps*(N-1)] = 5.*n[i] - qmpi1V(gasmix,m,0,i)*b10(gasmix) - 
                qmpi1V(gasmix,m,1,i)*b11(gasmix) ;
        }
    }
    QQ[eps*(N-1)][eps*(N-1)] = 0. ; 

    double detQQ = lu_det(QQ,QQ.size()) ; 
    double detQ  = lu_det(Q , Q.size()) ; 

    double muh = -0.5*kB*Th*(detQQ/detQ) ;

    return (mue + muh) * 1.e+3 ; 

}

double ZhangMurphyTP::ElCond ( GasMixture* gasmix, int order ) { 

    std::vector<int> Z ;
    for (int i = 0; i < N; i++) 
        Z.push_back((*gasmix)(i)->getCharge()) ;
      
    std::vector<std::vector<double>> DD ( N, std::vector<double>(N) ) ; 
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            DD[i][j] = Dij(gasmix,order,i,j) ;   

    double sigma = 0. ;
    // #pragma omp parallel for reduction (+: sigma)
    for (int i = 0; i < N; i++) {
        
        double Ti = T ;
        if (i == N-1)
            Ti = Te ;
        double sum = 0. ;
        // #pragma omp parallel for reduction (+: sum)
        for (int j = 0; j < N; j++) {
            
            if (std::isnan(DD[i][j])||std::isinf(DD[i][j]))
                continue;
                    
            sum += (mass[j]*n[j]*Z[j]*DD[i][j]) ;
        }
                
        sum *= (Z[i]/Ti) ;
        sigma += sum ; 
    }
    
    sigma *= (-pow(q_e,2.)*ntot)/(rho*kB) ; 

    return sigma * 1.e+33 ;    

} 

double ZhangMurphyTP::Dij (GasMixture* gasmix , int order , int i , int j ) { 
    
    double ni = n[i] ; 
    double mj = mass[j] ;
    double mi = mass[i] ;
    double Ti = T;

    if (i == N-1)
        Ti = Te ; 
    
    double CiZero = ci0(gasmix,order,j,i,i) ; 

    return 1.e-12 * ((ni*rho*kB*Ti)/(ntot*mj)) * sqrt((kB*Ti)/(2.*mi)) * CiZero ;    

}

double ZhangMurphyTP::DiT ( GasMixture* gasmix , int order , int i ) { 
    
    double Ti = T;

    if (dynamic_cast<Electron*>((*gasmix)(i)))
        Ti = Te ; 

    double ni = n[i] ;
    double mi = mass[i] ; 

    double aiZero = ai0(gasmix,order,i) ;
    
    return ni*mi*sqrt((kB*Ti)/(2.*mi))*aiZero * 1.e+3 ; 

}

double ZhangMurphyTP::ThermalCondEl ( GasMixture* gasmix, int order ) {

    double Th = T ;

    std::vector<double> lijD(N,0.);
    // #pragma omp parallel for 
    for (int j = 0; j < N; j++) {
        lijD[j] = lambdaijD(gasmix,order,N-1,j) ;
    }
    
    std::vector<double> DkT(N,0.);
    // #pragma omp parallel for 
    for (int j = 0; j < N; j++) 
        DkT[j] = DiT(gasmix, order, j) * 1.e-3;
    

    std::vector<std::vector<double>> Dmjk ( N,
        std::vector<double>(N) ) ;
    // #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++)         
            Dmjk[i][j] = Dij(gasmix,order,i,j) * 1.e+12 * mass[j] ;
    
    
    std::vector<std::vector<double>> Ejk ( N,
        std::vector<double>(N) ) ;
    lu_inv(Ejk,Dmjk,N) ;

    double Tk ; 
    double sum1 = 0. ;
    double sum2 = 0. ;
    // #pragma omp parallel for reduction(+:sum2)
    for (int j = 0; j < N; j++) {
        sum1 = 0. ;
        // #pragma omp parallel for reduction(+:sum1)
        for (int k = 0; k < N; k++) {
            
            Tk = T ; 
            if (k == N-1)
                Tk = Te ;
    
            sum1 += (lijD[j]*Tk*Ejk[j][k]*DkT[k])/(n[j]*mass[j]*mass[k]*Th) ;
    
        }    
        sum2 += sum1 ;
    }
    
    double liP = lambdaiPrime(gasmix,order,N-1) ;

    return (liP + (((rho*kB)/ntot)*sum2) ) * 1.e-13 * 1.e+4 ;
    
}

double ZhangMurphyTP::ThermalCondHeavy ( GasMixture* gasmix, int order ) {

    double Th = T ;

    std::vector<std::vector<double>> lijD (N,
        std::vector<double>(N)) ;
    std::vector<std::vector<double>> Dmjk ( N,
        std::vector<double>(N)) ;

    // #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            lijD[i][j] = lambdaijD(gasmix,order,i,j) ;
            Dmjk[i][j] = Dij(gasmix,order,i,j) * 1.e+12 * mass[j] ;
        }
    }

    std::vector<double> DkT(N,0.);
    // #pragma omp parallel for 
    for (int j = 0; j < N; j++) 
        DkT[j] = DiT(gasmix, order, j) * 1.e-3;

    std::vector<std::vector<double>> Ejk ( N,
        std::vector<double>(N) ) ;
    lu_inv(Ejk,Dmjk,N) ;

    std::vector<double> liP ( N-1 ) ;
    //outer cicle on heavy Species
    // #pragma omp parallel for
    for (int ii = 0; ii < N-1; ii++)
        liP[ii] = lambdaiPrime(gasmix,order,ii) ;
    
    double lambdaH = 0.;
    double Tk ;
    //outer cicle on heavy Species
    // #pragma omp parallel for reduction ( +:lambdaH )
    for (int ii = 0; ii < N-1; ii++) {
        double sum2 = 0. ;
        // #pragma omp parallel for reduction ( +:sum2 )
        for (int j = 0; j < N; j++) {
            double sum1 = 0. ;
            // #pragma omp parallel for reduction ( +:sum1 )
            for (int k = 0; k < N; k++) {
                
                Tk = T ; 
                if ( k == N-1 )
                    Tk = Te ;

                sum1 += (lijD[ii][j]*Tk*Ejk[j][k]*DkT[k])/(n[j]*mass[j]*mass[k]*Th) ;

            }
            sum2 += sum1 ;
        }
        lambdaH += liP[ii] + (( ( rho*kB ) / ntot ) * sum2 ) ;
    }
    
    return lambdaH * 1.e-13 * 1.e+4 ; 

}


std::string ZhangTpCsv::BuildFileName(const std::string& filename ) const  {

    return "TP_" + filename + ".csv";  

}

void ZhangTpCsv::PrepareHeader() {

    header = "Th [K], kₑ [W/(m·K)], kₕ [W/(m·K)], μ [Pa·s], σ [S/m], kₑθ [W/m], kₕθ [W/m]";

}

void ZhangTpCsv::PrintMessage(const std::string& filename)  { 

    std::cout << "Zhang, Murphy et.al. Transport Properties " << filename << " printed." << std::endl ;

}

ZhangTpCsv::ZhangTpCsv(CiBox* cbx) : ZhangMurphyTP(cbx) {};

ZhangTpCsv::ZhangTpCsv(CiBox* cbx, const std::string& folder ) : ZhangMurphyTP(cbx) { customFolder = folder; };

void ZhangTpCsv::PrepareData ( const std::vector<double>& temperatureRange, GasMixture* gasmix) {
    
    double T0 = temperatureRange[0] ;
    double theta = gasmix->theta->get() ;

    data.resize(temperatureRange.size(),std::vector<double>(7));

    for (int i = 0; i < temperatureRange.size(); i++)
    {
        gasmix->setT(temperatureRange[i]) ; 
        computeTransport(gasmix) ; 

        // Te
        data[i][0] = temperatureRange[i] ; 
        data[i][1] = Tp[0] ;
        data[i][2] = Tp[1] ;
        data[i][3] = Tp[2] ;
        data[i][4] = Tp[3] ;
        data[i][5] = Tp[4] ;
        data[i][6] = Tp[5] ;

    }
    
    // SETBACK
    gasmix->setT(T0);
    gasmix->restartComposition();
    computeTransport(gasmix);

}