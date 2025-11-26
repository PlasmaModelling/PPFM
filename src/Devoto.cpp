// PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "Devoto.h"
#include <chrono>
#include "GasMixture.h"
#include "Thermodynamics.h" // for reactive and internal contributions to thermal conductivities
#include <numbers>

void DevotoTP::computeTransport ( GasMixture* gasmix ) {
    
    auto start = std::chrono::high_resolution_clock::now();

    QtCalc(gasmix) ;

    D.resize( gasmix->getN() , std::vector<double>(gasmix->getN() )) ;
    #pragma omp parallel for collapse(2)
    for ( int i = 0; i < gasmix->getN(); i++ ) 
        for ( int j = 0; j < gasmix->getN(); j++ ) 
            D[i][j] = Dij ( gasmix, orders[0], i, j )  ;
    
    DT.resize( gasmix->getN() ,0.0 ) ; 
    #pragma omp parallel for
    for ( int i = 0 ; i < gasmix->getN() ; i++ )
        DT[i] = DiT ( gasmix , orders[1], i ) ; 
    
    Tp.resize(5,0.0) ; 
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            Tp[0] = ThermalCondEl(gasmix, orders[2]);
            #pragma omp section
            Tp[1] = ThermalCondHeavy(gasmix, orders[3]);
            #pragma omp section
            Tp[2] = Viscosity(gasmix, orders[4]);
            #pragma omp section
            Tp[3] = ElCond(gasmix, orders[5]);
            #pragma omp section
            Tp[4] = Qeh(gasmix);
        }
    }

    auto end = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> duration = end - start;
    // std::cout << "Devoto T = "<<gasmix->getTemperature()<<" elapsed in: "<< duration.count() << std::endl;
    
} ;

double DevotoTP::ThermalCondEl ( GasMixture* gasmix, int order ) {
    
    /* [#/micron^3] */
    std::vector<double> n = gasmix->getCompositionObj()->compositions(1.e-18) ; 
    /* [g] */
    std::vector<double> mass = gasmix->masses(1.e+3) ; 
    int N_SPC = gasmix->getN () ;
    double Te = gasmix->getTemperature() * gasmix->theta->get() ;

    int i, j ;
    double tcond ;
    std::vector<std::vector<double>> qq (order,std::vector<double>(order)) ;
    std::vector<std::vector<double>> QQ (order-1,std::vector<double>(order-1)) ;

    for (i = 1; i < order; i++) {
        for (j = 1; j < order; j++) {
            qq[i-1][j-1] = qsimpmpij( gasmix , i , j ) ;
        }
    }
    for (i = 2; i < order; i++) {
        for (j = 2; j < order; j++) {
            QQ[i-2][j-2] = qsimpmpij(  gasmix , i , j ) ;
        }
    }

    double detQ, detq ;
    
    switch (order) {
    case 2 :
        detQ = 1. ;  
        detq = qq[0][0] ;
        break;
    case 1 :
        detQ = 1;
        detq = 1;
    default:
        detQ = lu_det(QQ,order-2);
        detq = lu_det(qq,order-1);
        break;
    }

    tcond = (75. / 8.) * pow(n[N_SPC - 1], 2) * kB * sqrt(2. * PI * kB * Te / mass[N_SPC - 1]) * 
         /* (1. / (qq[0][0] - (pow(qq[0][1], 2) / qq[1][1]))) */ (detQ/detq) ;

    return tcond*1e-13*1e4 ;
}
    
double DevotoTP::ThermalCondHeavy ( GasMixture* gasmix, int order ) {

    // vecchio codice masse comuni per ioni di...
    // nuovo codice masse in Species da NIST con elettroni persi.
 
    int N_SPC = gasmix->getN() ;
    std::vector<double> n = gasmix->getCompositionObj()->compositions(1e-18) ; 
    std::vector<double> mass = gasmix->masses(1e+3) ;
    double T = gasmix->getTemperature() ;

    int i, j, m, p;
    double tcond;
    double det_Q, det_q;
    int N_ORDS = order ;

    std::vector<std::vector<double>> Q (N_ORDS * N_SPC + 1,
        std::vector<double>(N_ORDS * N_SPC + 1)) ;
    
    std::vector<std::vector<double>> q (N_ORDS * N_SPC, 
        std::vector<double>(N_ORDS * N_SPC)) ;

    for (m = 0; m < N_ORDS; m++)
    {
        for (p = 0; p < N_ORDS; p++)
        {
            for (i = 0; i < N_SPC; i++)
            {
                for (j = 0; j < N_SPC; j++)
                {
                    q[m * N_SPC + i][p * N_SPC + j] = qmpij( gasmix, m, p, i, j);
                }
            }
        }
    }

    for (m = 0; m < N_ORDS; m++)
    {
        for (p = 0; p < N_ORDS; p++)
        {
            for (i = 0; i < N_SPC; i++)
            {
                for (j = 0; j < N_SPC; j++)
                {
                    Q[m * N_SPC + i][p * N_SPC + j] = q[m * N_SPC + i][p * N_SPC + j];
                }
            }
        }
    }
    m = 1;
    for (i = 0; i < (N_SPC - 1); i++)
    { // electron last
        Q[m * N_SPC + i][N_ORDS * N_SPC] = n[i];
    }
    p = 1;
    for (j = 0; j < (N_SPC - 1); j++)
    { // electron last
        Q[N_ORDS * N_SPC][p * N_SPC + j] = n[j] / sqrt(mass[j]);
    }
    
    det_Q = DetLU(Q, N_ORDS * N_SPC + 1);
    det_q = DetLU(q, N_ORDS * N_SPC);

    tcond = -(75. / 8.) * kB * sqrt(2. * PI * kB * T) * (det_Q / det_q);

    return tcond * 1e-13 * 1e4 ;
}

double DevotoTP::Viscosity ( GasMixture* gasmix, int order ) {

    int N_SPC = gasmix->getN() ;
    std::vector<double> n = gasmix->getCompositionObj()->compositions(1e-18) ; 
    std::vector<double> mass = gasmix->masses(1e+3) ;
    double T = gasmix->getTemperature() ;
    
    double visc;
    int i, j, m, p;
    double det_qcap;
    int M_ORDV = order ;
    int N_Q = M_ORDV * N_SPC;

    std::vector<std::vector<double>> q_diff (N_Q+1,std::vector<double>(N_Q+1)) ;
    std::vector<std::vector<double>> qcap (N_Q,std::vector<double>(N_Q)) ;
    
    for (m = 0; m < M_ORDV; m++)
        for (p = 0; p < M_ORDV; p++)
            for (i = 0; i < N_SPC; i++)
                for (j = 0; j < N_SPC; j++)
                    qcap[m * N_SPC + i][p * N_SPC + j] = qcapmpij( gasmix, m, p, i, j);

    det_qcap = DetLU(qcap, N_Q);

    for (i = 0; i < N_Q; i++)
        for (j = 0; j < N_Q; j++)
            q_diff[i][j] = qcap[i][j];

    m = 0;

    for (i = 0; i < N_SPC; i++)
        q_diff[m * N_SPC + i][N_Q] = n[i] * sqrt(mass[i]);

    p = 0;

    for (j = 0; j < N_SPC; j++)
        q_diff[N_Q][p * N_SPC + j] = n[j];

    visc = -((5. / 2.) * sqrt(2. * PI * kB * T) * DetLU(q_diff, N_Q + 1) / det_qcap);

    return visc * 1e7 * 1e-4 ;
}

double DevotoTP::ElCond( GasMixture* gasmix, int order ) {

    int N_SPC = gasmix->getN() ;
    std::vector<double> n = gasmix->getCompositionObj()->compositions(1e-18) ; 
    std::vector<double> mass = gasmix->masses(1e+3) ;
    double T = gasmix->getTemperature() ;

    int i, j;
    double sigma = 0;
    double Tex;
    double det_Q, det_q;
    int N_ORDS = order ; 

    std::vector<std::vector<double>> QQ (order-1,std::vector<double>(order-1)) ;
    std::vector<std::vector<double>> qq (order,std::vector<double>(order)) ;
    
    Tex = T * gasmix->theta->get() ;

    for (i = 0; i < N_ORDS - 1; i++)
    {
        for (j = 0; j < N_ORDS - 1; j++)
        {
            QQ[i][j] = qsimpmpij( gasmix, i + 1, j + 1);
        }
    }
    det_Q = DetLU(QQ, N_ORDS - 1);

    for (i = 0; i < N_ORDS; i++)
    {
        for (j = 0; j < N_ORDS; j++)
        {
            qq[i][j] = qsimpmpij( gasmix, i, j);
        }
    }
    det_q = DetLU(qq, N_ORDS);

    sigma = 3. * q_e * q_e * pow(n[N_SPC - 1], 2) * sqrt(std::numbers::pi / (2. * kB * Tex * mass[N_SPC - 1])) * (det_Q / det_q);

    return sigma * 1e33*1e-12 ;
        
}

double DevotoTP::Qeh( GasMixture* gasmix ) {
    
    int N_SPC = gasmix->getN() ;    
    std::vector<double> n = gasmix->getCompositionObj()->compositions(1e-18) ; 
    std::vector<double> mass = gasmix->masses(1e+3) ;
    double T = gasmix->getTemperature() ;

    int j;
	double mu, gej, vej;
	double Qeh = 0.;
	for (j = 0; j < (N_SPC - 1); j++)
	{
		mu = mass[N_SPC - 1] * mass[j] / (mass[N_SPC - 1] + mass[j]);
		gej = sqrt(8 * kB * T /* * theta */ / (mu * PI));
		vej = n[N_SPC - 1] * n[j] * Qmpil( gasmix, 1, 1, (N_SPC - 1), j ) * gej;
		Qeh += (3. / 2.) * kB * (mu / (mass[N_SPC - 1] + mass[j])) * vej;
	}

	return Qeh * 1.0e3 ;
}


double DevotoTP::Dij ( GasMixture* gasmix , int order , int ii , int jj ) {

    int N = gasmix->getN() ; 

    std::vector<double> n = gasmix->getCompositionObj()->compositions(1.e-18) ;
    std::vector<double> mass = gasmix->masses(1.e+3) ;

    std::vector<std::vector<double>> QQ ( (N*order)+1 , 
        std::vector<double>((N*order)+1) ) ;
    
    std::vector<std::vector<double>> Q ( (N*order) , 
        std::vector<double>((N*order)) ) ;
         
    double Ti = gasmix->getTemperature() ; 
    double Tj = gasmix->getTemperature() ; 

    if ( ii == N-1 )
        Ti *= gasmix->theta->get() ; 
    if ( jj == N-1 )
        Tj *= gasmix->theta->get() ; 
    
    for( int m = 0 ; m < order ; m++) {
        for( int p = 0 ; p < order ; p++) {
            for( int i = 0 ; i < N ; i++) {
                for( int j = 0 ; j < N ; j++) {

                    QQ[m*N+i][p*N+j] = qmpij(gasmix,m,p,i,j) ;
                    Q [m*N+i][p*N+j] = QQ[m*N+i][p*N+j] ;
                }
            }
        }
    }

    int m = 0 ; 
    for (int i = 0; i < N; i++)
        QQ[m*N+i][N*order] = ((delta(i,jj)-delta(i,ii))) ;
    
    int p = 0 ;
    for (int j = 0; j < N; j++)
        QQ[N*order][p*N+j] = ((delta(j,ii)));
    
    double ntot = 0 ; double rho = 0 ; 
    for (int i = 0; i < N; i++) {
        ntot += n[i] ; 
        rho += n[i] * mass[i] ; 
    }
    
    double detQQ = DetLU(QQ,QQ.size()) ;
    double detQ = DetLU(Q,Q.size()) ;

    double D ;
    D = ( ( 3. * rho * n[ii] ) / ( 2. * ntot * mass[jj] ) ) * 
        sqrt ( ( 2. * PI * kB * Ti ) / ( mass[ii] ) ) * 
            ( detQQ / detQ )  ;
    
    return D * 1.e-12 ; 

}


double DevotoTP::DiT(GasMixture* gasmix, int order, int ii) {
    
    const std::vector<double> n = gasmix->getCompositionObj()->compositions(1.e-18);
    const std::vector<double> mass = gasmix->masses(1.e+3);  // in kg
    const int N = gasmix->getN();
    const double theta = gasmix->theta->get();

    double T = gasmix->getTemperature();
    if (ii == N - 1)
        T *= theta;

    if (ii == N - 1) {

        // Caso elettrone – determinante ridotto
        std::vector<std::vector<double>> Q(3, std::vector<double>(3));
        std::vector<std::vector<double>> q(4, std::vector<double>(4));

        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                q[i][j] = qsimpmpij(gasmix, i, j);

        Q[0][0] = qsimpmpij(gasmix, 0, 1);
        Q[0][1] = qsimpmpij(gasmix, 0, 2);
        Q[0][2] = qsimpmpij(gasmix, 0, 3);
        Q[1][0] = qsimpmpij(gasmix, 2, 1);
        Q[1][1] = qsimpmpij(gasmix, 2, 2);
        Q[1][2] = qsimpmpij(gasmix, 2, 3);
        Q[2][0] = qsimpmpij(gasmix, 3, 1);
        Q[2][1] = qsimpmpij(gasmix, 3, 2);
        Q[2][2] = qsimpmpij(gasmix, 3, 3);

        double detQ = DetLU(Q, 3);
        double detq = DetLU(q, 4);

        return (15. / 4.) * n[ii] * n[ii] * std::sqrt(2. * PI * mass[ii] * kB * T) * (detQ / detq) * 1.e+3;
        
    } else {
        // Caso specie pesante – matrice estesa
        const int Ntot = order * N;

        std::vector<std::vector<double>> Q(Ntot + 1, std::vector<double>(Ntot + 1));
        std::vector<std::vector<double>> q(Ntot, std::vector<double>(Ntot));

        for (int m = 0; m < order; ++m) {
            for (int p = 0; p < order; ++p) {
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        double val = qmpij(gasmix, m, p, i, j);
                        Q[m * N + i][p * N + j] = val;
                        q[m * N + i][p * N + j] = val;
                    }
                }
            }
        }

        // Vincolo: conservazione della massa (riga extra)
        int m = 1; // tipicamente m=1 (come nel codice originale)
        for (int i = 0; i < N; ++i)
            Q[m * N + i][Ntot] = n[i];

        // Vincolo: selezione della specie ii (colonna extra)
        int p = 0; // tipicamente p=0
        for (int j = 0; j < N; ++j)
            Q[Ntot][p * N + j] = (j == ii ? 1.0 : 0.0);

        double detQ = DetLU(Q, Ntot + 1);
        double detq = DetLU(q, Ntot);

        return (15. / 4.) * n[ii] * std::sqrt(2. * PI * mass[ii] * kB * T) * (detQ / detq) * 1.e+3;
    }
}

double DevotoLteLambdaR::Dbinij(int i, int j , GasMixture* gasmix){

    std::vector<double> masses = gasmix->masses(1.e+3) ; // g
    double mi = masses[i] ; 
    double mj = masses[j] ; 
    double mu = mi*mj/(mi+mj) ; 

    double T = gasmix->getTemperature() ;
    double P = gasmix->getPressure()*1.e-3 ; // g / (micron s^2)

    double ntot = gasmix->getCompositionObj()->ntot()*1.e-18 ; 

    // formula binary diffusion coefficients from Transport Properties of Ionized Monatomic Gases - R.S.Devoto 1965
    return (3./(16.* ntot)) * sqrt(2.*std::numbers::pi*kB*T/mu) * (1./Qmpil(gasmix,1,1,i,j)) ; // micron^2 / s
}

void DevotoLteLambdaR::computeTransport(GasMixture* gasmix){

    std::vector<double> nbefore = gasmix->getCompositionObj()->compositions(1.e-18) ; // #/micron^3

    double Tdbg = gasmix->getTemperature() ; 

    /* thermodynamics solve with the definition given to the composition
    object, so, if I give Rat, the composition will have just the electron 
    contribution in the debye length when computing thermodynamic properties and
    the Ghoroui definition composition HAS to be called if the composition used for 
    transport properties has to be with ghoroui. */
    // gasmix->getCompositionObj()->setDebyeModel("Rat2002Th") ;    
    /* composition computed (and left) with rat inside thermodynamics */
    thermo->computeThermodynamics(*gasmix) ; 

    // debug composition
    std::vector<double> nafter = gasmix->getCompositionObj()->compositions(1.e-18) ; // #/micron^3
    
    /* composition reset to ghouroui formulation */
    // gasmix->getCompositionObj()->setDebyeModel("Ghourui") ; 
    // gasmix->getCompositionObj()->CompositionSolve(gasmix,gasmix) ;
    
    DevotoTP::computeTransport(gasmix) ; 
    
    // Debug to see if reactive contribution has been computed with ideal composition
    // gasmix->setCompositionSolver(new GodinTrepSahaSolver(gasmix,gasmix,gasmix->getCompositionObj()->getPfBox())) ;
    // gasmix->getCompositionObj()->CompositionSolve(gasmix,gasmix) ; 

    nafter = gasmix->getCompositionObj()->compositions(1.e-18) ; // #/micron^3

    int N = gasmix->getN() ; 

    int M = gasmix->getM() ; 

    int L = N - M ;

    std::vector<double> mass = gasmix->masses(1.e+3) ; // g 

    std::vector<double> n = gasmix->getCompositionObj()->compositions(1.e-18) ; // #/micron^3

    double ntot = gasmix->getCompositionObj()->ntot() * 1.e-18 ; // #/micron^3

    double T = gasmix->getTemperature() ; // K

    double P = gasmix->getPressure()*1.e-3 ; // g / (micron s^2)

    std::vector<std::vector<double>> v = gasmix->getCompositionObj()->getReactionMatrixV() ; 
    std::vector<int> b = gasmix->getCompositionObj()->getBaseIndicesB();
    std::vector<int> bs = gasmix->getCompositionObj()->getNonBaseIndicesBs();

    std::vector<std::vector<double>> vv ( L, std::vector<double>(N,0.) );

    for (size_t i = 0; i < L; i++) {

        vv[i][bs[i]] = -1. ; 
        for (size_t j = 0; j < M; j++)
            vv[i][b[j]] = v[i][j] ;         
    
    }

    std::vector<double> hi (N,0.) ; 
    for (size_t i = 0; i < N; i++) 
        hi[i] = thermo->h(i) ; // J 
    
    for (size_t i = 0; i < N; i++) 
        hi[i] *= 1.e+15 ; // g micron^2/s^2
        
    std::vector<double> deltaHi (L,0.) ; 
    for (size_t i = 0; i < L; i++)
        for (size_t j = 0; j < N; j++)
            deltaHi[i] += vv[i][j]*hi[j] ; 
       
    std::vector<double> x(N,0.) ; 
    for (size_t i = 0; i < N; i++)
        x[i] = n[i]/ntot ; 
    
    std::vector<std::vector<double>> A (L,std::vector<double>(L,0.)) ; 
    for (size_t i = 0; i < L; i++)
        for (size_t j = 0; j < L; j++)
            for (size_t k = 0; k < N-1; k++)
                for (size_t l = k+1; l < N; l++)
                    A[i][j] += ((1.)/( ntot*Dbinij(k,l,gasmix) ) ) * x[k] * x[l] * ((vv[i][k]/x[k])-(vv[i][l]/x[l])) * ((vv[j][k]/x[k])-(vv[j][l]/x[l]));

    std::vector<std::vector<double>> AH(L+1,std::vector<double>(L+1,0.)) ; 
    for (size_t i = 0; i < L; i++) {

        AH[i][L] = deltaHi[i] ; 
        AH[L][i] = deltaHi[i] ; 

        for (size_t j = 0; j < L; j++)
            AH[i][j] = A[i][j] ; 
            
    }
    
    double num = DetLU(AH,L+1) ; 
    double den = DetLU(A,L) ; 

    double lambdaR = (-1./(kB*T*T))*(num/den) ;

    lambdaR *= 1e-13 * 1e4 ; // back to W/mK

    // std::cout<< T << "  " << lambdaR << std::endl ; 

    double dT = 20.;
    
    std::vector<double> n0 = n;
    for (size_t i = 0; i < N; i++)
        n0[i] *= 1.e+18 ; 
    
    gasmix->setT(T+dT) ; 
    thermo->computeThermodynamics(*gasmix) ; 
    std::vector<double> hf (N,0.) ; 
    for (size_t i = 0; i < N; i++)
        hf[i] = thermo->h(i)*1.e+15 ; 
    
    gasmix->getCompositionObj()->setn0(n0) ; 

    gasmix->setT(T-dT) ; 
    thermo->computeThermodynamics(*gasmix) ; 
    std::vector<double> hb (N,0.) ; 
    for (size_t i = 0; i < N; i++)
        hb[i] = thermo->h(i)*1.e+15 ; 
    
    gasmix->setT(T) ; 
    gasmix->getCompositionObj()->setn0(n0) ; 

    std::vector<double> cpintj (N,0.) ; 
    for (size_t i = 0; i < N; i++)
        cpintj[i] = ((hf[i]-hb[i])/(2.*dT))-(5./2.)*kB ; 
        
    double sum1 = 0.;
    double sum2 = 0.; 
    for (size_t j = 0; j < N; j++){

        sum2 = 0.;
        for (size_t i = 0; i < N; i++)
            sum2 += x[i]/Dbinij(i,j,gasmix);
        
        sum1 += x[j]*cpintj[j]/sum2 ; 
    }
    double lambdaInt = sum1*ntot*1.e-9;

    double mu    = Tp[2];
    double sigma = Tp[3];
    double Qeh   = Tp[4];

    Tp.resize(7) ;
    Tp[2] = lambdaR ;
    Tp[3] = lambdaInt ;
    Tp[4] = mu ;
    Tp[5] = sigma ;
    Tp[6] = Qeh ;
    
}
