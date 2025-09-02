 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include"ZMCoefficients.h"
#include"GasMixture.h"

double ZMCoefficients::b10 ( GasMixture* gasmix ) {

    double Q11, Q00, Q01, Q10 ; 
    
    Q11 = Q1[1][1] ;
    Q00 = Q1[0][0] ;
    Q01 = Q1[0][1] ;
    Q10 = Q1[1][0] ;

    return ( 5.*ne*Q11 ) / ( (Q00*Q11) - (Q10*Q01) ) ; 

}

double ZMCoefficients::b11 ( GasMixture* gasmix ) {

    double Q11, Q00, Q01, Q10 ; 
    
    Q11 = Q1[1][1] ;
    Q00 = Q1[0][0] ;
    Q01 = Q1[0][1] ;
    Q10 = Q1[1][0] ;

    return ( 5.*ne*Q10 ) / ( (Q10*Q01) - (Q00*Q11) ) ; 

}
 
double ZMCoefficients::c1psk ( GasMixture* gasmix , int epsilon , int s , int k , int pp ) {
    
    std::vector<std::vector<double>> QQ ( epsilon + 1 ,
        std::vector<double>(epsilon + 1) ) ; 
    
    std::vector<std::vector<double>> Q  ( epsilon , 
        std::vector<double>(epsilon) ) ; 

    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++){
        for (int p = 0; p < epsilon; p++){
            Q[m][p] = Q1[m][p] ; 
            QQ[m][p] = Q[m][p] ;
            QQ[epsilon][p] = delta(pp,p);
            if ( m == 0 )
                QQ[m][epsilon] = -3.*sqrt((kB*Te)/(2.*me))*((delta(N-1,s)-delta(N-1,k))/(kB*Te)) ;
            else 
                QQ[m][epsilon] = 0. ; 
        }
    }

    QQ[epsilon][epsilon] = 0. ; 

    double DetQQ = DetLU( QQ , QQ.size() ) ;
    double DetQ  = DetLU( Q  ,  Q.size() ) ; 

    return -(DetQQ/DetQ) ; 

}
double ZMCoefficients::ci0 ( GasMixture* gasmix , int epsilon , int s , int k , int i ) {
    
    if( i == N-1 )
        return c1psk ( gasmix , epsilon , s , k , 0 ) ;

    std::vector<std::vector<double>> QQ ( 1+(epsilon*(N-1)) , 
        std::vector<double>( 1+(epsilon*(N-1)) ) ) ;
    
    std::vector<std::vector<double>> Q ( epsilon*(N-1) , 
        std::vector<double>( epsilon*(N-1) ) )         ;

    #pragma omp parallel for collapse(4)
    for ( int m = 0 ; m < epsilon ;   m++ ) {
        for ( int p = 0 ; p < epsilon ; p++ ) {
            for ( int l = 0 ; l < N-1 ;   l++ ) {       
                for ( int j = 0 ; j < N-1 ; j++ ) {
                    Q[m*(N-1)+l][p*(N-1)+j] = Qtilde[m*(N-1)+l][p*(N-1)+j] ;   
                    QQ[m*(N-1)+l][p*(N-1)+j] = Q[m*(N-1)+l][p*(N-1)+j] ;
                    if ( p == 0 )
                        QQ[epsilon*(N-1)][p*(N-1)+j] = delta(i,j) ;
                    else 
                        QQ[epsilon*(N-1)][p*(N-1)+j] = 0. ;
                }
            }
        }
    }
    
    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++) {
        for (int l = 0; l < N-1; l++) {
            double sum = 0. ; 
            sum = 0. ;
            for (int pp = 0; pp < epsilon; pp++) {
                sum += qmpi1Bar(gasmix,m,pp,l)*c1psk(gasmix,epsilon,s,k,pp) ;
            }
            if ( m == 0 ) {
                QQ[m*(N-1)+l][epsilon*(N-1)] = -3.*(sqrt((kB*T)/(2.*mass[l]))*
                    ((delta(l,s)-delta(l,k))/(kB*T))) - sum ;
            }else {
                QQ[m*(N-1)+l][epsilon*(N-1)] = -sum ; 
            }
        }
    }
    
    QQ[epsilon*(N-1)][epsilon*(N-1)] = 0. ;

    double detQQ = DetLU ( QQ , QQ.size() ) ;    
    double detQ  = DetLU (  Q ,  Q.size() ) ;

    return -(detQQ/detQ) ; 

}

double ZMCoefficients::c11   ( GasMixture* gasmix , int epsilon , int s , int k ) {

    std::vector<std::vector<double>> QQ ( epsilon ,
        std::vector<double>(epsilon) ) ; 
    
    std::vector<std::vector<double>> Q  ( epsilon , 
        std::vector<double>(epsilon) ) ; 

    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++) {
        for (int p = 0; p < epsilon; p++) {
            Q[m][p] = Q1[m][p] ;
            QQ[m][p] = Q[m][p] ;
            if ( m == 0 )
                QQ[m][1] = -3.*sqrt((kB*Te)/(2.*me))*((delta(N-1,s)-delta(N-1,k))/(kB*Te)) ;
            else 
                QQ[m][1] = 0. ; 
        }
    }
    
    double detQQ = DetLU(QQ , QQ.size()) ;
    double detQ  = DetLU(Q  , Q.size() ) ;    

    return (detQQ)/(detQ) ;

}

double ZMCoefficients::ci1 ( GasMixture* gasmix , int epsilon , int s , int k , int i ) {
    
    if( i == N-1 )
        return c11 ( gasmix , epsilon , s , k ) ;

    std::vector<std::vector<double>> QQ ( 1+(epsilon*(N-1)) , 
        std::vector<double>( 1+(epsilon*(N-1)) ) ) ;
    
    std::vector<std::vector<double>> Q ( epsilon*(N-1) , 
        std::vector<double>( epsilon*(N-1) ) )         ;

    #pragma omp parallel for collapse(4)
    for ( int m = 0 ; m < epsilon ;   m++ ) {    
        for ( int p = 0 ; p < epsilon ; p++ ) {
            for ( int l = 0 ; l < N-1 ;   l++ ) {
                for ( int j = 0 ; j < N-1 ; j++ ) {                
                    Q[m*(N-1)+l][p*(N-1)+j] = Qtilde[m*(N-1)+l][p*(N-1)+j] ;
                    QQ[m*(N-1)+l][p*(N-1)+j] = Q[m*(N-1)+l][p*(N-1)+j] ;
                    if (p==1)
                        QQ[epsilon*(N-1)][p*(N-1)+j] = delta(i,j) ;
                    else 
                        QQ[epsilon*(N-1)][p*(N-1)+j] = 0. ;
                }
            }
        }   
    }

    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++) {
        for (int l = 0; l < N-1; l++) {
            double sum = 0. ; 
            sum = 0. ;
            for (int pp = 0; pp < epsilon; pp++) {
                sum += qmpi1Bar(gasmix,m,pp,l)*c1psk(gasmix,epsilon,s,k,pp) ;
            }
            if (m == 0) {
                QQ[m*(N-1)+l][epsilon*(N-1)] = -3.*(sqrt((kB*T)/(2.*mass[l]))*
                    ((delta(l,s)-delta(l,k))/(kB*T))) - sum ;
            }else {
                QQ[m*(N-1)+l][epsilon*(N-1)] = -sum ; 
            }
        }
    }
    
    QQ[epsilon*(N-1)][epsilon*(N-1)] = 0. ;

    double detQQ = DetLU ( QQ , QQ.size() ) ;    
    double detQ  = DetLU (  Q ,  Q.size() ) ;

    return -(detQQ/detQ) ; 

}

double ZMCoefficients::aip   ( GasMixture* gasmix , int epsilon , int i , int pp ) {

    std::vector<std::vector<double>> QQ ( epsilon + 1 ,
        std::vector<double>(epsilon + 1) ) ; 
    
    std::vector<std::vector<double>> Q  ( epsilon , 
        std::vector<double>(epsilon) ) ; 

    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++){
        for (int p = 0; p < epsilon; p++){
            Q[m][p] = Q1[m][p] ;
            QQ[m][p] = Q[m][p] ;
            QQ[epsilon][p] = delta(pp,p);
            if ( m == 1 )
                QQ[m][epsilon] = -((15./2.)*ne*sqrt((kB*Te)/(2.*me))) ;
            else 
                QQ[m][epsilon] = 0. ; 
        }   
    }
    
    QQ[epsilon][epsilon] = 0. ; 

    double DetQQ = DetLU( QQ , QQ.size() ) ;
    double DetQ  = DetLU( Q  ,  Q.size() ) ; 

    return -(DetQQ/DetQ) ; 

}

double ZMCoefficients::ai0   ( GasMixture* gasmix , int epsilon , int i ) {
    
    if( i == N-1 )
        return aip ( gasmix , epsilon , N-1 , 0 ) ;

    std::vector<std::vector<double>> QQ ( 1+(epsilon*(N-1)) , 
        std::vector<double>( 1+(epsilon*(N-1)) ) ) ;
    
    std::vector<std::vector<double>> Q ( epsilon*(N-1) , 
        std::vector<double>( epsilon*(N-1) ) )         ;

    #pragma omp parallel for collapse(4)
    for ( int m = 0 ; m < epsilon ;   m++ ) {
        for ( int p = 0 ; p < epsilon ; p++ ) {
            for ( int l = 0 ; l < N-1 ;   l++ ) {
                for ( int j = 0 ; j < N-1 ; j++ ) {
                    Q[m*(N-1)+l][p*(N-1)+j] = Qtilde[m*(N-1)+l][p*(N-1)+j] ;
                    QQ[m*(N-1)+l][p*(N-1)+j] = Q[m*(N-1)+l][p*(N-1)+j] ;
                if ( p==0 ) 
                    QQ[epsilon*(N-1)][p*(N-1)+j] = delta(i,j) ;
                else 
                    QQ[epsilon*(N-1)][p*(N-1)+j] = 0. ;
                }
            }
        }
    }
    
    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++) {
        for (int l = 0; l < N-1; l++) {
            double sum = 0. ; 
            sum = 0. ;
            double nl = n[l] ; 
            double ml = mass[l] ;
            for (int pp = 0; pp < epsilon; pp++) {
                sum += qmpi1Bar(gasmix,m,pp,l)*aip(gasmix,epsilon,l,pp) ;
            }
            if (m==1) {
                QQ[m*(N-1)+l][epsilon*(N-1)] = - ((15./2.) * nl * ( sqrt ((kB*T) / (2.*ml)) )) - sum;
            }else {
                QQ[m*(N-1)+l][epsilon*(N-1)] = - sum ; 
            }
        }
    }
    
    QQ[epsilon*(N-1)][epsilon*(N-1)] = 0. ;

    double detQQ = DetLU ( QQ , QQ.size() ) ;    
    double detQ  = DetLU (  Q ,  Q.size() ) ;

    return -(detQQ/detQ) ; 

}

double ZMCoefficients::a11 ( GasMixture* gasmix , int epsilon ) {


    std::vector<std::vector<double>> QQ ( epsilon ,
        std::vector<double>(epsilon) ) ; 
    
    std::vector<std::vector<double>> Q  ( epsilon , 
        std::vector<double>(epsilon) ) ; 
    
    #pragma omp parallel for collapse(2) 
    for (int m = 0; m < epsilon; m++){
        for (int p = 0; p < epsilon; p++){
            Q[m][p] = Q1[m][p] ;
            QQ[m][p] = Q[m][p] ;
            if ( m == 1 )
                QQ[m][1] = -((15./2.)*ne*sqrt((kB*Te)/(2.*me))) ;
            else 
                QQ[m][1] = 0. ; 
        }
    }
    
    double detQQ = DetLU(QQ , QQ.size()) ;
    double detQ  = DetLU(Q  , Q.size() ) ;    

    return (detQQ)/(detQ) ;

}

double ZMCoefficients::ai1 ( GasMixture* gasmix , int epsilon , int i ) {

    if( i == N-1 )
        return a11 ( gasmix , epsilon ) ;

    std::vector<std::vector<double>> QQ ( 1+(epsilon*(N-1)) , 
        std::vector<double>( 1+(epsilon*(N-1)) ) ) ;
    
    std::vector<std::vector<double>> Q ( epsilon*(N-1) , 
        std::vector<double>( epsilon*(N-1) ) )         ;

    #pragma omp parallel for collapse(4)
    for ( int m = 0 ; m < epsilon ;   m++ ) {
        for ( int p = 0 ; p < epsilon ; p++ ) { 
            for ( int l = 0 ; l < N-1 ;   l++ ) {
                for ( int j = 0 ; j < N-1 ; j++ ) {                         
                    Q[m*(N-1)+l][p*(N-1)+j] = Qtilde[m*(N-1)+l][p*(N-1)+j] ;
                    QQ[m*(N-1)+l][p*(N-1)+j] = Q[m*(N-1)+l][p*(N-1)+j] ;
                    if (p == 1)
                        QQ[epsilon*(N-1)][p*(N-1)+j] = delta(i,j) ;
                    else 
                        QQ[epsilon*(N-1)][p*(N-1)+j] = 0. ;
                }
            }
        }
    }
    
    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++) {
        for (int l = 0; l < N-1; l++) {
            double sum = 0. ; 
            sum = 0. ;
            for (int pp = 0; pp < epsilon; pp++) {
                sum += qmpi1Bar(gasmix,m,pp,l)*aip(gasmix,epsilon,l,pp) ;
            }
            double nl = n[l] ; 
            double ml = mass[l] ;
            if (m == 1) {
                QQ[m*(N-1)+l][epsilon*(N-1)] = -((15./2.)*nl*(sqrt((kB*T)/(2.*ml)))) - sum;
            }else {
                QQ[m*(N-1)+l][epsilon*(N-1)] = -sum ; 
            }
        }
    }
    
    QQ[epsilon*(N-1)][epsilon*(N-1)] = 0. ;

    double detQQ = DetLU ( QQ , QQ.size() ) ;    
    double detQ  = DetLU (  Q ,  Q.size() ) ;

    return -(detQQ/detQ) ; 

}

double ZMCoefficients::e1psk   ( GasMixture* gasmix , int epsilon , int s , int k , int pp ) {

    std::vector<std::vector<double>> QQ ( epsilon + 1 ,
        std::vector<double>(epsilon + 1) ) ; 
    
    std::vector<std::vector<double>> Q  ( epsilon , 
        std::vector<double>(epsilon) ) ; 

    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++){
        for (int p = 0; p < epsilon; p++){
            Q[m][p] = Q1[m][p] ;
            QQ[m][p] = Q[m][p] ;
            if ( m == 0 )
                QQ[m][epsilon] = -3.*sqrt((kB*Te)/(2.*me))*((delta(N-1,s)-delta(N-1,k))/(kB*T)) ;
            else 
                QQ[m][epsilon] = 0. ; 
            QQ[epsilon][p] = delta(pp,p);
        }
    }
    
    QQ[epsilon][epsilon] = 0. ; 

    double DetQQ = DetLU( QQ , QQ.size() ) ;
    double DetQ  = DetLU( Q  ,  Q.size() ) ; 

    return -(DetQQ/DetQ) ; 

}

double ZMCoefficients::ei0 ( GasMixture* gasmix , int epsilon , int s , int k , int ii ) {

    if( ii == N-1 )
        return e1psk ( gasmix , epsilon , s , k , 0 ) ;

    std::vector<std::vector<double>> QQ ( 1+(epsilon*(N-1)) , 
        std::vector<double>( 1+(epsilon*(N-1)) ) ) ;
    
    std::vector<std::vector<double>> Q ( epsilon*(N-1) , 
        std::vector<double>( epsilon*(N-1) ) )         ;

    #pragma omp parallel for collapse(4)
    for ( int m = 0 ; m < epsilon ;   m++ ) {    
        for ( int p = 0 ; p < epsilon ; p++ ) {
            for ( int l = 0 ; l < N-1 ;   l++ ) {   
                for ( int j = 0 ; j < N-1 ; j++ ) {                   
                    Q[m*(N-1)+l][p*(N-1)+j] = Qtilde[m*(N-1)+l][p*(N-1)+j] ;
                    QQ[m*(N-1)+l][p*(N-1)+j] = Q[m*(N-1)+l][p*(N-1)+j] ;
                    if ( p == 0 )
                        QQ[epsilon*(N-1)][p*(N-1)+j] = delta(ii,j) ;
                    else 
                        QQ[epsilon*(N-1)][p*(N-1)+j] = 0. ;
                }
            }
        }
    }

    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++) {
        for (int l = 0; l < N-1; l++) {
            double sum = 0. ; 
            sum = 0. ;
            for (int pp = 0; pp < epsilon; pp++) {
                sum += qmpi1Bar(gasmix,m,pp,l)*e1psk(gasmix,epsilon,s,k,pp) ; 
            }
            if ( m == 0 ) {
                QQ[m*(N-1)+l][epsilon*(N-1)] = -3.*(sqrt((kB*T)/(2.*mass[l])) *
                    (theta*(delta(l,s)-delta(l,k))/(kB*T))) - sum ;
            } else {
                QQ[m*(N-1)+l][epsilon*(N-1)] = -sum ; 
            }
        }
    }
    
    QQ[epsilon*(N-1)][epsilon*(N-1)] = 0. ;

    double detQQ = DetLU ( QQ , QQ.size() ) ;    
    double detQ  = DetLU (  Q ,  Q.size() ) ;

    return -(detQQ/detQ) ; 

}

double ZMCoefficients::e11sk ( GasMixture* gasmix , int epsilon , int s , int k ) {
    
    std::vector<std::vector<double>> QQ ( epsilon ,
        std::vector<double>(epsilon) ) ; 
    
    std::vector<std::vector<double>> Q  ( epsilon , 
        std::vector<double>(epsilon) ) ; 

    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++) {
        for (int p = 0; p < epsilon; p++) {
            Q[m][p] = Q1[m][p] ;
            QQ[m][p] = Q[m][p] ;
            if ( m == 0 )
                QQ[m][1] = -3.*sqrt((kB*Te)/(2.*me))*((delta(N-1,s)-delta(N-1,k))/(kB*T)) ;
            else 
                QQ[m][1] = 0. ; 
        }
    }

    double detQQ = DetLU(QQ , QQ.size()) ;    
    double detQ  = DetLU(Q  , Q.size() ) ;    

    return (detQQ)/(detQ) ;

}

double ZMCoefficients::ei1   ( GasMixture* gasmix , int epsilon , int s , int k , int ii ) {
    
    if( ii == N-1 )
        return e11sk ( gasmix , epsilon , s , k ) ;

    std::vector<std::vector<double>> QQ ( 1+(epsilon*(N-1)) , 
        std::vector<double>( 1+(epsilon*(N-1)) ) ) ;
    
    std::vector<std::vector<double>> Q ( epsilon*(N-1) , 
        std::vector<double>( epsilon*(N-1) ) )         ;

    #pragma omp parallel for collapse(4)
    for ( int m = 0 ; m < epsilon ;   m++ ) {  
        for ( int p = 0 ; p < epsilon ; p++ ) {   
            for ( int l = 0 ; l < N-1 ;   l++ ) {          
                for ( int j = 0 ; j < N-1 ; j++ ) {                         
                    Q[m*(N-1)+l][p*(N-1)+j] = Qtilde[m*(N-1)+l][p*(N-1)+j] ;
                    QQ[m*(N-1)+l][p*(N-1)+j] = Q[m*(N-1)+l][p*(N-1)+j] ;
                    if ( p == 1 )
                        QQ[epsilon*(N-1)][p*(N-1)+j] = delta(ii,j) ;
                    else 
                        QQ[epsilon*(N-1)][p*(N-1)+j] = 0. ;
                }
            }
        }
    }
            
    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++) {
        for (int l = 0; l < N-1; l++) {
            double sum = 0. ; 
            sum = 0. ;
            for (int pp = 0; pp < epsilon; pp++) {
                sum += qmpi1Bar(gasmix,m,pp,l)*e1psk(gasmix,epsilon,s,k,pp) ;
            }
            if (m==0) {
                QQ[m*(N-1)+l][epsilon*(N-1)] = -3.*(sqrt((kB*T)/(2.*mass[l]))*
                    (theta*(delta(l,s)-delta(l,k))/(kB*T))) - sum ;
            }else {
                QQ[m*(N-1)+l][epsilon*(N-1)] = -sum ; 
            }
        }
    }
    
    QQ[epsilon*(N-1)][epsilon*(N-1)] = 0. ;

    double detQQ = DetLU ( QQ , QQ.size() ) ;    
    double detQ  = DetLU (  Q ,  Q.size() ) ;

    return -(detQQ/detQ) ; 

}

double ZMCoefficients::f10 ( GasMixture* gasmix , int epsilon ) {

    std::vector<std::vector<double>> QQ( epsilon , std::vector<double>(epsilon) ) ;
    std::vector<std::vector<double>> Q ( epsilon , std::vector<double>(epsilon) ) ;

    double C = -(15./2.)*ne*sqrt((kB*Te)/(2.*me)) ;
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < epsilon; i++) {
        for (int j = 0; j < epsilon; j++) {
            
            Q[i][j] = Q1[i][j] ; 

            if ( j == 0 ) {
                if ( i == 1 )
                    QQ[i][j] = C ;          
                else 
                    QQ[i][j] = 0.;   
            } else {
                QQ[i][j] = Q[i][j] ; 
            } 
        }
    }

    double detQQ = DetLU ( QQ , QQ.size() ) ;
    double detQ  = DetLU ( Q  ,  Q.size() ) ;    

    return detQQ/detQ ;   

}

double ZMCoefficients::f11 ( GasMixture* gasmix , int epsilon ) {

    std::vector<std::vector<double>> QQ( epsilon , std::vector<double>(epsilon) ) ;
    std::vector<std::vector<double>> Q ( epsilon , std::vector<double>(epsilon) ) ;

    double C = -(15./2.)*ne*sqrt((kB*Te)/(2.*me)) ;

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < epsilon; i++) {
        for (int j = 0; j < epsilon; j++) {
            
            Q[i][j] = Q1[i][j] ;

            if (j==1) {
                if (i==1)
                    QQ[i][j] = C ;          
                else 
                    QQ[i][j] = 0.;   
            } else {
                QQ[i][j] = Q[i][j] ; 
            } 
        }
    }

    double detQQ = DetLU ( QQ , QQ.size() ) ;
    double detQ  = DetLU ( Q  ,  Q.size() ) ;    

    return detQQ/detQ ;   
}

double ZMCoefficients::fi0 ( GasMixture* gasmix , int epsilon , int ii ) {

    if( ii == N-1 )
        return f10 ( gasmix , epsilon ) ;

    std::vector<std::vector<double>> QQ ( 1+(epsilon*(N-1)) , 
        std::vector<double>( 1+(epsilon*(N-1)) ) ) ;
    
    std::vector<std::vector<double>> Q ( epsilon*(N-1) , 
        std::vector<double>( epsilon*(N-1) ) )         ;

    #pragma omp parallel for collapse(4)
    for ( int m = 0 ; m < epsilon ;   m++ )    {
        for ( int p = 0 ; p < epsilon ; p++ )    {
            for ( int l = 0 ; l < N-1 ;   l++ )     {     
                for ( int j = 0 ; j < N-1 ; j++ )    {                     
                    Q[m*(N-1)+l][p*(N-1)+j] = Qtilde[m*(N-1)+l][p*(N-1)+j] ;
                    QQ[m*(N-1)+l][p*(N-1)+j] = Q[m*(N-1)+l][p*(N-1)+j] ;
                    if ( p==0 ) 
                        QQ[epsilon*(N-1)][p*(N-1)+j] = delta(ii,j) ;
                    else 
                        QQ[epsilon*(N-1)][p*(N-1)+j] = 0. ;
                }
            }
        }
    }

    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++) {
        for (int l = 0; l < N-1; l++) {
            
            double sum = 0. ; 
            sum = 0. ;
            for (int pp = 0; pp < epsilon; pp++) {
                sum += qmpi1Bar(gasmix,m,pp,l)*aip(gasmix,epsilon,l,pp) ;
            }
    
            QQ[m*(N-1)+l][epsilon*(N-1)] = - sum;
            
        }
    }
    
    QQ[epsilon*(N-1)][epsilon*(N-1)] = 0. ;

    double detQQ = DetLU ( QQ , QQ.size() ) ;    
    double detQ  = DetLU (  Q ,  Q.size() ) ;

    return -(detQQ/detQ) ;

}

double ZMCoefficients::fi1 ( GasMixture* gasmix , int epsilon , int ii ) {
    
    if( ii == N-1 )
        return f11 ( gasmix , epsilon ) ;

    std::vector<std::vector<double>> QQ ( 1+(epsilon*(N-1)) , 
        std::vector<double>( 1+(epsilon*(N-1)) ) ) ;
    
    std::vector<std::vector<double>> Q ( epsilon*(N-1) , 
        std::vector<double>( epsilon*(N-1) ) )         ;

    #pragma omp parallel for collapse(4)
    for ( int m = 0 ; m < epsilon ;   m++ )    {
        for ( int p = 0 ; p < epsilon ; p++ )    {
            for ( int l = 0 ; l < N-1 ;   l++ )    {      
                for ( int j = 0 ; j < N-1 ; j++ )   {                      
                    Q[m*(N-1)+l][p*(N-1)+j] = Qtilde[m*(N-1)+l][p*(N-1)+j] ;
                    QQ[m*(N-1)+l][p*(N-1)+j] = Q[m*(N-1)+l][p*(N-1)+j] ;
                    if (p == 1)
                        QQ[epsilon*(N-1)][p*(N-1)+j] = delta(ii,j) ;
                    else 
                        QQ[epsilon*(N-1)][p*(N-1)+j] = 0. ;
                }
            }
        }
    }

    #pragma omp parallel for collapse(2)
    for (int m = 0; m < epsilon; m++) {
        for (int l = 0; l < N-1; l++) {
            double sum = 0. ; 
            sum = 0. ;
            for (int pp = 0; pp < epsilon; pp++) {
                sum += qmpi1Bar(gasmix,m,pp,l)*aip(gasmix,epsilon,l,pp) ;
            }
                
            QQ[m*(N-1)+l][epsilon*(N-1)] = -sum ; 
        }
    }
    
    QQ[epsilon*(N-1)][epsilon*(N-1)] = 0. ;

    double detQQ = DetLU ( QQ , QQ.size() ) ;    
    double detQ  = DetLU (  Q ,  Q.size() ) ;

    return -(detQQ/detQ) ; 

}

double ZMCoefficients::wi ( GasMixture* gasmix , int i ) {
    
    //[p] = Pa = kg / m s^2
    double p = gasmix->getPressure() ; 
    //[p] = g / micron s^2
    p *= 1.e+3 * 1.e-6 ; 

    double x1 = n.back() / ntot ;
    double xi = n[i] / ntot ;

    double D = 1 + x1*( theta - 1 ) ; 

    double w1 = x1*p*(1-x1)*pow(D,-2.) ;
    double wi = -xi*x1*p*pow(D,-2.) ;

    if ( i == N-1 )
        return w1 ; 
    else
        return wi ;
    
}

double ZMCoefficients::lambdaiPrime ( GasMixture* gasmix, int epsilon, int i ) {
    
    double Ti = T ;
    double Th = Ti ;
    
    if (i == N-1)
        Ti *= theta ;
    
    double mi = mass[i] ;
    double ni = n[i] ; 

    double aiUno = ai1(gasmix,epsilon,i) ;

    double l = (-5./4.)*kB*(Ti/Th)*ni*sqrt((2.*kB*Ti)/mi) * aiUno ;
    
    return l ;

}

double ZMCoefficients::lambdaijD ( GasMixture* gasmix, int epsilon, int i , int j ) {

    double ni = n[i] ;
    double nj = n[j] ;
    double mi = mass[i] ;
    double mj = mass[j] ;

    double Ti = T ; 
    if ( i == N-1 )
        Ti *= theta ;
    
    double ciUno = ci1(gasmix,epsilon,j,i,i) ;

    double lijD = (5./4.)*ni*nj*mj*kB*Ti*sqrt((2.*kB*Ti)/mi)*ciUno ; 

    return lijD ;  

}


double ZMCoefficients::lambdaiPrimeTheta ( GasMixture* gasmix, int epsilon, int ii ) {
    
    double ni = n[ii] ;
    double mi = mass[ii] ; 

    double Ti = T ;
    if ( ii == N-1 )
        Ti *= theta ; 
    
    double Fi1 = fi1(gasmix,epsilon,ii) ;

    double sum = 0. ;
    for (int j = 0; j < N; j++){
        sum += ei1(gasmix,epsilon,j,ii,ii) * wi(gasmix,j) ;
    }

    return -(5./4.)*kB*ni*sqrt((2.*kB*Ti)/(mi))*(Fi1 - sum) ;    

}

double ZMCoefficients::DiThetaStar ( GasMixture* gasmix, int epsilon, int ii ) {
    
    double ni = n[ii] ;
    double mi = mass[ii] ;
    double Ti = T ;

    if (ii == N-1)
        Ti = Te ; 
    
    double D = ni*mi*sqrt((kB*Ti)/(2.*mi))*fi0(gasmix,epsilon,ii) ;
    
    return D ;

}

void ZMCoefficients::init ( GasMixture* gasmix ) {

    N = gasmix->getN() ;
    n = gasmix->Comp->compositions(1.e-18) ; 
    mass = gasmix->masses(1.e+3) ;
    T = gasmix->getTemperature() ;
    theta = gasmix->theta->get() ;
    Te = T*theta ;
    me = mass.back() ;
    ne = n.back() ;

    for (int i = 0; i < N; i++){
        ntot += n[i];
        rho += n[i] * mass[i] ;
    }
    
    Qtilde.resize ( 4*(N-1), std::vector<double>(4*(N-1)) ) ; 
    Q1.resize (4,std::vector<double>(4)) ;
    
    #pragma omp parallel for collapse(4) 
    for ( int m = 0 ; m < 4 ;   m++ )    
        for ( int p = 0 ; p < 4 ; p++ )    
            for ( int l = 0 ; l < N-1 ;   l++ )          
                for ( int j = 0 ; j < N-1 ; j++ )                         
                    Qtilde[m*(N-1)+l][p*(N-1)+j] = qmpijBar(gasmix,m,p,l,j) ;
    
    #pragma omp parallel for collapse(2)
    for (int m = 0; m < 4; m++)
        for (int p = 0; p < 4; p++)
            Q1[m][p] = qsimpmpij( gasmix, m, p ) ;

};
