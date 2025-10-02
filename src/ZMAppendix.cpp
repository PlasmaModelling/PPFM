 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include"ZMAppendix.h"
#include"GasMixture.h"

double ZMAppendix::Qmpil ( GasMixture* gasmix , int l , int s , int i , int j ) {

    double Ti, Tj, mi, mj ; 
    
    double Th = gasmix->getTemperature() ;

    Ti = gasmix->getTemperature() ;
    if (dynamic_cast<Electron*>((*gasmix)(i)))
        Ti *= gasmix->theta->get() ; 
    
    Tj = gasmix->getTemperature() ;
    if (dynamic_cast<Electron*>((*gasmix)(j)))
        Tj *= gasmix->theta->get() ; 
    
    mi = (*gasmix)(i)->getMass() * 1.e+03 ;
    mj = (*gasmix)(j)->getMass() * 1.e+03 ;

    double muij = ( mi * mj ) / (mi + mj) ; 
    double Tij = pow( ( (1./(mi+mj)) * ((mi/Tj) + (mj/Ti)) ) , (-1.) ) ; 

    double dim = ((factorial(s+1) * (2*l + 1 - pow(-1,l))) / (4.*(l+1))) *
        sqrt( ( kB*Tij ) / ( 2*PI*muij ) ) ; 
    double q = Appendix::Qmpil(gasmix,l,s,i,j) ; 
    
    return q*dim ;

}
double ZMAppendix::qmpi1V ( GasMixture* gasmix , int m , int p , int i ) {
    
    double sum = 0.;
    double M1,Mi ; 

    double n1 = gasmix->getCompositionObj()->compositions(1.e-18).back() ;
    double ni = gasmix->getCompositionObj()->compositions(1.e-18)[i] ;

    double m1 = gasmix->masses(1.e+3).back() ;
    double mi = gasmix->masses(1.e+3)[i] ;

    int N = gasmix->getN() ;
    int e = N-1 ;

    Mi = mi / (mi+m1) ;
    M1 = m1 / (mi+m1) ;

    switch (m) {
    case 0:
        switch (p) {
        case 0:
            sum += (Mi*M1) * ( (-80./3.)*Qmpil(gasmix,1,1,i,e) + 8.*Qmpil(gasmix,2,2,i,e) );
        break;
        case 1:
            sum +=  (pow(Mi,2)*M1) * ( (-280./3.)*Qmpil(gasmix,1,1,i,e) + 
                (112./3.)*Qmpil(gasmix,1,2,i,e) + 28.*Qmpil(gasmix,2,2,i,e) -
                    8.*Qmpil(gasmix,2,3,i,e)  
            );
        break;
    
        
        }
        break;
    case 1:
        switch (p) {
        case 0:
            sum += ( pow((Mi/M1) , -1) * (pow(Mi,2)*M1) ) * ( (-280./3.)*Qmpil(gasmix,1,1,i,e) + 
                (112./3.)*Qmpil(gasmix,1,2,i,e) + 28.*Qmpil(gasmix,2,2,i,e) -
                    8.*Qmpil(gasmix,2,3,i,e)  
            );
        break;
        case 1:
            sum += (pow(Mi,2.)*pow(M1,2.)) * ( (-1540./3.)*Qmpil(gasmix,1,1,i,e) +
                (784./3.)*Qmpil(gasmix,1,2,i,e) - (128./3.)*Qmpil(gasmix,1,3,i,e) +
                    (602./3.)*Qmpil(gasmix,2,2,i,e) - (56.)*Qmpil(gasmix,2,3,i,e) +
                        8.*Qmpil(gasmix,2,4,i,e)-16.*Qmpil(gasmix,3,3,i,e) 
            );
        break;
    
        
        }
        break;
    }
    return sum*ni*n1 ;
}

double ZMAppendix::qcapmpij( GasMixture* gasmix, int m, int p, int i, int j) {

    std::vector<double> n = gasmix->getCompositionObj()->compositions(1.e-18) ; 
    std::vector<double> mass = gasmix->masses(1.e+03) ;
    int N = gasmix->getN() ;

    double ni = n[i] ; 
    double nk ;
    double sum = 0.;

    double mi = mass[i];
    double mk ;
    double Mi ; 
    double Mk ;
    
    switch (m) {

     case 0:
        switch (p) {
         
         case 0:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += nk * ni *( 
                    
                    (delta(i,j) * 
                        ((80./3.)*(Mi*Mk)*Qmpil(gasmix,1,1,i,k) + 
                        8.*pow(Mk,2.)*Qmpil(gasmix,2,2,i,k) ) ) +
                    
                    (delta(j,k) * (Mi*Mk) *
                        (-(80./3.)*Qmpil(gasmix,1,1,i,k) + 8.*Qmpil(gasmix,2,2,i,k) )
                    )
                ) ;  
                 
            }
            
         break;
         
         case 1:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += nk * ni *( 

                    (delta(i,j) * 
                        ((280./3.)*(Mi*pow(Mk,2.))*Qmpil(gasmix,1,1,i,k) - 
                        (112./3.)*(Mi*pow(Mk,2.))*Qmpil(gasmix,1,2,i,k) + 
                        28.*pow(Mk,3.)*Qmpil(gasmix,2,2,i,k) - 
                        8.*pow(Mk,3.)*Qmpil(gasmix,2,3,i,k)) ) +
                    
                    (delta(j,k) * (pow(Mi,2.)*Mk) * 
                        (-(280./3.)*Qmpil(gasmix,1,1,i,k) + (112./3.)*Qmpil(gasmix,1,2,i,k) +
                        28.*Qmpil(gasmix,2,2,i,k) - 8.*Qmpil(gasmix,2,3,i,k))
                    )
                ) ;      
            }
         break;
        }
     break;

     case 1:
        switch (p) {
         
         case 0:
            for (int k = 0; k < N - 1 ; k++) {
              
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += nk * ni *( 

                    (delta(i,j) * 
                        ((280./3.)*(Mi*pow(Mk,2.))*Qmpil(gasmix,1,1,i,k) - 
                        (112./3.)*(Mi*pow(Mk,2.))*Qmpil(gasmix,1,2,i,k) + 
                        28.*pow(Mk,3.)*Qmpil(gasmix,2,2,i,k) - 
                        8.*pow(Mk,3.)*Qmpil(gasmix,2,3,i,k)) ) +

                    (delta(j,k) * pow((Mi/Mk),-1) * (pow(Mi,2.)*Mk) * 
                        (-(280./3.)*Qmpil(gasmix,1,1,i,k) + (112./3.)*Qmpil(gasmix,1,2,i,k) +
                        28.*Qmpil(gasmix,2,2,i,k) - 8.*Qmpil(gasmix,2,3,i,k))
                    )
                


                ) ;      
            }
            
         break;
         
         case 1:
            for (int k = 0; k < N - 1 ; k++) {
              
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += nk * ni *( 
                    
                    (delta(i,j) * 
                        ( ((1540./3.)*( (4./11.)*pow(Mi,3.)*Mk + (7./11.)*pow(Mk,3.)*Mi )*Qmpil(gasmix,1,1,i,k)) - 
                        ((784./3.)*(Mi*pow(Mk,3.))*Qmpil(gasmix,1,2,i,k)) + 
                        ((128./3.)*(Mi*pow(Mk,3.))*Qmpil(gasmix,1,3,i,k) )+  
                        ((602./3.)*((22./43.)*pow(Mi,2)*pow(Mk,2.) + 
                        (21./43.)*pow(Mk,4.))*Qmpil(gasmix,2,2,i,k)) -
                        (56.*pow(Mk,4.)*Qmpil(gasmix,2,3,i,k)) + 
                        (8.*pow(Mk,4.)*Qmpil(gasmix,2,4,i,k))  +
                        (16.*(Mi*pow(Mk,3.))*Qmpil(gasmix,3,3,i,k) ))) +
                    
                    (delta(j,k) * (pow(Mi,2.)*pow(Mk,2.)) * ((-1540./3.)*Qmpil(gasmix,1,1,i,k) +
                       (784./3.)*Qmpil(gasmix,1,2,i,k) - (128./3.)*Qmpil(gasmix,1,3,i,k) +
                       (602./3.)*Qmpil(gasmix,2,2,i,k) - (56.)*Qmpil(gasmix,2,3,i,k) +
                       8.*Qmpil(gasmix,2,4,i,k)-16.*Qmpil(gasmix,3,3,i,k)))
                ) ;         
            }
            
         break;
        }
     break;
    }
    
    return sum ;
};

double ZMAppendix::qmpi1 ( GasMixture* gasmix , int m , int p , int i ) {
    
    double sum = 0.;
    double M1,Mi ; 

    double n1 = gasmix->getCompositionObj()->compositions(1.e-18).back() ;
    double ni = gasmix->getCompositionObj()->compositions(1.e-18)[i] ;

    double m1 = gasmix->masses(1.e+3).back() ;
    double mi = gasmix->masses(1.e+3)[i] ;

    int N = gasmix->getN() ;
    int e = N-1 ;

    Mi = mi / (mi+m1) ;
    M1 = m1 / (mi+m1) ;

    switch (m) {
     
     case 0:
        
        switch (p) {
     
         case 0:
             
            sum += -8.*Qmpil(gasmix,1,1,i,e) ; 
            sum *= (pow(M1,0.5)*pow(Mi,0.5)) ;
            break;
         
         case 1:
             
            sum += -20.*Qmpil(gasmix,1,1,i,e) +
                8.*Qmpil(gasmix,1,2,i,e); 
            sum *= (pow(Mi,(3./2.))*pow(M1,0.5)) ;
            break;
    
         case 2:
             
            sum += -35.*Qmpil(gasmix,1,1,i,e) +
            28.*Qmpil(gasmix,1,2,i,e) - 
            4.*Qmpil(gasmix,1,3,i,e) ; 
            sum *= (pow(Mi,(5./2.))*pow(M1,0.5)) ; 
            break;

         case 3:
            
            sum += -(105./2.)*Qmpil(gasmix,1,1,i,e) +
            63.*Qmpil(gasmix,1,2,i,e) - 
            18.*Qmpil(gasmix,1,3,i,e) +
            (4./3.)*Qmpil(gasmix,1,4,i,e) ; 
            sum *= (pow(Mi,(7./2.))*pow(M1,0.5)) ; 
            break;
         
         
        }
        break;
     
     case 1:
     
        switch (p) {
     
         case 0:
            sum += -20.*Qmpil(gasmix,1,1,i,e) +
                8.*Qmpil(gasmix,1,2,i,e); 
            sum *= (pow(Mi,(3./2.))*pow(M1,0.5)) ;
            sum *= pow(Mi/M1,-1.) ;
            break;
     
         case 1:
            
            sum += -110.*Qmpil(gasmix,1,1,i,e) +
            40.*Qmpil(gasmix,1,2,i,e) - 
            8.*Qmpil(gasmix,1,3,i,e) +
            16.*Qmpil(gasmix,2,2,i,e) ; 
            sum *= (pow(Mi,(3./2.))*pow(M1,(3./2.))) ; 
            break;

         case 2:
             
            sum += -(595./2.)*Qmpil(gasmix,1,1,i,e) +
            189.*Qmpil(gasmix,1,2,i,e) - 
            38.*Qmpil(gasmix,1,3,i,e) +
            4.*Qmpil(gasmix,1,4,i,e) +
            56.*Qmpil(gasmix,2,2,i,e) -
            16.*Qmpil(gasmix,2,3,i,e) ; 
            sum *= (pow(Mi,(5./2.))*pow(M1,(3./2.))) ; 
            break;

         case 3:
            
            sum += -(2415./4.)*Qmpil(gasmix,1,1,i,e) +
            588.*Qmpil(gasmix,1,2,i,e) - 
            162.*Qmpil(gasmix,1,3,i,e) + 
            (64./3.)*Qmpil(gasmix,1,4,i,e) - 
            (4./3.)*Qmpil(gasmix,1,5,i,e) + 
            126.*Qmpil(gasmix,2,2,i,e) - 
            72.*Qmpil(gasmix,2,3,i,e) +
            8.*Qmpil(gasmix,2,4,i,e) ; 
            sum *= (pow(Mi,(7./2.))*pow(M1,(3./2.))) ; 
            break;
         
         
        }
        break;
     
     case 2:
     
        switch (p) {
     
         case 0:
             
            sum += -35.*Qmpil(gasmix,1,1,i,e) +
            28.*Qmpil(gasmix,1,2,i,e) - 
            4.*Qmpil(gasmix,1,3,i,e) ; 
            sum *= (pow(Mi,(5./2.))*pow(M1,0.5)) ; 
            sum *= pow(Mi/M1,-2.) ;
            break;
     
         case 1:

            sum += -(595./2.)*Qmpil(gasmix,1,1,i,e) +
            189.*Qmpil(gasmix,1,2,i,e) - 
            38.*Qmpil(gasmix,1,3,i,e) +
            4.*Qmpil(gasmix,1,4,i,e) +
            56.*Qmpil(gasmix,2,2,i,e) -
            16.*Qmpil(gasmix,2,3,i,e) ; 
            sum *= (pow(Mi,(5./2.))*pow(M1,(3./2.))) ; 
            sum *= pow(Mi/M1,-1.) ;
            break;
     
         case 2:
             
            sum += -(8505./8.)*Qmpil(gasmix,1,1,i,e) +
            833.*Qmpil(gasmix,1,2,i,e) - 
            241.*Qmpil(gasmix,1,3,i,e) + 
            28.*Qmpil(gasmix,1,4,i,e) - 
            2.*Qmpil(gasmix,1,5,i,e) + 
            308.*Qmpil(gasmix,2,2,i,e) - 
            112.*Qmpil(gasmix,2,3,i,e) +
            16.*Qmpil(gasmix,2,4,i,e) - 
            16.*Qmpil(gasmix,3,3,i,e) ; 
            sum *= (pow(Mi,(5./2.))*pow(M1,(5./2.))) ; 
            break;

         case 3:
             
            sum += -(42735./16.)*Qmpil(gasmix,1,1,i,e) +
            (22071./8.)*Qmpil(gasmix,1,2,i,e) - 
            (2001./2.)*Qmpil(gasmix,1,3,i,e) + 
            (499./3.)*Qmpil(gasmix,1,4,i,e) - 
            (41./3.)*Qmpil(gasmix,1,5,i,e) + 
            (2./3.)*Qmpil(gasmix,1,6,i,e) + 
            945.*Qmpil(gasmix,2,2,i,e) -
            522.*Qmpil(gasmix,2,3,i,e) + 
            100.*Qmpil(gasmix,2,4,i,e) - 
            8.*Qmpil(gasmix,2,5,i,e) - 
            72.*Qmpil(gasmix,3,3,i,e) +
            16.*Qmpil(gasmix,3,4,i,e) ; 
            sum *= (pow(Mi,(7./2.))*pow(M1,(5./2.))) ; 
            break;
         
         
        }
        break;

     case 3:
        
        switch (p) {
        
         case 0:
        
            sum += -(105./2.)*Qmpil(gasmix,1,1,i,e) +
            63.*Qmpil(gasmix,1,2,i,e) - 
            18.*Qmpil(gasmix,1,3,i,e) +
            (4./3.)*Qmpil(gasmix,1,4,i,e) ; 
            sum *= (pow(Mi,(7./2.))*pow(M1,0.5)) ;
            sum *= pow((Mi/M1),-3.) ;
            break;
        
         case 1:

            sum += -(2415./4.)*Qmpil(gasmix,1,1,i,e) +
            588.*Qmpil(gasmix,1,2,i,e) - 
            162.*Qmpil(gasmix,1,3,i,e) + 
            (64./3.)*Qmpil(gasmix,1,4,i,e) - 
            (4./3.)*Qmpil(gasmix,1,5,i,e) + 
            126.*Qmpil(gasmix,2,2,i,e) - 
            72.*Qmpil(gasmix,2,3,i,e) +
            8.*Qmpil(gasmix,2,4,i,e) ; 
            sum *= (pow(Mi,(7./2.))*pow(M1,(3./2.))) ; 
            sum *= pow(Mi/M1,-2.) ;
            break;
         
         case 2:

            sum += -(42735./16.)*Qmpil(gasmix,1,1,i,e) +
            (22071./8.)*Qmpil(gasmix,1,2,i,e) - 
            (2001./2.)*Qmpil(gasmix,1,3,i,e) + 
            (499./3.)*Qmpil(gasmix,1,4,i,e) - 
            (41./3.)*Qmpil(gasmix,1,5,i,e) + 
            (2./3.)*Qmpil(gasmix,1,6,i,e) + 
            945.*Qmpil(gasmix,2,2,i,e) -
            522.*Qmpil(gasmix,2,3,i,e) + 
            100.*Qmpil(gasmix,2,4,i,e) - 
            8.*Qmpil(gasmix,2,5,i,e) - 
            72.*Qmpil(gasmix,3,3,i,e) +
            16.*Qmpil(gasmix,3,4,i,e) ; 
            sum *= (pow(Mi,(7./2.))*pow(M1,(5./2.))) ; 
            sum *= pow(Mi/M1,-1.) ;
            break;
         
         case 3:
             
            sum += -(255255./32.)*Qmpil(gasmix,1,1,i,e) +
            (76923./8.)*Qmpil(gasmix,1,2,i,e) - 
            (34119./8.)*Qmpil(gasmix,1,3,i,e) + 
            895.*Qmpil(gasmix,1,4,i,e) - 
            (201./2.)*Qmpil(gasmix,1,5,i,e) + 
            6.*Qmpil(gasmix,1,6,i,e) - 
            (2./9.)*Qmpil(gasmix,1,7,i,e) +
            (14553./4.)*Qmpil(gasmix,2,2,i,e) -
            2430.*Qmpil(gasmix,2,3,i,e) + 
            626.*Qmpil(gasmix,2,4,i,e) - 
            72.*Qmpil(gasmix,2,5,i,e) +
            4.*Qmpil(gasmix,2,6,i,e) - 
            444.*Qmpil(gasmix,3,3,i,e) +
            144.*Qmpil(gasmix,3,4,i,e) -
            16.*Qmpil(gasmix,3,5,i,e) +
            (32./3.)*Qmpil(gasmix,4,4,i,e) ; 
            sum *= (pow(Mi,(7./2.))*pow(M1,(7./2.))) ;
            break;
              
        }
        break;
    
    }
    return sum * ni * n1 ; 
}

double ZMAppendix::qsimpmpij ( GasMixture* gasmix , int m , int p ) {

    std::vector<double> n = gasmix->getCompositionObj()->compositions(1.e-18) ;
    std::vector<double> mass = gasmix->masses(1.e+03) ;
    int N_SPC = gasmix->getN() ;

    // e- 
    int e = N_SPC-1;
    double m1 = mass.back() ;
    double n1 = n.back() ; 


    double mj;
    double nj;
    double M1;
    double Mj;

    double sum = 0.;

    switch (m)
    {
     case 0:
        switch (p)
        {
         case 0:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    8.*Mj*Qmpil(gasmix,1,1,e,j) 
                ) ;

            }
            
            sum += 0. ;
            break;
         case 1:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    (20. * pow(Mj, 2) * Qmpil(gasmix, 1, 1, e, j)) - 
                    (8. * pow(Mj, 2) * Qmpil(gasmix, 1, 2, e, j))
                ) ;

            }
            
            sum += 0. ;
            break;
         case 2:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    (35. * pow(Mj, 3) * Qmpil(gasmix, 1, 1, e, j)) - 
                    (28. * pow(Mj, 3) * Qmpil(gasmix, 1, 2, e, j)) + 
                    (4. * pow(Mj, 3) * Qmpil(gasmix, 1, 3, e, j))

                ) ;

            }
            
            sum += 0. ;
            break;

         case 3:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    ((105. / 2.) * pow(Mj, 4) * Qmpil(gasmix, 1, 1, e, j)) - 
                    (63. * pow(Mj, 4) * Qmpil(gasmix, 1, 2, e, j)) + 
                    (18. * pow(Mj, 4) * Qmpil(gasmix, 1, 3, e, j)) - 
                    ((4. / 3.) * pow(Mj, 4) * Qmpil(gasmix, 1, 4, e, j))
                ) ;

            }
            
            sum += 0. ;
            break;
        }
        break;

     case 1:
        switch (p)
        {
         case 0:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    (20. * pow(Mj, 2) * Qmpil(gasmix, 1, 1, e, j)) - 
                    (8. * pow(Mj, 2) * Qmpil(gasmix, 1, 2, e, j))
                ) ;

            }
            
            sum += 0. ;
            break;
         case 1:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    (110. * ((6. / 11.) * pow(M1, 2) * Mj + (5. / 11.) * pow(Mj, 3)) * Qmpil(gasmix, 1, 1, e, j)) - 
                    (40. * pow(Mj, 3) * Qmpil(gasmix, 1, 2, e, j)) + 
                    (8. * pow(Mj, 3) * Qmpil(gasmix, 1, 3, e, j)) 
                ) ;

            }
            
            sum += pow(n1,2.) * (4.*Qmpil(gasmix,2,2,e,e)) ;
            break;
         case 2:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    ((595. / 2.) * ((12. / 17.) * pow(M1, 2) * pow(Mj, 2) + (5. / 17.) * pow(Mj, 4)) * Qmpil(gasmix, 1, 1, e, j)) -                         
                    (189. * ((4. / 9.) * pow(M1, 2) * pow(Mj, 2) + (5. / 9.) * pow(Mj, 4)) * Qmpil(gasmix, 1, 2, e, j)) +                         
                    (38. * pow(Mj, 4) * Qmpil(gasmix, 1, 3, e, j)) -                         
                    (4. * pow(Mj, 4) * Qmpil(gasmix, 1, 4, e, j)) 
                ) ;

            }
            
            sum += pow(n1,2.) * (7.*Qmpil(gasmix,2,2,e,e)-2.*Qmpil(gasmix,2,3,e,e)) ;
            break;
         case 3:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    ((2415. / 4.) * ((18. / 23.) * pow(M1, 2) * pow(Mj, 3) + (5. / 23.) * pow(Mj, 5)) * Qmpil(gasmix, 1, 1, e, j)) - 
                    (588. * ((9. / 14.) * pow(M1, 2) * pow(Mj, 3) + (5. / 14.) * pow(Mj, 5)) * Qmpil(gasmix, 1, 2, e, j)) + 
                    (162. * ((1. / 3.) * pow(M1, 2) * pow(Mj, 3) + (2. / 3.) * pow(Mj, 5)) * Qmpil(gasmix, 1, 3, e, j)) - 
                    ((64. / 3.) * pow(Mj, 5) * Qmpil(gasmix, 1, 4, e, j)) + 
                    ((4. / 3.) * pow(Mj, 5) * Qmpil(gasmix, 1, 5, e, j)) 
                ) ;

            }
            
            sum += pow(n1,2.) * (
                (63./8.) * Qmpil(gasmix,2,2,e,e) -
                (9./2.) * Qmpil(gasmix,2,3,e,e) +
                (1./2.) * Qmpil(gasmix,2,4,e,e) 
            ) ;
            break;
        }
        break;

     case 2:
        switch (p)
        {
         case 0:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    (35. * pow(Mj, 3) * Qmpil(gasmix, 1, 1, e, j)) - 
                    (28. * pow(Mj, 3) * Qmpil(gasmix, 1, 2, e, j)) + 
                    (4. * pow(Mj, 3) * Qmpil(gasmix, 1, 3, e, j))
                ) ;

            }
            
            sum += 0. ;
            break;
         case 1:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    ((595. / 2.) * ((12. / 17.) * pow(M1, 2) * pow(Mj, 2) + (5. / 17.) * pow(Mj, 4)) * Qmpil(gasmix, 1, 1, e, j)) -                         
                    (189. * ((4. / 9.) * pow(M1, 2) * pow(Mj, 2) + (5. / 9.) * pow(Mj, 4)) * Qmpil(gasmix, 1, 2, e, j)) +                         
                    (38. * pow(Mj, 4) * Qmpil(gasmix, 1, 3, e, j)) -                         
                    (4. * pow(Mj, 4) * Qmpil(gasmix, 1, 4, e, j)) 
                ) ;

            }
            
            sum += pow(n1,2.) * (7.*Qmpil(gasmix,2,2,e,e)-2.*Qmpil(gasmix,2,3,e,e)) ;
            break;
         case 2:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    ((8505. / 8.) * ((40. / 243.) * pow(M1, 4) * Mj + (168. / 243.) * pow(M1, 2) * pow(Mj, 3) + (35. / 243.) * pow(Mj, 5)) * Qmpil(gasmix, 1, 1, e, j)) - 
                    (833. * ((12. / 17.) * pow(M1, 2) * pow(Mj, 3) + (5. / 17.) * pow(Mj, 5)) * Qmpil(gasmix, 1, 2, e, j)) + 
                    (241. * ((108. / 241.) * pow(M1, 2) * pow(Mj, 3) + (133. / 241.) * pow(Mj, 5)) * Qmpil(gasmix, 1, 3, e, j)) - 
                    (28. * pow(Mj, 5) * Qmpil(gasmix, 1, 4, e, j)) + 
                    (2. * pow(Mj, 5) * Qmpil(gasmix, 1, 5, e, j)) 
                ) ;

            }
            
            sum += pow(n1,2.) * (
                (77./4.)*Qmpil(gasmix,2,2,e,e)-7.*Qmpil(gasmix,2,3,e,e)+Qmpil(gasmix,2,4,e,e)
            ) ;
            break;
         case 3:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    ((42735. / 16.) * ((120. / 407.) * pow(M1, 4) * pow(Mj, 2) + (252. / 407.) * pow(M1, 2) * pow(Mj, 4) + (35. / 407.) * pow(Mj, 6)) * Qmpil(gasmix, 1, 1, e, j)) - 
                    ((22071. / 8.) * ((120. / 1051.) * pow(M1, 4) * pow(Mj, 2) + (756. / 1051.) * pow(M1, 2.) * pow(Mj, 4.) + (175. / 1051.) * pow(Mj, 6.)) * Qmpil(gasmix, 1, 2, e, j)) + 
                    ((2001. / 2.) * ((450. / 667.) * pow(M1, 2) * pow(Mj, 4) + (217. / 667.) * pow(Mj, 6)) * Qmpil(gasmix, 1, 3, e, j)) - 
                    ((499. / 3.) * ((198. / 499.) * pow(M1, 2) * pow(Mj, 4) + (301. / 499.) * pow(Mj, 6)) * Qmpil(gasmix, 1, 4, e, j)) + 
                    ((41. / 3.) * pow(Mj, 6) * Qmpil(gasmix, 1, 5, e, j)) - 
                    ((2. / 3.) * pow(Mj, 6) * Qmpil(gasmix, 1, 6, e, j)) 

                ) ;

            }
            
            sum += pow(n1,2.) * (
                (945./32.)*Qmpil(gasmix,2,2,e,e)-(261./16.)*Qmpil(gasmix,2,3,e,e) +
                (25./8.)*Qmpil(gasmix,2,4,e,e)-(1./4.)*Qmpil(gasmix,2,5,e,e)
            ) ;
            break;
        }
        break;

     case 3:
        switch (p)
        {
         case 0:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    ((105. / 2.) * pow(Mj, 4) * Qmpil(gasmix, 1, 1, e, j)) - 
                    (63. * pow(Mj, 4) * Qmpil(gasmix, 1, 2, e, j)) + 
                    (18. * pow(Mj, 4) * Qmpil(gasmix, 1, 3, e, j)) - 
                    ((4. / 3.) * pow(Mj, 4) * Qmpil(gasmix, 1, 4, e, j))
                ) ;

            }
            
            sum += 0. ;
            break;
         case 1:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    ((2415. / 4.) * ((18. / 23.) * pow(M1, 2) * pow(Mj, 3) + (5. / 23.) * pow(Mj, 5)) * Qmpil(gasmix, 1, 1, e, j)) - 
                    (588. * ((9. / 14.) * pow(M1, 2) * pow(Mj, 3) + (5. / 14.) * pow(Mj, 5)) * Qmpil(gasmix, 1, 2, e, j)) + 
                    (162. * ((1. / 3.) * pow(M1, 2) * pow(Mj, 3) + (2. / 3.) * pow(Mj, 5)) * Qmpil(gasmix, 1, 3, e, j)) - 
                    ((64. / 3.) * pow(Mj, 5) * Qmpil(gasmix, 1, 4, e, j)) + 
                    ((4. / 3.) * pow(Mj, 5) * Qmpil(gasmix, 1, 5, e, j)) 
                ) ;

            }
            
            sum += pow(n1,2.) * (
                (63./8.) * Qmpil(gasmix,2,2,e,e) -
                (9./2.) * Qmpil(gasmix,2,3,e,e) +
                (1./2.) * Qmpil(gasmix,2,4,e,e) ) ;
            break;
         case 2:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    ((42735. / 16.) * ((120. / 407.) * pow(M1, 4) * pow(Mj, 2) + (252. / 407.) * pow(M1, 2) * pow(Mj, 4) + (35. / 407.) * pow(Mj, 6)) * Qmpil(gasmix, 1, 1, e, j)) - 
                    ((22071. / 8.) * ((120. / 1051.) * pow(M1, 4) * pow(Mj, 2) + (756. / 1051.) * pow(M1, 2.) * pow(Mj, 4.) + (175. / 1051.) * pow(Mj, 6.)) * Qmpil(gasmix, 1, 2, e, j)) + 
                    ((2001. / 2.) * ((450. / 667.) * pow(M1, 2) * pow(Mj, 4) + (217. / 667.) * pow(Mj, 6)) * Qmpil(gasmix, 1, 3, e, j)) - 
                    ((499. / 3.) * ((198. / 499.) * pow(M1, 2) * pow(Mj, 4) + (301. / 499.) * pow(Mj, 6)) * Qmpil(gasmix, 1, 4, e, j)) + 
                    ((41. / 3.) * pow(Mj, 6) * Qmpil(gasmix, 1, 5, e, j)) - 
                    ((2. / 3.) * pow(Mj, 6) * Qmpil(gasmix, 1, 6, e, j))
                ) ;

            }
            
            sum += pow(n1,2.) * (
                (945./32.)*Qmpil(gasmix,2,2,e,e)-(261./16.)*Qmpil(gasmix,2,3,e,e) +
                (25./8.)*Qmpil(gasmix,2,4,e,e)-(1./4.)*Qmpil(gasmix,2,5,e,e)
            ) ;
            break;
         case 3:
            for (int j  = 0; j < N_SPC - 1 ; j++)
            {
                nj = n[j] ; 
                mj = mass[j] ; 
                M1 = m1 / (m1+mj) ; 
                Mj = mj / (m1+mj) ; 

                sum += n1 * nj * (
                    ((255255. / 32.) * ((112. / 2431.) * pow(M1, 6) * Mj + (1080. / 2431.) * pow(M1, 4) * pow(Mj, 3) + (1134. / 2431.) * pow(M1, 2) * pow(Mj, 5) + (105. / 2431.) * pow(Mj, 7)) * Qmpil(gasmix, 1, 1, e, j)) - 
                    ((76923. / 8.) * ((120. / 407.) * pow(M1, 4) * pow(Mj, 3) + (252. / 407.) * pow(M1, 2) * pow(Mj, 5) + (35. / 407.) * pow(Mj, 7)) * Qmpil(gasmix, 1, 2, e, j)) + 
                    ((34119./8.)*((440./3791.)*pow(M1,4.)*pow(Mj,3.) + (2700. / 3791.) * pow(M1, 2) * pow(Mj, 5) + (651. / 3791.) * pow(Mj, 7)) * Qmpil(gasmix, 1, 3, e, j)) - 
                    ((895.) * ((594. / 895.) * pow(M1, 2) * pow(Mj, 5) + (301. / 895.) * pow(Mj, 7)) * Qmpil(gasmix, 1, 4, e, j)) + 
                    ((201. / 2.) * ((26. / 67.) * pow(M1, 2) * pow(Mj, 5) + (41. / 67.) * pow(Mj, 7)) * Qmpil(gasmix, 1, 5, e, j)) - 
                    (6. * pow(Mj, 7) * Qmpil(gasmix, 1, 6, e, j)) + 
                    ((2. / 9.) * pow(Mj, 7) * Qmpil(gasmix, 1, 7, e, j)) 
                ) ;

            }
            
            sum += pow(n1,2.) * (
                ((14553. / 256.) * Qmpil(gasmix, 2, 2, e, e)) - 
                ((1215. / 32.) * Qmpil(gasmix, 2, 3, e, e)) + 
                ((313. / 32.) * Qmpil(gasmix, 2, 4, e, e)) - 
                ((9. / 8.) * Qmpil(gasmix, 2, 5, e, e)) + 
                ((1. / 16.) * Qmpil(gasmix, 2, 6, e, e)) + 
                ((1. / 6.) * Qmpil(gasmix, 4, 4, e, e))
            ) ;
            break;
        }
        break;
    }
    return sum ;
}


double ZMAppendix::qmpij ( GasMixture* gasmix , int m , int p , int i , int j ) {
    
    std::vector<double> n = gasmix->getCompositionObj()->compositions(1.e-18) ; 
    std::vector<double> mass = gasmix->masses(1.e+03) ;
    int N = gasmix->getN() ;

    double ni = n[i] ; 
    double nk ;
    double sum = 0.;

    double mi = mass[i];
    double mk ;
    double Mi ; 
    double Mk ;
    
    switch (m) {

     case 0:
        
        switch (p) {
         
         case 0:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        8.*Mk*Qmpil(gasmix,1,1,i,k) 
                    ) ) + 
                    ( delta(j,k) * (pow(Mi,0.5)*pow(Mk,0.5)) * (
                        -8.*Qmpil(gasmix,1,1,i,k)
                    ) )
                ); 
                 
            }
            
         break;
         
         case 1:
            for (int k = 0; k < N - 1  ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        20.*pow(Mk,2.)*Qmpil(gasmix,1,1,i,k) -
                        8.*pow(Mk,2.)*Qmpil(gasmix,1,2,i,k)
                    ) ) + 
                    ( delta(j,k) * (pow(Mi,3./2.)*pow(Mk,0.5)) * (
                        -20.*Qmpil(gasmix,1,1,i,k) +
                        8.*Qmpil(gasmix,1,2,i,k)
                    ) )
                );   
            }
         break;
         
         case 2:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        35.*pow(Mk,3)*Qmpil(gasmix,1,1,i,k) -
                        28.*pow(Mk,3)*Qmpil(gasmix,1,2,i,k) +
                         4.*pow(Mk,3)*Qmpil(gasmix,1,3,i,k)
                    ) ) + 
                    ( delta(j,k) * (pow(Mi,5./2.)*pow(Mk,0.5)) * (
                        -35.*Qmpil(gasmix,1,1,i,k) +
                        28.*Qmpil(gasmix,1,2,i,k) -
                         4.*Qmpil(gasmix,1,3,i,k)
                    ) )
                );   
            }
         break;

         case 3:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        (105./2.)*pow(Mk,4.)*Qmpil(gasmix,1,1,i,k) - 
                        (63.)*pow(Mk,4.)*Qmpil(gasmix,1,2,i,k) +
                        (18.)*pow(Mk,4.)*Qmpil(gasmix,1,3,i,k) -
                        (4./3.)*pow(Mk,4.)*Qmpil(gasmix,1,4,i,k)
                    ) ) + 
                    ( delta(j,k) * (pow(Mi,7./2.)*pow(Mk,0.5)) * (
                        - (105./2.)*Qmpil(gasmix,1,1,i,k) +
                        (63.)*Qmpil(gasmix,1,2,i,k) -
                        (18.)*Qmpil(gasmix,1,3,i,k) +
                        (4./3.)*Qmpil(gasmix,1,4,i,k)
                    ) )
                );   
            }
         break;
        }
     break;

     case 1:
        
        switch (p) {
         
         case 0:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        20.*pow(Mk,2.)*Qmpil(gasmix,1,1,i,k) -
                        8.*pow(Mk,2.)*Qmpil(gasmix,1,2,i,k)
                    ) ) + 
                    ( delta(j,k) * (pow(Mi/Mk,-1)) * (pow(Mi,3./2.)*pow(Mk,0.5)) * (
                        -20.*Qmpil(gasmix,1,1,i,k) +
                        8.*Qmpil(gasmix,1,2,i,k) 
                    ) )
                ); 
                 
            }
            
         break;
         
         case 1:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        (110.*((6./11.)*pow(Mi,2.)*Mk+(5./11.)*pow(Mk,3.))*Qmpil(gasmix,1,1,i,k)) -
                        (40.*pow(Mk,3.)*Qmpil(gasmix,1,2,i,k)) +
                        (8.*pow(Mk,3.)*Qmpil(gasmix,1,3,i,k)) +
                        (16.*Mi*pow(Mk,2.)*Qmpil(gasmix,2,2,i,k)) 
                    ) ) + 
                    ( delta(j,k) * (pow(Mi,3./2.)*pow(Mk,3./2.)) * (
                        (-110. * Qmpil(gasmix, 1, 1, i, k)) + 
                        (40. * Qmpil(gasmix, 1, 2, i, k)) - 
                        (8. * Qmpil(gasmix, 1, 3, i, k)) + 
                        (16. * Qmpil(gasmix, 2, 2, i, k))
                    ) )
                );   
            }
         break;
         
         case 2:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        ((595./2.)*((12./17.)*pow(Mi,2.)*pow(Mk,2.)+(5./17.)*pow(Mk,4.))*Qmpil(gasmix,1,1,i,k)) - 
                        ((189.)*((4./9.)*pow(Mi,2.)*pow(Mk,2.)+(5./9.)*pow(Mk,4.))*Qmpil(gasmix,1,2,i,k)) +
                        ((38.)*pow(Mk,4.)*Qmpil(gasmix,1,3,i,k)) - 
                        ((4.)*pow(Mk,4.)*Qmpil(gasmix,1,4,i,k)) + 
                        ((56.)*Mi*pow(Mk,3.)*Qmpil(gasmix,2,2,i,k)) - 
                        ((16.)*Mi*pow(Mk,3.)*Qmpil(gasmix,2,3,i,k))
                    ) ) + 
                    ( delta(j,k) * (pow(Mi,5./2.)*pow(Mk,3./2.)) * (
                        (-(595. / 2.) * Qmpil(gasmix, 1, 1, i, k)) + 
                        (189. * Qmpil(gasmix, 1, 2, i, k)) - 
                        (38. * Qmpil(gasmix, 1, 3, i, k)) + 
                        (4. * Qmpil(gasmix, 1, 4, i, k)) + 
                        (56. * Qmpil(gasmix, 2, 2, i, k)) - 
                        (16. * Qmpil(gasmix, 2, 3, i, k))

                    ) )
                );   
            }
         break;

         case 3:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        ((2415./4.)*((18./23.)*pow(Mi,2.)*pow(Mk,3.)+(5./23.)*pow(Mk,5.))*Qmpil(gasmix,1,1,i,k)) - 
                        ((588.)*((9./14.)*pow(Mi,2.)*pow(Mk,3.)+(5./14.)*pow(Mk,5.))*Qmpil(gasmix,1,2,i,k)) +
                        ((162.)*((1./3.)*pow(Mi,2.)*pow(Mk,3.)+(2./3.)*pow(Mk,5.))*Qmpil(gasmix,1,3,i,k)) - 
                        ((64./3.)*pow(Mk,5.)*Qmpil(gasmix,1,4,i,k)) + 
                        ((4./3.)*pow(Mk,5.)*Qmpil(gasmix,1,5,i,k)) + 
                        ((126.)*Mi*pow(Mk,4.)*Qmpil(gasmix,2,2,i,k)) - 
                        ((72.)*Mi*pow(Mk,4.)*Qmpil(gasmix,2,3,i,k)) +
                        ((8.)*Mi*pow(Mk,4.)*Qmpil(gasmix,2,4,i,k))
                    ) ) + 
                    ( delta(j,k) * (pow(Mi,7./2.)*pow(Mk,3./2.)) * (
                        (-(2415. / 4.) * Qmpil(gasmix, 1, 1, i, k)) + 
                        (588. * Qmpil(gasmix, 1, 2, i, k)) - 
                        (162. * Qmpil(gasmix, 1, 3, i, k)) + 
                        ((64. / 3.) * Qmpil(gasmix, 1, 4, i, k)) - 
                        ((4. / 3. )* Qmpil(gasmix, 1, 5, i, k)) + 
                        (126. * Qmpil(gasmix, 2, 2, i, k)) - 
                        (72. * Qmpil(gasmix, 2, 3, i, k)) + 
                        (8. * Qmpil(gasmix, 2, 4, i, k))

                    ) )
                );   
            }
         break;
        }
     break;

     case 2:
        
        switch (p) {
         
         case 0:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        35.*pow(Mk,3)*Qmpil(gasmix,1,1,i,k) -
                        28.*pow(Mk,3)*Qmpil(gasmix,1,2,i,k) +
                         4.*pow(Mk,3)*Qmpil(gasmix,1,3,i,k)
                    ) ) + 
                    ( delta(j,k) * pow(Mi/Mk,-2.) * (pow(Mi,5./2.)*pow(Mk,0.5)) * (
                        -35.*Qmpil(gasmix,1,1,i,k) +
                        28.*Qmpil(gasmix,1,2,i,k) -
                         4.*Qmpil(gasmix,1,3,i,k) 
                    ))
                ); 
                 
            }
            
         break;
         
         case 1:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        ((595./2.)*((12./17.)*pow(Mi,2.)*pow(Mk,2.)+(5./17.)*pow(Mk,4.))*Qmpil(gasmix,1,1,i,k)) - 
                        ((189.)*((4./9.)*pow(Mi,2.)*pow(Mk,2.)+(5./9.)*pow(Mk,4.))*Qmpil(gasmix,1,2,i,k)) +
                        ((38.)*pow(Mk,4.)*Qmpil(gasmix,1,3,i,k)) - 
                        ((4.)*pow(Mk,4.)*Qmpil(gasmix,1,4,i,k)) + 
                        ((56.)*Mi*pow(Mk,3.)*Qmpil(gasmix,2,2,i,k)) - 
                        ((16.)*Mi*pow(Mk,3.)*Qmpil(gasmix,2,3,i,k))
                    ) ) + 
                    ( delta(j,k) * pow(Mi/Mk,-1.) * (pow(Mi,5./2.)*pow(Mk,3./2.)) * (
                        (-(595. / 2.) * Qmpil(gasmix, 1, 1, i, k)) + 
                        (189. * Qmpil(gasmix, 1, 2, i, k)) - 
                        (38. * Qmpil(gasmix, 1, 3, i, k)) + 
                        (4. * Qmpil(gasmix, 1, 4, i, k)) + 
                        (56. * Qmpil(gasmix, 2, 2, i, k)) - 
                        (16. * Qmpil(gasmix, 2, 3, i, k))

                    ) )
                );   
            }
         break;
         
         case 2:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        ((8505./8.)*((40./243.)*pow(Mi,4.)*Mk+(168./243.)*pow(Mi,2.)*pow(Mk,3.)+(35./243.)*pow(Mk,5.))*Qmpil(gasmix,1,1,i,k)) -
                        ((833.)*((12./17.)*pow(Mi,2.)*pow(Mk,3.)+(5./17.)*pow(Mk,5.) )*Qmpil(gasmix,1,2,i,k)) +
                        ((241.)*((108./241.)*pow(Mi,2.)*pow(Mk,3.)+(133./241.)*pow(Mk,5.) )*Qmpil(gasmix,1,3,i,k)) -
                        ((28.)*pow(Mk,5.)*Qmpil(gasmix,1,4,i,k)) +
                        ((2.)*pow(Mk,5.)*Qmpil(gasmix,1,5,i,k)) +
                        ((308.)*((4./11.)*pow(Mi,3.)*pow(Mk,2.)+(7./11.)*Mi*pow(Mk,4.) )*Qmpil(gasmix,2,2,i,k)) -
                        ((112.)*Mi*pow(Mk,4.)*Qmpil(gasmix,2,3,i,k)) +
                        ((16.)*Mi*pow(Mk,4.)*Qmpil(gasmix,2,4,i,k)) +
                        ((16.)*pow(Mi,2.)*pow(Mk,3.)*Qmpil(gasmix,3,3,i,k)) 
                    ) ) + 
                    ( delta(j,k) * (pow(Mi,5./2.)*pow(Mk,5./2.)) * (
                        ((-8505. / 8.) * Qmpil(gasmix, 1, 1, i, k)) + 
                        (833. * Qmpil(gasmix, 1, 2, i, k)) - 
                        (241. * Qmpil(gasmix, 1, 3, i, k)) + 
                        (28. * Qmpil(gasmix, 1, 4, i, k)) - 
                        (2. * Qmpil(gasmix, 1, 5, i, k)) + 
                        (308. * Qmpil(gasmix, 2, 2, i, k)) - 
                        (112. * Qmpil(gasmix, 2, 3, i, k)) + 
                        (16. * Qmpil(gasmix, 2, 4, i, k)) - 
                        (16. * Qmpil(gasmix, 3, 3, i, k))
                    ) )
                );   
            }
         break;

         case 3:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        ((42735./16.)*((120./407.)*pow(Mi,4.)*pow(Mk,2.)+(252./407.)*pow(Mi,2.)*pow(Mk,4.)+(35./407.)*pow(Mk,6.) )*Qmpil(gasmix,1,1,i,k)) -
                        ((22071./8.)*((120./1051.)*pow(Mi,4.)*pow(Mk,2.)+(756./1051.)*pow(Mi,2)*pow(Mk,4.)+(175./1051.)*pow(Mk,6.) )*Qmpil(gasmix,1,2,i,k)) +
                        ((2001./2.)*((450./667.)*pow(Mi,2.)*pow(Mk,4.)+(217./667.)*pow(Mk,6.) )*Qmpil(gasmix,1,3,i,k)) -
                        ((499./3.)*((198/499)*pow(Mi,2.)*pow(Mk,4.)+(301./499)*pow(Mk,6.))*Qmpil(gasmix,1,4,i,k)) +
                        ((41./3.)*pow(Mk,6.)*Qmpil(gasmix,1,5,i,k)) -
                        ((2./3.)*pow(Mk,6.)*Qmpil(gasmix,1,6,i,k)) +
                        ((945.)*((8./15.)*pow(Mi,3.)*pow(Mk,3.)+(7./15.)*Mi*pow(Mk,5.))*Qmpil(gasmix,2,2,i,k)) -
                        ((522.)*((8./29.)*pow(Mi,3.)*pow(Mk,3.)+(21./29.)*Mi*pow(Mk,5.)) *Qmpil(gasmix,2,3,i,k)) +
                        ((100.)*pow(Mi,1.)*pow(Mk,5.)*Qmpil(gasmix,2,4,i,k)) -
                        ((8.)*pow(Mi,1.)*pow(Mk,5.)*Qmpil(gasmix,2,5,i,k)) +
                        ((72.)*pow(Mi,2.)*pow(Mk,4.)*Qmpil(gasmix,3,3,i,k)) -
                        ((16.)*pow(Mi,2.)*pow(Mk,4.)*Qmpil(gasmix,3,4,i,k)) 
                    ) ) + 
                    ( delta(j,k) * (pow(Mi,7./2.)*pow(Mk,5./2.)) * (
                        ((-42735. / 16.) * Qmpil(gasmix, 1, 1, i, k)) + 
                        ((22071. / 8.) * Qmpil(gasmix, 1, 2, i, k)) - 
                        ((2001. / 2.) * Qmpil(gasmix, 1, 3, i, k)) + 
                        ((499. / 3.) * Qmpil(gasmix, 1, 4, i, k)) - 
                        ((41. / 3.) * Qmpil(gasmix, 1, 5, i, k)) + 
                        ((2. / 3.) * Qmpil(gasmix, 1, 6, i, k)) + 
                        (945. * Qmpil(gasmix, 2, 2, i, k)) - 
                        (522. * Qmpil(gasmix, 2, 3, i, k)) + 
                        (100. * Qmpil(gasmix, 2, 4, i, k)) - 
                        (8. * Qmpil(gasmix, 2, 5, i, k)) - 
                        (72. * Qmpil(gasmix, 3, 3, i, k)) + 
                        (16. * Qmpil(gasmix, 3, 4, i, k))

                    ) )
                );   
            }
         break;
        }
     break;
    
     case 3:
        
        switch (p) {
         
         case 0:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        (105./2.)*pow(Mk,4.)*Qmpil(gasmix,1,1,i,k) - 
                        (63.)*pow(Mk,4.)*Qmpil(gasmix,1,2,i,k) +
                        (18.)*pow(Mk,4.)*Qmpil(gasmix,1,3,i,k) -
                        (4./3.)*pow(Mk,4.)*Qmpil(gasmix,1,4,i,k)
                    ) ) + 
                    ( delta(j,k) * pow(Mi/Mk,-3.) * (pow(Mi,7./2.)*pow(Mk,0.5)) * (
                        - (105./2.)*Qmpil(gasmix,1,1,i,k) +
                        (63.)*Qmpil(gasmix,1,2,i,k) -
                        (18.)*Qmpil(gasmix,1,3,i,k) +
                        (4./3.)*Qmpil(gasmix,1,4,i,k)
                    ) )
                ); 
                 
            }
            
         break;
         
         case 1:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        ((2415./4.)*((18./23.)*pow(Mi,2.)*pow(Mk,3.)+(5./23.)*pow(Mk,5.) )*Qmpil(gasmix,1,1,i,k)) - 
                        ((588.)*((9./14.)*pow(Mi,2.)*pow(Mk,3.)+(5./14.)*pow(Mk,5.))*Qmpil(gasmix,1,2,i,k)) +
                        ((162.)*((1./3.)*pow(Mi,2.)*pow(Mk,3.)+(2./3.)*pow(Mk,5.))*Qmpil(gasmix,1,3,i,k)) - 
                        ((64./3.)*pow(Mk,5.)*Qmpil(gasmix,1,4,i,k)) + 
                        ((4./3.)*pow(Mk,5.)*Qmpil(gasmix,1,5,i,k)) + 
                        ((126.)*Mi*pow(Mk,4.)*Qmpil(gasmix,2,2,i,k)) - 
                        ((72.)*Mi*pow(Mk,4.)*Qmpil(gasmix,2,3,i,k)) +
                        ((8.)*Mi*pow(Mk,4.)*Qmpil(gasmix,2,4,i,k))
                    ) ) + 
                    ( delta(j,k) * pow(Mi/Mk,-2.) * (pow(Mi,7./2.)*pow(Mk,3./2.)) * (
                        (-(2415. / 4.) * Qmpil(gasmix, 1, 1, i, k)) + 
                        (588. * Qmpil(gasmix, 1, 2, i, k)) - 
                        (162. * Qmpil(gasmix, 1, 3, i, k)) + 
                        ((64. / 3.) * Qmpil(gasmix, 1, 4, i, k)) - 
                        ((4. / 3. )* Qmpil(gasmix, 1, 5, i, k)) + 
                        (126. * Qmpil(gasmix, 2, 2, i, k)) - 
                        (72. * Qmpil(gasmix, 2, 3, i, k)) + 
                        (8. * Qmpil(gasmix, 2, 4, i, k))
                    ))
                );   
            }
         break;
         
         case 2:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        ((42735./16.)*((120./407.)*pow(Mi,4.)*pow(Mk,2.)+(252./407.)*pow(Mi,2.)*pow(Mk,4.)+(35./407.)*pow(Mk,6.) )*Qmpil(gasmix,1,1,i,k)) -
                        ((22071./8.)*((120./1051.)*pow(Mi,4.)*pow(Mk,2.)+(756./1051.)*pow(Mi,2)*pow(Mk,4.)+(175./1051.)*pow(Mk,6.) )*Qmpil(gasmix,1,2,i,k)) +
                        ((2001./2.)*((450./667.)*pow(Mi,2.)*pow(Mk,4.)+(217./667.)*pow(Mk,6.) )*Qmpil(gasmix,1,3,i,k)) -
                        ((499./3.)*((198/499)*pow(Mi,2.)*pow(Mk,4.)+(301./499)*pow(Mk,6.))*Qmpil(gasmix,1,4,i,k)) +
                        ((41./3.)*pow(Mk,6.)*Qmpil(gasmix,1,5,i,k)) -
                        ((2./3.)*pow(Mk,6.)*Qmpil(gasmix,1,6,i,k)) +
                        ((945.)*((8./15.)*pow(Mi,3.)*pow(Mk,3.)+(7./15.)*Mi*pow(Mk,5.))*Qmpil(gasmix,2,2,i,k)) -
                        ((522.)*((8./29.)*pow(Mi,3.)*pow(Mk,3.)+(21./29.)*Mi*pow(Mk,5.)) *Qmpil(gasmix,2,3,i,k)) +
                        ((100.)*pow(Mi,1.)*pow(Mk,5.)*Qmpil(gasmix,2,4,i,k)) -
                        ((8.)*pow(Mi,1.)*pow(Mk,5.)*Qmpil(gasmix,2,5,i,k)) +
                        ((72.)*pow(Mi,2.)*pow(Mk,4.)*Qmpil(gasmix,3,3,i,k)) -
                        ((16.)*pow(Mi,2.)*pow(Mk,4.)*Qmpil(gasmix,3,4,i,k)) 
                    ) ) + 
                    ( delta(j,k) * pow(Mi/Mk,-1.) * (pow(Mi,7./2.)*pow(Mk,5./2.)) * (
                        ((-42735. / 16.) * Qmpil(gasmix, 1, 1, i, k)) + 
                        ((22071. / 8.) * Qmpil(gasmix, 1, 2, i, k)) - 
                        ((2001. / 2.) * Qmpil(gasmix, 1, 3, i, k)) + 
                        ((499. / 3.) * Qmpil(gasmix, 1, 4, i, k)) - 
                        ((41. / 3.) * Qmpil(gasmix, 1, 5, i, k)) + 
                        ((2. / 3.) * Qmpil(gasmix, 1, 6, i, k)) + 
                        (945. * Qmpil(gasmix, 2, 2, i, k)) - 
                        (522. * Qmpil(gasmix, 2, 3, i, k)) + 
                        (100. * Qmpil(gasmix, 2, 4, i, k)) - 
                        (8. * Qmpil(gasmix, 2, 5, i, k)) - 
                        (72. * Qmpil(gasmix, 3, 3, i, k)) + 
                        (16. * Qmpil(gasmix, 3, 4, i, k))
                    ))
                );   
            }
         break;

         case 3:
            for (int k = 0; k < N - 1 ; k++) {
                
                nk = n[k];
                mk = mass[k] ;
                Mi = mi/(mi+mk) ;
                Mk = mk/(mi+mk) ;  

                sum += ni* nk * (
                    ( delta(i,j) * (
                        ((255255./32.)*((112./2431.)*pow(Mi,6.)*pow(Mk,1.)+(1080./2431.)*pow(Mi,4.)*pow(Mk,3.)+(1134./2431.)*pow(Mi,2.)*pow(Mk,5.)+(105./2431.)*pow(Mi,0.)*pow(Mk,7.))*Qmpil(gasmix,1,1,i,k)) -
                        ((76923./8.)*((120./407.)*pow(Mi,4.)*pow(Mk,3.)+(252./407.)*pow(Mi,2)*pow(Mk,5.)+(35./407.)*pow(Mk,7.))*Qmpil(gasmix,1,2,i,k)) +
                        ((34119./8.)*((440./3791.)*pow(Mi,4.)*pow(Mk,3.)+(2700./3791.)*pow(Mi,2)*pow(Mk,5.)+(651./3791.)*pow(Mk,7.))*Qmpil(gasmix,1,3,i,k)) -
                        ((895.)*((594./895.)*pow(Mi,2.)*pow(Mk,5.)+(301./895.)*pow(Mk,7.))*Qmpil(gasmix,1,4,i,k)) +
                        ((201./2.)*((26./67.)*pow(Mi,2.)*pow(Mk,5.)+(41./67.)*pow(Mk,7.))*Qmpil(gasmix,1,5,i,k)) -
                        (6.*pow(Mk,7.)*Qmpil(gasmix,1,6,i,k)) +
                        ((2./9.)*pow(Mk,7.)*Qmpil(gasmix,1,7,i,k)) +
                        ((14553./4.)*((8./77.)*pow(Mi,5.)*pow(Mk,2.)+(48./77.)*pow(Mi,3.)*pow(Mk,4.)+(21./77.)*Mi*pow(Mk,6.))*Qmpil(gasmix,2,2,i,k)) -
                        ((2430.)*((8./15.)*pow(Mi,3.)*pow(Mk,4.)+(7./15.)*Mi*pow(Mk,6.))*Qmpil(gasmix,2,3,i,k)) +
                        ((626.)*((176./626.)*pow(Mi,3.)*pow(Mk,4.)+(450./626.)*Mi*pow(Mk,6.))*Qmpil(gasmix,2,4,i,k)) -
                        (72.*Mi*pow(Mk,6.)*Qmpil(gasmix,2,5,i,k)) +
                        (4.*Mi*pow(Mk,6.)*Qmpil(gasmix,2,6,i,k)) +
                        (444.*((120./444.)*pow(Mi,4.)*pow(Mk,3.)+(324./444.)*pow(Mi,2.)*pow(Mk,5.))*Qmpil(gasmix,3,3,i,k)) -
                        (144.*pow(Mi,2.)*pow(Mk,5.)*Qmpil(gasmix,3,4,i,k)) +
                        (16.*pow(Mi,2.)*pow(Mk,5.)*Qmpil(gasmix,3,5,i,k)) +
                        ((32./3.)*pow(Mi,3.)*pow(Mk,4.)*Qmpil(gasmix,4,4,i,k))
                    ) ) + 
                    ( delta(j,k) * (pow(Mi,7./2.)*pow(Mk,7./2.)) * (
                        ((-255255. / 32.) * Qmpil(gasmix, 1, 1, i, k)) + 
                        ((76923. / 8.) * Qmpil(gasmix, 1, 2, i, k)) - 
                        ((34119. / 8.) * Qmpil(gasmix, 1, 3, i, k)) + 
                        (895. * Qmpil(gasmix, 1, 4, i, k)) - 
                        ((201. / 2.) * Qmpil(gasmix, 1, 5, i, k)) + 
                        (6. * Qmpil(gasmix, 1, 6, i, k)) - 
                        ((2. / 9.) * Qmpil(gasmix, 1, 7, i, k)) + 
                        ((14553. / 4.) * Qmpil(gasmix, 2, 2, i, k)) - 
                        (2430. * Qmpil(gasmix, 2, 3, i, k)) + 
                        (626. * Qmpil(gasmix, 2, 4, i, k)) - 
                        (72. * Qmpil(gasmix, 2, 5, i, k)) + 
                        (4. * Qmpil(gasmix, 2, 6, i, k)) - 
                        (444. * Qmpil(gasmix, 3, 3, i, k)) + 
                        (144. * Qmpil(gasmix, 3, 4, i, k)) - 
                        (16. * Qmpil(gasmix, 3, 5, i, k)) + 
                        ((32. / 3.) * Qmpil(gasmix, 4, 4, i, k))
 
                    ) )
                );   
            }
         break;
        }
     break;

    }
    
    return sum ;

}

double ZMAppendix::qmpijBar ( GasMixture* gasmix, int m , int p , int i , int j ) {
    
    int N = gasmix->getN() ;
    double ni = gasmix->getCompositionObj()->compositions(1.e-18)[i] ;
    double mi = gasmix->masses(1.e+3)[i] ;
    double T = gasmix->getTemperature() ; 
    double mj = gasmix->masses(1.e+3)[j] ;
    double nj = gasmix->getCompositionObj()->compositions(1.e-18)[j] ;

    double qij = qmpij(gasmix,m,p,i,j) ;
    double qii = qmpij(gasmix,m,p,i,i) ;

    return qij - (((nj*sqrt(mj*T))/(ni*sqrt(mi*T)))*qii*delta(m,0)*delta(p,0)) ;

}
double ZMAppendix::qmpi1Bar ( GasMixture* gasmix, int m , int p , int i ) {
    
    int N = gasmix->getN() ;
    double ni = gasmix->getCompositionObj()->compositions(1.e-18)[i] ;
    double mi = gasmix->masses(1.e+3)[i] ;
    double n1 = gasmix->getCompositionObj()->compositions(1.e-18).back() ;
    double m1 = gasmix->masses(1.e+3).back() ;
    double T = gasmix->getTemperature() ; 
    double Te = T*gasmix->theta->get() ;

    double qi1 = qmpi1(gasmix,m,p,i) ;
    double qii = qmpij(gasmix,m,p,i,i) ;

    return qi1 - (((n1*sqrt(m1*Te))/(ni*sqrt(mi*T)))*qii*delta(m,0)*delta(p,0)) ;

}
