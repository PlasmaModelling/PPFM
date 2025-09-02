 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include"Transport.h"
#include <chrono>
#include "CiBox.h"
#include "CollisionIntegral.h"
#include "GasMixture.h"

Transport::Transport ( GasMixture* mix ) { Ci = new CiBox(mix) ; }

void Transport::QtCalc ( GasMixture* gasmix) {

    double Th = gasmix->getTemperature() ; 
    double Te = gasmix->theta->get() * Th ;
    double lambda = gasmix->Comp->getDebyeLength(Te) ;
    
    int Ninteractions = Ci->InteractionsNumber() ;

    Qt.resize ( 16, std::vector<double>(Ninteractions) ) ;

    Ci->computeCollisionIntegrals( Te, Th, lambda ) ;    

    Qt.resize ( 16, std::vector<double>(Ninteractions) ) ;

    for ( int i = 0 ; i < Ninteractions ; i++ ) {
        
        std::vector<double> Q = (*Ci)[i]->omega4th ;
        
        for ( int j = 0 ; j < Q.size() ; j++ ) 
            /* CollisionIntegrals in micron^2  */
            Qt[j][i] = Q[j] * 1e-08 ;
    }   
}

double Appendix::Qmpil ( GasMixture* gasmix , int m , int p , int i , int l ) {

    int N_specs = gasmix->getN() ; 

    int il = -1;
    int mp = 0;
    double value = 0.;

    m -= 1;
    p -= 1;

    if (m == 0)
        mp = p;
    if (m == 1)
        mp = 7 + (p - 1);
    if (m == 2)
        mp = 12 + (p - 2);
    if (m == 3)
        mp = 15;

    if (l < i)
    {
        int tmp;
        tmp = l;
        l = i;
        i = tmp;
    }
    for (int I = 0; I < N_specs; I++) {
        for (int L = I; L < N_specs; L++) {
            il += 1;
            if (i == I && l == L) {
                value = Qt[mp][il];
                break;
            } 
            if (value!=0)
                break;
        }
    }
    return value;
};

double Appendix::qmpij ( GasMixture* gasmix , int m , int p , int i , int j ) {

    std::vector<double> n = gasmix->Comp->compositions(1.e-18) ;
    std::vector<double> mass = gasmix->masses(1.e+03) ;
    int N_SPC = gasmix->getN() ; 

    int l;
    double sum = 0.;

    double mi = mass[i];
    double mj = mass[j];
    double mj2 = mj * mj;
    double mj4 = mj2 * mj2;

	switch(m) {
	case 0:
		switch(p) {
		case 0:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				sum += 8.*(n[l]*sqrt(mi)/sqrt(mi+ml))*Qmpil ( gasmix , 1, 1, i, l )  * (n[i]*sqrt(ml/mj)*(delta(i,j)-delta(j,l)) - 
						n[j]*sqrt(ml*mj)*(1-delta(i,l))/mi);
			}
			break;
		case 1:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				sum += 8.*(n[l]*pow(ml,1.5)/pow(mi+ml,1.5)) *
						  (2.5*Qmpil ( gasmix , 1, 1, i, l ) -3.*Qmpil ( gasmix , 1, 2, i, l ) ) *
				  		  (delta(i,j)-delta(j,l));
			}
			sum *= n[i]*pow(mi/mj,1.5);
			break;
		case 2:
			for(l=0;l<N_SPC;l++){		
				double ml = mass[l];
				sum += 8.*(n[l]*pow(ml,2.5)/pow(mi+ml,2.5))*(delta(i,j)-delta(j,l)) *
				  		  ((35./8.)*Qmpil ( gasmix , 1, 1, i, l )  - (21./2.)*Qmpil ( gasmix , 1, 2, i, l )  + 6.*Qmpil ( gasmix , 1, 3, i, l ) );
			}
			sum *= n[i]*pow(mi/mj,2.5);
			break;
		case 3:
			for(l=0;l<N_SPC;l++){		
				double ml = mass[l];
				sum += 8.*(n[l]*pow(ml,3.5)/pow(mi+ml,3.5))*(delta(i,j)-delta(j,l)) *
						  ((105./16.)*Qmpil ( gasmix , 1, 1, i, l )  - (189./8.)*Qmpil ( gasmix , 1, 2, i, l )  + 27.*Qmpil ( gasmix , 1, 3, i, l )  - 10.*Qmpil ( gasmix , 1, 4, i, l ) );
			}
			sum *= n[i]*pow(mi/mj,3.5);
			break;
		}
		break;
		
	case 1:
		switch(p) {
		case 0:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				sum += 8.*(n[l]*pow(ml,1.5)/pow(mi+ml,1.5)) *
					  (2.5*Qmpil ( gasmix , 1, 1, i, l ) -3.*Qmpil ( gasmix , 1, 2, i, l ) ) *
					  (delta(i,j)-delta(j,l));
			}
			sum *= n[i]*sqrt(mi/mj);
			break;
		case 1:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				double ml2 = ml*ml;
		
				sum += 8.*(n[l]*sqrt(ml)/pow(mi+ml,2.5)) *
			  			  ((delta(i,j)-delta(j,l))*(1.25*(6.*mj2+5.*ml2)*Qmpil ( gasmix , 1, 1, i, l )  -
				  					15.*ml2*Qmpil ( gasmix , 1, 2, i, l ) +12.*ml2*Qmpil ( gasmix , 1, 3, i, l ) )+
				  					(delta(i,j)+delta(j,l))*4.*mj*ml*Qmpil ( gasmix , 2, 2, i, l ) );
			}
			sum *= n[i]*pow(mi/mj,1.5);
			break;
		case 2:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				double ml2 = ml*ml;
				
				sum += 8.*(n[l]*pow(ml,1.5)/pow(mi+ml,3.5))*
						  ((delta(i,j)-delta(j,l))*((35./16.)*(12.*mj2+5.*ml2)*Qmpil ( gasmix , 1, 1, i, l )  -
						  							(63./2.)*(mj2+(5./4.)*ml2)*Qmpil ( gasmix , 1, 2, i, l )  +
						  							 57.*ml2*Qmpil ( gasmix , 1, 3, i, l )  - 30.*ml2*Qmpil ( gasmix , 1, 4, i, l ) ) +
						   (delta(i,j)+delta(j,l))*(14.*mj*ml*Qmpil ( gasmix , 2, 2, i, l )  - 16.*mj*ml*Qmpil ( gasmix , 2, 3, i, l ) ));
			}
			sum *= n[i]*pow(mi/mj,2.5);
			break;
		case 3:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				double ml2 = ml*ml;
		
				sum += 8.*(n[l]*pow(ml,2.5)/pow(mi+ml,4.5))*
						  ((delta(i,j)-delta(j,l))*((105./32.)*(18.*mj2+5.*ml2)*Qmpil ( gasmix , 1, 1, i, l )  -
						  							(63./4.)*(9.*mj2+5.*ml2)*Qmpil ( gasmix , 1, 2, i, l )  +
						  							 81.*(mj2+2.*ml2)*Qmpil ( gasmix , 1, 3, i, l )  -
						  							 160.*ml2*Qmpil ( gasmix , 1, 4, i, l )  + 60.*ml2*Qmpil ( gasmix , 1, 5, i, l ) ) +
						   (delta(i,j)+delta(j,l))*mj*ml*((63./2.)*Qmpil ( gasmix , 2, 2, i, l )  - 72.*Qmpil ( gasmix , 2, 3, i, l )  + 40.*Qmpil ( gasmix , 2, 4, i, l ) ));
			}
			sum *= n[i]*pow(mi/mj,3.5);
			break;
		}
		break;
		
	case 2:
		switch(p) {
		case 0:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
		
				sum += 8.*(n[l]*pow(ml,2.5)/pow(mi+ml,2.5))*(delta(i,j)-delta(j,l)) *
						  ((35./8.)*Qmpil ( gasmix , 1, 1, i, l )  - (21./2.)*Qmpil ( gasmix , 1, 2, i, l )  + 6.*Qmpil ( gasmix , 1, 3, i, l ) );
			}
			sum *= n[i]*sqrt(mi/mj);
			break;
		case 1:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				double ml2 = ml*ml;
		
				sum += 8.*(n[l]*pow(ml,1.5)/pow(mi+ml,3.5))*
						  ((delta(i,j)-delta(j,l))*((35./16.)*(12.*mj2+5.*ml2)*Qmpil ( gasmix , 1, 1, i, l )  -
						  							(63./2.)*(mj2+(5./4.)*ml2)*Qmpil ( gasmix , 1, 2, i, l )  +
						  							 57.*ml2*Qmpil ( gasmix , 1, 3, i, l )  - 30.*ml2*Qmpil ( gasmix , 1, 4, i, l ) ) +
						   (delta(i,j)+delta(j,l))*(14.*mj*ml*Qmpil ( gasmix , 2, 2, i, l )  - 16.*mj*ml*Qmpil ( gasmix , 2, 3, i, l ) ));
			}
			sum *= n[i]*pow(mi/mj,1.5);
			break;
		case 2:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				double ml2 = ml*ml;
				double ml3 = ml2*ml;
				double ml4 = ml2*ml2;
				
				sum += 8.*(n[l]*sqrt(ml)/pow(mi+ml,4.5))*
						  ((delta(i,j)-delta(j,l))*((35./64.)*(40.*mj4+168.*mj2*ml2+35.*ml4)*Qmpil ( gasmix , 1, 1, i, l )  -
						  							(21./8.)*ml2*(84.*mj2+35.*ml2)*Qmpil ( gasmix , 1, 2, i, l )  +
						  							 1.5*ml2*(108.*mj2+133.*ml2)*Qmpil ( gasmix , 1, 3, i, l )  -
						  							 210.*ml4*Qmpil ( gasmix , 1, 4, i, l )  + 90.*ml4*Qmpil ( gasmix , 1, 5, i, l )  +
						  							 24.*mj2*ml2*Qmpil ( gasmix , 3, 3, i, l ) ) +
						   (delta(i,j)+delta(j,l))*(7.*mj*ml*(4.*mj2+7.*ml2)*Qmpil ( gasmix , 2, 2, i, l )  -
						   						    112.*mj*ml3*Qmpil ( gasmix , 2, 3, i, l )  + 80*mj*ml3*Qmpil ( gasmix , 2, 4, i, l ) ));
			}
			sum *= n[i]*pow(mi/mj,2.5);
			break;
		case 3:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				double ml2 = ml*ml;
				double ml3 = ml2*ml;
				double ml4 = ml2*ml2;
				
				sum += 8.*(n[l]*pow(ml,1.5)/pow(mi+ml,5.5))*
						  ((delta(i,j)-delta(j,l))*((105./128.)*(120.*mj4+252.*mj2*ml2+35.*ml4)*Qmpil ( gasmix , 1, 1, i, l )  -
						  							(63./64.)*(120.*mj4+756.*mj2*ml2+175.*ml4)*Qmpil ( gasmix , 1, 2, i, l )  +
						  							(9./4.)*ml2*(450.*mj2+217.*ml2)*Qmpil ( gasmix , 1, 3, i, l )  -
						  							(5./2.)*ml2*(198.*mj2+301.*ml2)*Qmpil ( gasmix , 1, 4, i, l )  +
						  							 615.*ml4*Qmpil ( gasmix , 1, 5, i, l )  - 210.*ml4*Qmpil ( gasmix , 1, 6, i, l )  +
						  							 108.*mj2*ml2*Qmpil ( gasmix , 3, 3, i, l )  - 120.*mj2*ml2*Qmpil ( gasmix , 3, 4, i, l ) ) +
						   (delta(i,j)+delta(j,l))*((63./4.)*mj*ml*(8.*mj2+7.*ml2)*Qmpil ( gasmix , 2, 2, i, l )  -
						   						    18.*mj*ml*(8.*mj2+21.*ml2)*Qmpil ( gasmix , 2, 3, i, l )  +
						   						    500.*mj*ml3*Qmpil ( gasmix , 2, 4, i, l )  - 240.*mj*ml3*Qmpil ( gasmix , 2, 5, i, l ) ));
			}
			sum *= n[i]*pow(mi/mj,3.5);
			break;
		}
		break;

	case 3:
		switch(p) {
		case 0:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				
				sum += 8.*(n[l]*pow(ml,3.5)/pow(mi+ml,3.5))*(delta(i,j)-delta(j,l)) *
						  ((105./16.)*Qmpil ( gasmix , 1, 1, i, l )  - (189./8.)*Qmpil ( gasmix , 1, 2, i, l )  + 27.*Qmpil ( gasmix , 1, 3, i, l )  - 10.*Qmpil ( gasmix , 1, 4, i, l ) );
			}
			sum *= n[i]*sqrt(mi/mj);
			break;
		case 1:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				double ml2 = ml*ml;
				
				sum += 8.*(n[l]*pow(ml,2.5)/pow(mi+ml,4.5))*
						  ((delta(i,j)-delta(j,l))*((105./32.)*(18.*mj2+5.*ml2)*Qmpil ( gasmix , 1, 1, i, l )  -
						  							(63./4.)*(9.*mj2+5.*ml2)*Qmpil ( gasmix , 1, 2, i, l )  +
						  							 81.*(mj2+2.*ml2)*Qmpil ( gasmix , 1, 3, i, l )  -
						  							 160.*ml2*Qmpil ( gasmix , 1, 4, i, l )  + 60.*ml2*Qmpil ( gasmix , 1, 5, i, l ) ) +
						   (delta(i,j)+delta(j,l))*mj*ml*((63./2.)*Qmpil ( gasmix , 2, 2, i, l )  - 72.*Qmpil ( gasmix , 2, 3, i, l )  + 40.*Qmpil ( gasmix , 2, 4, i, l ) ));
			}
			sum *= n[i]*pow(mi/mj,1.5);
			break;
		case 2:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				double ml2 = ml*ml;
				double ml3 = ml2*ml;
				double ml4 = ml2*ml2;
				
				sum += 8.*(n[l]*pow(ml,1.5)/pow(mi+ml,5.5))*
						  ((delta(i,j)-delta(j,l))*((105./128.)*(120.*mj4+252.*mj2*ml2+35.*ml4)*Qmpil ( gasmix , 1, 1, i, l )  -
						  							(63./64.)*(120.*mj4+756.*mj2*ml2+175.*ml4)*Qmpil ( gasmix , 1, 2, i, l )  +
						  							(9./4.)*ml2*(450.*mj2+217.*ml2)*Qmpil ( gasmix , 1, 3, i, l )  -
						  							(5./2.)*ml2*(198.*mj2+301.*ml2)*Qmpil ( gasmix , 1, 4, i, l )  +
						  							 615.*ml4*Qmpil ( gasmix , 1, 5, i, l )  - 210.*ml4*Qmpil ( gasmix , 1, 6, i, l )  +
						  							 108.*mj2*ml2*Qmpil ( gasmix , 3, 3, i, l )  - 120.*mj2*ml2*Qmpil ( gasmix , 3, 4, i, l ) ) +
						   (delta(i,j)+delta(j,l))*((63./4.)*mj*ml*(8.*mj2+7.*ml2)*Qmpil ( gasmix , 2, 2, i, l )  -
						   						    18.*mj*ml*(8.*mj2+21.*ml2)*Qmpil ( gasmix , 2, 3, i, l )  +
						   						    500.*mj*ml3*Qmpil ( gasmix , 2, 4, i, l )  - 240.*mj*ml3*Qmpil ( gasmix , 2, 5, i, l ) ));
			}
			sum *= n[i]*pow(mi/mj,2.5);
			break;
		case 3:
			for(l=0;l<N_SPC;l++){
				double ml = mass[l];
				double ml2 = ml*ml;
				double ml3 = ml2*ml;
				double ml4 = ml2*ml2;
		
				sum += 8.*(n[l]*sqrt(ml)/pow(mi+ml,13./2.))*
						  ((delta(i,j)-delta(j,l))*((105./256.)*(112.*pow(mj,6)+1080.*ml2*mj4+
						  										 1134.*mj2*ml4+105.*pow(ml,6))*Qmpil ( gasmix , 1, 1, i, l )  -
						  							(567./64.)*ml2*(120.*mj4+252.*mj2*ml2+35.*ml4)*Qmpil ( gasmix , 1, 2, i, l )  +
						  							(27./16.)*ml2*(440.*mj4+2700.*mj2*ml2+651.*ml4)*Qmpil ( gasmix , 1, 3, i, l )  -
						  							(15./2.)*ml4*(594.*mj2+301.*ml2)*Qmpil ( gasmix , 1, 4, i, l )  +
						  							(135./2.)*ml4*(26.*mj2+41.*ml2)*Qmpil ( gasmix , 1, 5, i, l )  - 1890.*pow(ml,6)*Qmpil ( gasmix , 1, 6, i, l )  +
						  							 560.*pow(ml,6)*Qmpil ( gasmix , 1, 7, i, l )  + 18.*mj2*ml2*(10.*mj2+27.*ml2)*Qmpil ( gasmix , 3, 3, i, l )  -
						  							 1080.*mj2*ml4*Qmpil ( gasmix , 3, 4, i, l )  + 720.*mj2*ml4*Qmpil ( gasmix , 3, 5, i, l ) ) +
						   (delta(i,j)+delta(j,l))*((189./16.)*mj*ml*(8.*mj4+48.*pow(ml*mj,2)+21.*ml4)*Qmpil ( gasmix , 2, 2, i, l )  -
						   						    162.*mj*ml3*(8.*mj2+7.*ml2)*Qmpil ( gasmix , 2, 3, i, l )  +
						   						    10.*mj*ml3*(88.*mj2+225.*ml2)*Qmpil ( gasmix , 2, 4, i, l )  -
						   						    2160.*mj*pow(ml,5)*Qmpil ( gasmix , 2, 5, i, l )  + 840.*mj*pow(ml,5)*Qmpil ( gasmix , 2, 6, i, l )  +
						   						    64.*pow(mj*ml,3)*Qmpil ( gasmix , 4, 4, i, l ) ));
			}
			sum *= n[i]*pow(mi/mj,3.5);
			break;
		}
		break;		
	}
	
	return sum;
}

double Appendix::qsimpmpij( GasMixture* gasmix , int m , int p ) {

    std::vector<double> n = gasmix->Comp->compositions(1.e-18) ;
    std::vector<double> mass = gasmix->masses(1.e+03) ;
    int N_SPC = gasmix->getN() ;

    int l;
    double sum = 0.;

    switch (m)
    {
    case 0:
        switch (p)
        {
        case 0:
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * Qmpil ( gasmix , 1, 1, N_SPC - 1, l ) ;
            }
            sum *= 8 * n[N_SPC - 1];
            break;
        case 1:
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((5. / 2.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - 3. * Qmpil ( gasmix , 1, 2, N_SPC - 1, l ) );
            }
            sum *= 8 * n[N_SPC - 1];
            break;
        case 2:
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((35. / 8.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (21. / 2.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + 6. * Qmpil ( gasmix , 1, 3, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            break;

        case 3:
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((105. / 16.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (189. / 8.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + 27. * Qmpil ( gasmix , 1, 3, N_SPC - 1, l )  - 10. * Qmpil ( gasmix , 1, 4, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            break;
        }
        break;

    case 1:
        switch (p)
        {
        case 0:
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((5. / 2.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - 3. * Qmpil ( gasmix , 1, 2, N_SPC - 1, l ) );
            }
            sum *= 8 * n[N_SPC - 1];
            break;
        case 1:
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((25. / 4.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (15.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + 12. * Qmpil ( gasmix , 1, 3, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            sum += 8. * sqrt(2.) * pow(n[N_SPC - 1], 2) * Qmpil ( gasmix , 2, 2, N_SPC - 1, N_SPC - 1 ) ;
            break;
        case 2:
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((175. / 16.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (315. / 8.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + 57. * Qmpil ( gasmix , 1, 3, N_SPC - 1, l )  - 30. * Qmpil ( gasmix , 1, 4, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            sum += 8. * sqrt(2.) * pow(n[N_SPC - 1], 2) * ((7. / 4.) * Qmpil ( gasmix , 2, 2, N_SPC - 1, N_SPC - 1 )  - 2. * Qmpil ( gasmix , 2, 3, N_SPC - 1, N_SPC - 1 ) );
            break;
        case 3:
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((525. / 32.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (315. / 4.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + 162. * Qmpil ( gasmix , 1, 3, N_SPC - 1, l )  - 160. * Qmpil ( gasmix , 1, 4, N_SPC - 1, l )  + 60. * Qmpil ( gasmix , 1, 5, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            sum = 8. * sqrt(2.) * pow(n[N_SPC - 1], 2) * ((63. / 32.) * Qmpil ( gasmix , 2, 2, N_SPC - 1, N_SPC - 1 )  - (9. / 2.) * Qmpil ( gasmix , 2, 3, N_SPC - 1, N_SPC - 1 )  + (5. / 2.) * Qmpil ( gasmix , 2, 4, N_SPC - 1, N_SPC - 1 ) ) + sum;
            break;
        }
        break;

    case 2:
        switch (p)
        {
        case 0:
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((35. / 8.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (21. / 2.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + 6. * Qmpil ( gasmix , 1, 3, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            break;
        case 1:
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((175. / 16.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (315. / 8.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + 57. * Qmpil ( gasmix , 1, 3, N_SPC - 1, l )  - 30. * Qmpil ( gasmix , 1, 4, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            sum += 8. * sqrt(2.) * pow(n[N_SPC - 1], 2) * ((7. / 4.) * Qmpil ( gasmix , 2, 2, N_SPC - 1, N_SPC - 1 )  - 2. * Qmpil ( gasmix , 2, 3, N_SPC - 1, N_SPC - 1 ) );
            break;
        case 2:
            
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((1225. / 64.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (735. / 8.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + (399. / 2.) * Qmpil ( gasmix , 1, 3, N_SPC - 1, l )  - 210. * Qmpil ( gasmix , 1, 4, N_SPC - 1, l )  + 90. * Qmpil ( gasmix , 1, 5, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            sum += 8. * sqrt(2.) * pow(n[N_SPC - 1], 2) * ((77. / 16.) * Qmpil ( gasmix , 2, 2, N_SPC - 1, N_SPC - 1 )  - 7. * Qmpil ( gasmix , 2, 3, N_SPC - 1, N_SPC - 1 )  + 5. * Qmpil ( gasmix , 2, 4, N_SPC - 1, N_SPC - 1 ) );
            break;
        case 3:
            
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((3675. / 128.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (11025. / 64.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + (1953. / 4.) * Qmpil ( gasmix , 1, 3, N_SPC - 1, l )  - (1505 / 2.) * Qmpil ( gasmix , 1, 4, N_SPC - 1, l )  + 615. * Qmpil ( gasmix , 1, 5, N_SPC - 1, l )  - 210. * Qmpil ( gasmix , 1, 6, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            sum += 8. * sqrt(2.) * pow(n[N_SPC - 1], 2) * ((945. / 128.) * Qmpil ( gasmix , 2, 2, N_SPC - 1, N_SPC - 1 )  - (261. / 16.) * Qmpil ( gasmix , 2, 3, N_SPC - 1, N_SPC - 1 )  + (125. / 8.) * Qmpil ( gasmix , 2, 4, N_SPC - 1, N_SPC - 1 )  - (15. / 2.) * Qmpil ( gasmix , 2, 5, N_SPC - 1, N_SPC - 1 ) );
            break;
        }
        break;

    case 3:
        switch (p)
        {
        case 0:
            
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((105. / 16.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (189. / 8.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + 27. * Qmpil ( gasmix , 1, 3, N_SPC - 1, l )  - 10. * Qmpil ( gasmix , 1, 4, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            break;
        case 1:
            
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((525. / 32.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (315. / 4.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + 162. * Qmpil ( gasmix , 1, 3, N_SPC - 1, l )  - 160. * Qmpil ( gasmix , 1, 4, N_SPC - 1, l )  + 60. * Qmpil ( gasmix , 1, 5, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            sum = 8. * sqrt(2.) * pow(n[N_SPC - 1], 2) * ((63. / 32.) * Qmpil ( gasmix , 2, 2, N_SPC - 1, N_SPC - 1 )  - (9. / 2.) * Qmpil ( gasmix , 2, 3, N_SPC - 1, N_SPC - 1 )  + (5. / 2.) * Qmpil ( gasmix , 2, 4, N_SPC - 1, N_SPC - 1 ) ) + sum;
            break;
        case 2:
            
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((3675. / 128.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (11025. / 64.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + (1953. / 4.) * Qmpil ( gasmix , 1, 3, N_SPC - 1, l )  - (1505 / 2.) * Qmpil ( gasmix , 1, 4, N_SPC - 1, l )  + 615. * Qmpil ( gasmix , 1, 5, N_SPC - 1, l )  - 210. * Qmpil ( gasmix , 1, 6, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            sum += 8. * sqrt(2.) * pow(n[N_SPC - 1], 2) * ((945. / 128.) * Qmpil ( gasmix , 2, 2, N_SPC - 1, N_SPC - 1 )  - (261. / 16.) * Qmpil ( gasmix , 2, 3, N_SPC - 1, N_SPC - 1 )  + (125. / 8.) * Qmpil ( gasmix , 2, 4, N_SPC - 1, N_SPC - 1 )  - (15. / 2.) * Qmpil ( gasmix , 2, 5, N_SPC - 1, N_SPC - 1 ) );
            break;
        case 3:
            
            for (l = 0; l < N_SPC - 1; l++)
            {
                sum += n[l] * ((11025. / 256.) * Qmpil ( gasmix , 1, 1, N_SPC - 1, l )  - (19845. / 64.) * Qmpil ( gasmix , 1, 2, N_SPC - 1, l )  + (17577. / 16.) * Qmpil ( gasmix , 1, 3, N_SPC - 1, l )  - (4515 / 2.) * Qmpil ( gasmix , 1, 4, N_SPC - 1, l )  + (5535. / 2.) * Qmpil ( gasmix , 1, 5, N_SPC - 1, l )  - 1890. * Qmpil ( gasmix , 1, 6, N_SPC - 1, l )  + 560. * Qmpil ( gasmix , 1, 7, N_SPC - 1, l ) );
            }
            sum *= 8. * n[N_SPC - 1];
            sum += 8. * sqrt(2.) * pow(n[N_SPC - 1], 2) * ((14553. / 1024.) * Qmpil ( gasmix , 2, 2, N_SPC - 1, N_SPC - 1 )  - (1215. / 32.) * Qmpil ( gasmix , 2, 3, N_SPC - 1, N_SPC - 1 )  + (1565. / 32.) * Qmpil ( gasmix , 2, 4, (N_SPC - 1), N_SPC - 1 )  - (135. / 4.) * Qmpil ( gasmix , 2, 5, N_SPC - 1, N_SPC - 1 )  + (105. / 8) * Qmpil ( gasmix , 2, 6, N_SPC - 1, N_SPC - 1 )  + Qmpil ( gasmix , 4, 4, N_SPC - 1, N_SPC - 1 ) );
            break;
        }
        break;
    }
    return sum;
};

double Appendix::qcapmpij( GasMixture* gasmix, int m, int p, int i, int j) {

std::vector<double> n = gasmix->Comp->compositions(1.e-18) ; 
std::vector<double> mass = gasmix->masses(1.e+03) ;
int N_SPC = gasmix->getN() ;

int l;
double sum = 0.;

double mi = mass[i];
double mj = mass[j];
double mj2 = mj * mj;

switch (m)
{
case 0:
    switch (p)
    {
    case 0:

        for (l = 0; l < N_SPC; l++)
        {
            double ml = mass[l];
            sum += 8. * (n[l] * sqrt(ml) / pow(mi + ml, 1.5)) * ((10. / 3) * Qmpil ( gasmix , 1, 1, i, l )  * (delta(i, j) - delta(j, l)) * mj + 2. * ml * Qmpil ( gasmix , 2, 2, i, l )  * (delta(i, j) + delta(j, l)));
        }
        sum *= n[i] * (mi / mj);
        break;
    case 1:

        for (l = 0; l < N_SPC; l++)
        {
            double ml = mass[l];
            sum += 8. * (n[l] * pow(ml, 1.5) / pow(mi + ml, 2.5)) * ((delta(i, j) - delta(j, l)) * mj * ((35. / 3) * Qmpil ( gasmix , 1, 1, i, l )  - 14. * Qmpil ( gasmix , 1, 2, i, l ) ) + (delta(i, j) + delta(j, l)) * ml * (7. * Qmpil ( gasmix , 2, 2, i, l )  - 8. * Qmpil ( gasmix , 2, 3, i, l ) ));
        }
        sum *= n[i] * pow(mi / mj, 2);
        break;
    }
    break;
case 1:
    switch (p)
    {
    case 0:

        for (l = 0; l < N_SPC; l++)
        {
            double ml = mass[l];
            sum += 8. * (n[l] * pow(ml, 1.5) / pow(mi + ml, 2.5)) * ((delta(i, j) - delta(j, l)) * mj * ((35. / 3) * Qmpil ( gasmix , 1, 1, i, l )  - 14. * Qmpil ( gasmix , 1, 2, i, l ) ) + (delta(i, j) + delta(j, l )) * ml * (7. * Qmpil ( gasmix , 2, 2, i, l )  - 8. * Qmpil ( gasmix , 2, 3, i, l ) ));
        }
        sum *= n[i] * (mi / mj);
        break;
    case 1:

        for (l = 0; l < N_SPC; l++)
        {
            double ml = mass[l];
            double ml2 = ml * ml;

            sum += 8. * (n[l] * sqrt(ml) / pow(mi + ml, 3.5)) *
                    ((delta(i, j) - delta(j, l)) * mj * ((1. / 6.) * (140. * mj2 + 245. * ml2) * Qmpil ( gasmix , 1, 1, i, l )  - ml2 * (98. * Qmpil ( gasmix , 1, 2, i, l )  - 64. * Qmpil ( gasmix , 1, 3, i, l )  - 24. * Qmpil ( gasmix , 3, 3, i, l ) )) +
                    (delta(i, j) + delta(j, l)) * ml * ((1. / 6.) * (154 * mj2 + 147. * ml2) * Qmpil ( gasmix , 2, 2, i, l )  - ml2 * (56. * Qmpil ( gasmix , 2, 3, i, l )  - 40. * Qmpil ( gasmix , 2, 4, i, l ) )));
        }
        sum *= n[i] * pow(mi / mj, 2);
        break;
    }
    break;
}
return sum;
};

