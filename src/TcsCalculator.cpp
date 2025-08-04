
#include "Interaction.h"
#include "Potential.h"
#include "TcsCalculator.h"
#include <stdexcept>
#include <numbers>
#include <filesystem>
#include <regex>
#include <set>

#ifdef PPFM_USE_CURL
    #include <curl/curl.h>
    void DcsLoader::DownloadAndSaveData(const std::string& hyperref) {

        CURL* curl = curl_easy_init();

        if (!curl) 
            throw std::runtime_error("Failed to initialize CURL.");
        
        std::stringstream raw_data;
        curl_easy_setopt(curl, CURLOPT_URL, hyperref.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &raw_data);

        CURLcode res = curl_easy_perform(curl);
        if (res != CURLE_OK) {
            
            curl_easy_cleanup(curl);
            throw std::runtime_error("CURL error: " + std::string(curl_easy_strerror(res)));
        
        }
        
        curl_easy_cleanup(curl);

        if (raw_data.str().empty()) 
            throw std::runtime_error("Downloaded content is empty.");
        
        std::filesystem::path datadir = BuildFilePath("Differential_Cross_Sections");
        std::filesystem::create_directories(datadir);

        std::ofstream output_file(datadir / BuildFileName(Name));
        if (!output_file.is_open()) 
            throw std::runtime_error("Unable to open output file.");
        
        output_file << raw_data.str();
        output_file.close();
    }
#else
    void DcsLoader::DownloadAndSaveData(const std::string& hyperref) {
    
            std::cerr << "[PPFM] libcurl not available: Download not available skipped.\n";

    }
#endif

CsCalculator::CsCalculator( InteractionInterface* i ) { Name = i->InteractionName(); } 

std::vector<double> CsCalculator::operator()(int l) {
    
    l--;
    return (l>-1 && l<4)? column(Q,l):throw std::invalid_argument("Invalid l");

}

// ______________________ Implementation CsHolder _____________________

/// @brief Construct a new object
CsHolder::CsHolder( CsCalculator* Qe , CsCalculator* Qin ) : Qe(Qe), Qin(Qin) {Name = Qe->Name;}

/// @brief Construct a new object
CsHolder::CsHolder( InteractionInterface* i ) : CsCalculator(i) {}

void CsHolder::Compute()  { Qe->Compute() ; Qin->Compute() ; computed = true; }

// ______________________ Implementation MultiCs _____________________

// Constructor
MultiCs::MultiCs( InteractionInterface* i, std::vector<CsCalculator*> c, std::vector<double> gs ) : CsCalculator(i) , Qs(c), statesG(gs) {}  

MultiCs::MultiCs( InteractionInterface* i, std::vector<CsCalculator*> c ) : CsCalculator(i) , Qs(c) {}  

CsCalculator*& MultiCs::operator[](int i) {return Qs[i];}

// std::vector<double> MultiCs::operator() ( int s, int l ) { return (*Qs[s])(l); };

int MultiCs::Size() { 

    if ( Qs.size() == statesG.size() ) 
        return Qs.size() ;  
    else 
        throw std::invalid_argument("Invalid Multiple Cross Sections,"
        "\n different sizes between state potentials and degeneracies ") ;

}

void MultiCs::Compute() { for ( int i = 0; i < Qs.size(); i++ ) { Qs[i]->Compute(); } computed = true; }

// ______________________ Implementation ThresholdCs _____________________

// Constructor
ThresholdCs::ThresholdCs ( InteractionInterface* i, std::vector<CsCalculator*> c, std::vector<double> Elim ) : 
    MultiCs(i,c), Elim(Elim) {
        
    if( Qs.size() != (Elim.size()+1) ){
        throw std::runtime_error("Bad threshold cross section initialization. \n"
        "CrossSection calculator must be one larger than the desired thresholds");
    }

}

void ThresholdCs::Compute() {
    
    if (Qs.empty() || Elim.size() != Qs.size() - 1) 
        throw std::runtime_error("Invalid state: empty calculators or incorrect sizes\n"
                                 "calculators must be equals to limit energies plus one.");
    
    double previousLimit = 0.0;
    
    for (int i = 0; i < Qs.size(); ++i) {
        
        if (Qs[i]->E.empty()) 
            throw std::runtime_error("Threshold calculator: " + std::to_string(i) + " has empty energy range.");
        
        int it = 0;
        while (it < Qs[i]->E.size() && Qs[i]->E[it] < previousLimit) 
            ++it;
        
        int itEnd = it;
        while (itEnd < Qs[i]->E.size() && (i == Qs.size() - 1 || Qs[i]->E[itEnd] < Elim[i]))
            ++itEnd;
        
        Qs[i]->E.assign(Qs[i]->E.begin() + it, Qs[i]->E.begin() + itEnd);
        
        if (i < Elim.size()) 
            previousLimit = Elim[i];
    }

    for (int i = 0; i < Qs.size(); ++i) {
        if (Qs[i]->Q.empty()) 
            Qs[i]->Compute();  
    }

    E.clear();
    E.shrink_to_fit();
    for (const auto& q : Qs) 
        E = concatenate(E, q->E);

    Q.clear();
    Q.shrink_to_fit();
    for (const auto& q : Qs) 
        Q = concatenate(Q, q->Q);

    computed = true ;
}

// __________________________ Implementation parser ___________________________

void TcsParser::TcsParse( std::ifstream& file, 
    std::vector<double>& E , std::vector<std::vector<double>>& Q ) {

    E.clear();
    Q.clear();
    E.shrink_to_fit();
    Q.shrink_to_fit();

    std::string line ; 

    std::getline(file, line);

    while (std::getline(file,line)) {
    
        std::istringstream lineStream(line);
        double energy, qmE;
        char separator;

        if (lineStream >> energy >> separator >> qmE && separator == ',') {

            E.push_back(energy);

            // Creazione della riga con valori Q^1 = Q^2 = Q^3 = Q^4 = Q^m(E)
            /// DA IMPLEMENTARE : in base quante colonne del file tot l=1,2,3,4... 
            Q.push_back(std::vector<double>(4, qmE));

        }
    }

    // Verifica della validità dei dati
    if (E.empty() || Q.size() != E.size()) 
        throw std::runtime_error("Error parsing file: data mismatch.");
}

// ______________________ Implementation ElasticLoader  _____________________

ElasticLoader::ElasticLoader( InteractionInterface* i ) : CsCalculator(i) { Init(); } ;
    
/// @brief Construct a new object
ElasticLoader::ElasticLoader( const std::string& prefix, InteractionInterface* interaction ) : 
    CsCalculator(interaction) { this->customPrefix = prefix + "_"; Init(); }

void ElasticLoader::Compute()  { Init(); computed = true ; }

void ElasticLoader::Init()  { if(!loaded) LoadData("Momentum_Transport_Cross_Section",Name); }

std::string ElasticLoader::BuildFileName(const std::string& name)  {

    return "TCS_"+customPrefix+name+".csv";

}

void ElasticLoader::ParseFile(std::ifstream& file)  { TcsParse(file,E,Q); }

// ______________________ Implementation InelasticLoader  _____________________


InelasticLoader::InelasticLoader( InteractionInterface* i ) : CsCalculator(i) { Init(); } ;

/// @brief Construct a new object
InelasticLoader::InelasticLoader( const std::string& prefix, InteractionInterface* interaction ) : 
    CsCalculator(interaction) { this->customPrefix = prefix + "_"; Init(); }

void InelasticLoader::Compute()  { Init(); computed = true; }

void InelasticLoader::Init()  { if( !loaded ) LoadData("Momentum_Transport_Cross_Section",Name); }

std::string InelasticLoader::BuildFileName( const std::string& name ) {

    return "TCS_In_"+customPrefix+name+".csv" ; 

}

void InelasticLoader::ParseFile(std::ifstream& file)  { TcsParse(file,E,Q); }

//______________________ PhaseShiftsLoader implementations ______________________

void PhaseShiftsLoader::ComputeFromPhaseShifts() {

    Q.resize(E.size(), std::vector<double>(4,0.)) ;    

    /* Larger phase shift matrix as terated in original code for Ar_e- interaction.
    delta[l] > delta[l+1] > delta[l+2]... delta[infinity] faints.  */
    std::vector<std::vector<double>> eta = etaL ; 
    for (auto& row : eta)
        row.insert(row.end(), 5, 0.0) ;     

    Species* sp1 = inT->GetSp1() ;
    Species* sp2 = inT->GetSp2() ;

    double mi,mj,mij ; 
    mi = sp1->getMass();
    mj = sp2->getMass();

    mij = mi * mj / ( mi + mj ) ; 

    #pragma omp parallel for
    for (int i = 0; i < E.size(); i++) {

        double gij = sqrt((2. * E[i] * qe) / mij);
        double k = (gij * mij * 2. * std::numbers::pi) / hPlanck;

        double q0 = 0., q1 = 0., q2 = 0., q3 = 0.;

        #pragma omp parallel for reduction(+:q0, q1, q2, q3)
        for (int l = 0; l < 4; l++) {

            q0 += (l + 1.) * pow(sin(eta[i][l] - eta[i][l + 1]), 2.);

            q1 += ((l + 1.) * (l + 2.) / (2. * l + 3.)) * pow(sin(eta[i][l] - eta[i][l + 2]), 2.);

            q2 += ((l + 1.) / (2. * l + 5.)) *
                ((((l + 3.) * (l + 2.) / (2. * l + 3.)) * pow(sin(eta[i][l] - eta[i][l + 3]), 2.)) +
                ((2. * (l * l + 2. * l - 1.) / (2. * l - 1.)) * pow(sin(eta[i][l] - eta[i][l + 1]), 2.)));

            q3 += (((l + 2.) * (l + 1.)) / ((2. * l + 7.) * (2. * l + 3.))) *
                ((((l + 4.) * (l + 3.)) / (2. * l + 5.)) * pow(sin(eta[i][l] - eta[i][l + 4]), 2.)) +
                ((2. * (2. * l * l + 6. * l - 3. * l) / (2. * l - 1.)) * pow(sin(eta[i][l] - eta[i][l + 1]), 2.));
        }

        const double scale = (4. * std::numbers::pi) * pow(k, -2.) * 1e+20;

        Q[i][0] = q0 * scale;
        Q[i][1] = q1 * scale;
        Q[i][2] = q2 * scale;
        Q[i][3] = q3 * scale;
    }
}

void PhaseShiftsLoader::ParseFile(std::ifstream& file) {

    std::string line;
    
    // Jumps the first row
    if (!std::getline(file, line)) return;

    // parsing temporary variables
    std::vector<std::vector<double>> tempEtaL;
    std::vector<double> tempE;

    // necessary number of columns
    int maxColumns = 4; 
    int actualColumns = 0;

    while (std::getline(file, line)) {

        std::stringstream ss(line);
        std::vector<double> rowData;
        double value;

        // rowdata storing
        while (ss >> value) {
    
            rowData.push_back(value);
            if (ss.peek() == ',') ss.ignore(); 
    
        }

        // minimum dimension needed
        if (rowData.size() < 2) continue; // Serve almeno una colonna per E e una per eta

        // maximum available columns in datafile
        int numEtaColumns = rowData.size() - 1;
        if (numEtaColumns > actualColumns) actualColumns = numEtaColumns;

        // Save first column as the energy
        tempE.push_back(rowData[0]);

        // Save etaL columns
        std::vector<double> etaRow(rowData.begin() + 1, rowData.end());

        // if nColumns < 4 repeat the last available value 
        while (etaRow.size() < maxColumns) 
            etaRow.push_back(etaRow.back()); 
        

        tempEtaL.push_back(etaRow);
    
    }

    int finalColumns = (actualColumns >= maxColumns) ? actualColumns : maxColumns;
    
    E = std::move(tempE);
 
    etaL.resize(E.size(), std::vector<double>(finalColumns, 0.0));
    #pragma omp parallel for collapse(2)
    for ( int i = 0; i < tempEtaL.size(); ++i) 
        for (int j = 0; j < tempEtaL[0].size(); ++j) 
            etaL[i][j] = tempEtaL[i][j];

}

PhaseShiftsLoader::PhaseShiftsLoader( InteractionInterface* i ) : CsCalculator(i) , inT(i) { Init(); } ;

/// @brief Construct a new object
PhaseShiftsLoader::PhaseShiftsLoader( const std::string& prefix, InteractionInterface* i ) : 
    CsCalculator(i) , inT(i) { this->customPrefix = prefix + "_" ; Init(); }

void PhaseShiftsLoader::Compute()  {

    Init();

    ComputeFromPhaseShifts() ; 

    computed = true ;

}

void PhaseShiftsLoader::Init()  { if( !loaded ) LoadData("Phase_Shifts",Name); }

std::string PhaseShiftsLoader::BuildFileName( const std::string& name ) { 
    return "PS_"+customPrefix+name+".csv"; 
}

// ______________________ Implementation AdaptDeflAngle _____________________

// Integrando di formula 1 COLONNA
double AdaptDeflAngle::IntegrandKi(double r,double b,double E) {
    
    if (this->pot == nullptr) {
        // Se è nullo, restituisci un messaggio di errore
        throw std::runtime_error("Pointer phi for " + this->Name + " is null.");
    }

    return 2.*b*(pow( std::abs(1.- (pow(b,2) / pow(r,2)) - (this->pot->Pot(r)/E)), 
        -0.5 ))/pow(r,2);
    
}

// Adatta step di integrazione
std::tuple<double, bool> AdaptDeflAngle::AdaptStep(
        double dx1, double F_bar_x1, double F_x1)       {
    
    const double Tol_rel = 1e-3;

    // restituita per c.c. formula 5 COLONNA
    bool control = true; 

    // err relativo formula 5 COLONNA
    double eps = std::abs((F_bar_x1 - F_x1) / F_x1);
    
    // formula 4 COLONNA n = 1
    double q = 0.9 * std::pow(std::abs(Tol_rel / eps), 2);
    if (q < 0.5) q = 0.5;
    if (q > 1.5) q = 1.5;
    
    dx1 *= q;
    if (eps > Tol_rel) {
    
        // formula 5 non sotto soglia con step attuale
        control = false;
    }
    
    return std::make_tuple(dx1,control);

}

// estremi integrale formula 1 COLONNA
double AdaptDeflAngle::AnalyzeIntegrandKiMin( double b, double E) {
    
    try {
    
        double Tol_dx = 1.e-06;
        double rmin = 1.e-2 ;

        double XIN;

        if (b>1.)
            XIN = b*5.;
        else 
            XIN = 5.;
        
        double XFIN = 0.;
        double dxIN = 0.1;
        const double segno = -1.;
        double x0 = XIN;
        double dx1 = dxIN;
        
        double F_x0 = IntegrandKi(x0,b,E);

        std::vector<double> p {0.,F_x0};
        
        while (segno*x0 < segno*XFIN){
            
            bool control = false;
            double x1;
            double F_x1;
            
            while (control == false){
                x1 = x0+dx1*segno;
                F_x1 = IntegrandKi(x1,b,E);
                p = interpCoeff({x0 , x1},{F_x0 , F_x1});
                /* Quando f approssimata col punto intermedio 
                non aggiunge precisione ha finito 
                (formula 5 controllo in adaptStep) */
                double Fbar_x01 = Fbar((x1+x0)/2.,p);
                double F_x01 = IntegrandKi((x1+x0)/2.,b,E);
                std::tie(dx1,control) = AdaptStep(dx1,Fbar_x01,F_x01);
            }
            x0 = x1;
            F_x0 = F_x1;
            if (dx1<Tol_dx){
                rmin = x0 ;
                break ;
            }
        
        }
        return rmin;

    } catch(const std::exception& e){

        throw ; 

    }
}

// estremi integrale formula 1
double AdaptDeflAngle::AnalyzeIntegrandKiMax(double b) {
    return b > 1 ? b * 1.e+02 : 1.e+02;
}


// integrazione frattale adattiva come da articolo di COLONNA formula1
double AdaptDeflAngle::FractalKi(double rmin, double rmax, double b, double E) {
    
    try{

        double Tol_I = 1e-05;
        int Lx = 3;
        bool control = true ;
        std::vector<double> x = linspace (rmin, rmax, Lx);

        std::vector<double> ff;
        for (int i = 0; i < Lx; i++)
            ff.push_back(IntegrandKi(x[i],b,E));

        int iOK = 0;
        while (control == true){

            control = false;
            for (int i = iOK; i < Lx  -2; i += 2 ){
        
                double I1 = 0.5*std::abs(x[i] - x[i+2])*(ff[i]+ff[i+2]);
                double I2 = 0.5*std::abs(x[i] - x[i+2])*(0.5*ff[i]+0.5*ff[i+2]+ff[i+1]);
                double EI = std::abs(I2-I1);
        
                if (EI>Tol_I) {
                    
                    control = true;
                    std::vector<double> xin = std::vector<double> ( x.begin(), x.begin() + i + 1) ;
                    std::vector<double> xfin = std::vector<double> ( x.begin() + i + 2, x.end()) ;
                    std::vector<double> ffin = std::vector<double> ( ff.begin(), ff.begin() + i + 1);
                    std::vector<double> fffin = std::vector<double> ( ff.begin() + i + 2, ff.end());
                    double xa = x[i];
                    double xb = x[i+2];
                    std::vector<double> xc {x[i+1]};
                    
                    // concatenazione vettori 
                    x = concatenate(xin, std::vector<double> {0.5*(xa+xc[0])}, xc, std::vector<double>  
                        {0.5*(xc[0]+xb)}, xfin);
                    
                    ff = concatenate(ffin, std::vector<double> {IntegrandKi(0.5*(xa+xc[0]),b,E)} , 
                        std::vector<double> {ff[i+1]} , std::vector<double> {IntegrandKi(0.5*(xc[0]+xb),b,E)} , 
                            fffin);
                    // Lx numero di punti del dominio di integrazione aumenta
                    Lx = x.size();
                    break ;
                
                }
                
                else{
                
                    iOK = i+2;
                
                }
            }
        }

        double Itot = 0;
        #pragma omp parallel for reduction(+:Itot)
        for (int i = 0; i < Lx-2; i+=2)
            Itot += 0.5*std::abs( x[i] - x[i+2] )*( 0.5*ff[i] + 0.5*ff[i+2] + ff[i+1] );
        
        double I2rmax = 2*asin(b/rmax);

        return Itot + I2rmax ;
    
    } catch(const std::exception& e){
    
        throw ; 
    }
}

// Angolo di deflessione di formula 1 COLONNA
double AdaptDeflAngle::deflectionAngle ( double b, double E ) {
    
    try {
    
        double rmin = AnalyzeIntegrandKiMin(b,E) ;
        double rmax = AnalyzeIntegrandKiMax(b) ;        

        double ki = FractalKi(rmin,rmax,b,E) ;
        
        theta = std::abs(std::numbers::pi - ki) ;
        
        return theta ; 
    
    } catch (const std::exception& e ){
       
        throw ;
        return 0. ; 
    }
}

// __________________________________ Implementation AvrgDeflAngle ____________________________

AvrgDeflAngle::AvrgDeflAngle( InteractionInterface* i, Potential* pot ) : AbInitioTcsIntegration(i,pot) {

    InitE() ; 

    // interparticle distance 1/R interval
    ws = logspace ( log10(1. / 1000.), log10(1. / 0.0001), N) ; 

    // potential psi(w) = phi(1/R) 
    v.resize(ws.size());
    #pragma omp parallel for
    for (int i = 0; i < ws.size(); i++)
        v[i] = pot->Pot(r0/ws[i]) ;
    
    ws[0] = 0.;
    
}

double AvrgDeflAngle::deflectionAngle ( double Bs, double Gst ) {
        
    double Y,Yold;
	double c = 0.;
	int i;
	
	Yold = 1. - ((v[0]/e0)/Gst) - (Bs*pow(ws[0],2));
	for (i=1;i<ws.size();i++){
        
        double num = 0.;
        double den = 0.;

 		Y = 1. - ((v[i]/e0)/Gst) - (Bs*pow(ws[i],2.));
    	
        if(Y>0.){	
    	
            num = 2.*(ws[i] - ws[i-1]) ; 
            den = sqrt(Y) + sqrt(Yold) ; 
        	c += num / den ;
    	
        } 
		else {

            num = 2.*(ws[i] - ws[i-1])*sqrt(Yold) ;  
    		den = (Yold - Y) ;
            c += num / den ;
    		c *= 2.*sqrt(Bs) ;

            double cbcbcb = std::numbers::pi - c ;
            return std::numbers::pi - c ;
    	}
    	Yold = Y;
    }

    c = std::numbers::pi - 2.*sqrt(Bs)*c;
	return c;
}
// __________________________________ Implementation AvrgChiIntegrator ____________________________

std::tuple < double, double > AvrgChiIntegrator::AsymptoticBmax ( double Bmax0 , double Gsi ) {
    
    double c, Bs;

    // Tolerance for effective deflection
    const double chiTol = 1e-2 ; 
    // Number of iterations in which the tolerance must be met for asymptotic behavior
    const int asymptoticTail = 20 ;

    double Bs0 = 0.0;
    double c0 = deflectionAngle ( Bs0, Gsi ) ; 
    double BsStep = Bmax0 / N;
    
    double Bmax = Bmax0;
    int count = 0;
    for (Bs = Bs0 + BsStep; Bs < Bmax; Bs += BsStep) {
        
        c = deflectionAngle(Bs, Gsi);

        if ( fabs(c) < chiTol) {
        
            count++;
        
            if (count >= asymptoticTail) {
                
                // The Bmax is that found once the tolerance is met twice
                Bmax = Bs ;
                // Finer step for integration
                BsStep = Bmax / N;

                break;        
        
            }
        
        } else {
        
            count = 0;  // Reset count if condition fails
        
        }
        
    }

    // if ( count == 0 )throw std::runtime_error("No Bmax found within the given range");
    return std::make_tuple(Bmax, BsStep);
}

std::vector<double> AvrgChiIntegrator::IntegrateChi ( std::tuple<double, double> B , double Gsi ) {

    double rm2 = r0 * r0;
    
    double maxx = std::get<0>(B) ; 
    double step = std::get<1>(B) ; 

    std::vector<double> q ( 4, 0. ) ;
    double c0 = deflectionAngle ( 0., Gsi ) ; 
    for ( double Bs = 0; Bs < maxx; Bs += step ) {
        
        double c = deflectionAngle ( Bs, Gsi );

        q[0] += std::numbers::pi * rm2 * ( step ) * (( 1.0 - cos(     c     )) + (1.0 - cos(     c0     ))) / 2. ;
        q[1] += std::numbers::pi * rm2 * ( step ) * (( 1.0 - pow( cos(c), 2 )) + (1.0 - pow( cos(c0), 2 ))) / 2. ;
        q[2] += std::numbers::pi * rm2 * ( step ) * (( 1.0 - pow( cos(c), 3 )) + (1.0 - pow( cos(c0), 3 ))) / 2. ;
        q[3] += std::numbers::pi * rm2 * ( step ) * (( 1.0 - pow( cos(c), 4 )) + (1.0 - pow( cos(c0), 4 ))) / 2. ;

        c0 = c;
    
    }
    return q ;  
}

void AvrgChiIntegrator::Compute() {
        
    Q.resize(E.size(), std::vector<double>(4));

    // Adimensional energy ratio
    Gs.resize(E.size());
    #pragma omp parallel for
    for (int i = 0; i < E.size(); i++)
        Gs[i] = E[i] / e0 ; 

    try {
        
        #pragma omp parallel for
        for (int i = 0; i < Gs.size(); i++) {
           
            std::tuple<double,double> BB = AsymptoticBmax ( 2500., Gs[i] ) ;  
            Q[i] = IntegrateChi( BB , Gs[i] ) ;
        
        }
    
    } catch(const std::exception& e) {

        throw;
    
    }

    computed = true ; 

} ; 

// __________________________________ Implementation AdaptChiIntegrator ____________________________

std::vector<double> AdaptChiIntegrator::IntegrandQ( double b, double E ) {
    
    try {
        
        double f1,f2,f3,f4;

        f1 = 2.*std::numbers::pi * b * ( 1. - pow( cos( deflectionAngle ( b, E ) ) , 1. )) ;
        f2 = 2.*std::numbers::pi * b * ( 1. - pow( cos( deflectionAngle ( b, E ) ) , 2. )) ;
        f3 = 2.*std::numbers::pi * b * ( 1. - pow( cos( deflectionAngle ( b, E ) ) , 3. )) ;
        f4 = 2.*std::numbers::pi * b * ( 1. - pow( cos( deflectionAngle ( b, E ) ) , 4. )) ;

        std::vector<double> q ({f1,f2,f3,f4}) ;
        
        return q ; 

    } catch (const std::exception& e ){
    
        throw ;
    
    }
}

std::vector<double> AdaptChiIntegrator::FractalQ( double bmin, double bmax, double Ei ) {

    try {
    
        double Tol_I_rel = 1.e-02;
        double I_min = 1.e-03;

        int Lx = 3;
        std::vector<std::vector<double>> ff(Lx, std::vector<double>(4, 0.0));
        bool control = true;
        std::vector<double> x = linspace(bmin, bmax, Lx);

        for (int i = 0; i < Lx; i++) {
            ff[i] = IntegrandQ(x[i], Ei);
        }

        // Iterazione frattale
        int iOK = 0;
        while (control) {
            control = false;
            for (int i = iOK; i < Lx - 2; i += 2) {
                double I1 = 0.5 * std::abs(x[i] - x[i + 2]) * (ff[i][0] + ff[i + 2][0]);
                double I2 = 0.5 * std::abs(x[i] - x[i + 2]) * (0.5 * ff[i][0] + 0.5 * ff[i + 2][0] + ff[i + 1][0]);
                double EI_rel = std::abs((I2 - I1) / I2);

                if (EI_rel > Tol_I_rel && I2 > I_min) {
                    control = true;

                    std::vector<double> xin(x.begin(), x.begin() + i + 1);
                    std::vector<double> xfin(x.begin() + i + 2, x.end());

                    double xa = x[i], xb = x[i + 2], xc = x[i + 1];
                    std::vector<double> f_ac = IntegrandQ(0.5 * (xa + xc), Ei);
                    std::vector<double> f_cb = IntegrandQ(0.5 * (xc + xb), Ei);

                    x = concatenate(xin, {0.5 * (xa + xc), xc, 0.5 * (xc + xb)}, xfin);

                    std::vector<std::vector<double>> ff_in(ff.begin(), ff.begin() + i + 1);
                    std::vector<std::vector<double>> ff_fin(ff.begin() + i + 2, ff.end());

                    ff = concatenate(ff_in, {f_ac, ff[i + 1], f_cb}, ff_fin);
                    Lx = x.size();
                    break;
                } else {
                    iOK = i + 2;
                }
            }
        }

        // Integrazione finale
        std::vector<double> q(4, 0.0);
        for (int i = 0; i < Lx - 2; i += 2) {
            for (int j = 0; j < 4; j++) {

                if (std::isnan(ff[i][j]) || std::isnan(ff[i + 2][j]) || std::isnan(ff[i + 1][j]))
                    continue;
                
                q[j] += 0.5 * std::abs(x[i] - x[i + 2]) *
                        (0.5 * ff[i][j] + 0.5 * ff[i + 2][j] + ff[i + 1][j]);
            }
        }

        return q;

    } catch (const std::exception& e) {
        
        throw; 
    
    }
}

void AdaptChiIntegrator::Compute(){
    
    Q.resize(E.size(), std::vector<double>(4));
    
    try {
        #pragma omp parallel for
        for (int i = 0; i < E.size(); i++) 
            Q[i] = FractalQ( 1.e-3 , 50., E[i] ) ;       
    
    } catch(const std::exception& e) {
        throw;
    }

    computed = true;

}

// __________________________________ Implementation DcsLoader ____________________________

// Constructor
DcsLoader::DcsLoader(InteractionInterface* interaction) : 
    CsHolder(interaction) {

    Qe = new ElasticLoader() ; 
    Qin = new InelasticLoader() ;
    Name = interaction->InteractionName();
    Init();

}

DcsLoader::DcsLoader( const std::string& prefix, InteractionInterface* i ) : 
    CsHolder(i) { 
    
    customPrefix = prefix + "_" ; 
    Init(); 

}

// Constructor
DcsLoader::DcsLoader(InteractionInterface* interaction, const std::string& hyperref) : 
    CsHolder( interaction ){
        
    DownloadAndSaveData(hyperref);
    Init();

}

void DcsLoader::IntegrateDifferentialCrossSections() {

    Qe->Q.resize(Qe->E.size(), std::vector<double>(4));
    Qin->Q.resize(Qin->E.size(), std::vector<double>(4));

    if (anglesElastic.empty() || sigmaElastic.empty()) {
        // Se non ci sono dati Elastic, imposta Qe->Q a zero
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Qe->E.size(); ++i) {
            for (int j = 0; j < 4; ++j) {
                Qe->Q[i][j] = 0.0;
            }
        }
    } else {
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Qe->E.size(); ++i) {
            for (int j = 0; j < 4; ++j) {

                std::vector<double> integrand(anglesElastic.size());

                for (int k = 0; k < anglesElastic.size(); ++k) {

                    integrand[k] = sigmaElastic[k][i] *
                        (pow(1. - std::cos(anglesElastic[k] * std::numbers::pi / 180.0), j+1)) *
                            std::sin(anglesElastic[k] * std::numbers::pi / 180.0);

                }

                Qe->Q[i][j] = trapz(anglesElastic, integrand);

            }
        }
    }


    if (anglesInlastic.empty() || sigmaInelastic.empty()) {
        // Se non ci sono dati Inelastic, imposta Qin->Q a zero
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Qin->E.size(); ++i) {
            for (int j = 0; j < 4; ++j) {
                Qin->Q[i][j] = 0.0;
            }
        }
    } else {
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Qin->E.size(); ++i) {
            for (int j = 0; j < 4; ++j) {

                std::vector<double> integrand(anglesInlastic.size());

                for (int k = 0; k < anglesInlastic.size(); ++k) {

                    integrand[k] = sigmaInelastic[k][i] *
                        (pow(1. - std::cos(anglesInlastic[k] * std::numbers::pi / 180.0), j+1)) *
                            std::sin(anglesInlastic[k] * std::numbers::pi / 180.0);

                }

                Qin->Q[i][j] = trapz(anglesInlastic, integrand);

            }
        }
    }
}

void DcsLoader::ParseFile(std::ifstream& file) {
    
    std::regex elasticRegex(R"(.*Elastic.*)");
    std::regex inelasticRegex(R"(.*Inelastic.*)");
    std::string line;

    bool isElastic = false;
    bool isInelastic = false;
    bool elasticFound = false; 
    bool inelasticFound = false; 

    while (std::getline(file, line)) {
        
        if (std::regex_match(line, elasticRegex)) {



            isElastic = true;
            isInelastic = false;

            elasticFound = true;
            continue;
        
        }

        if (std::regex_match(line, inelasticRegex)) {
        
            isElastic = false;
            isInelastic = true;
            inelasticFound = true;
            continue;
        
        }

        
        if (line.find("Energies[eV]") != std::string::npos) {
       

       
            std::istringstream iss(line.substr(line.find(",") + 1));
            std::vector<double> energies;
       
            std::string value;

            while (std::getline(iss, value, ',')) 
                energies.push_back(std::stod(value));

            if (isElastic) 
                Qe->E = energies; 
            else if (isInelastic) 
                Qin->E = energies;
            
            continue;
        }

        E = Qe->E ;
        
        if (line.find("Angle[deg]") != std::string::npos) {


        
            while (std::getline(file, line)) {
              


                if (line.empty() || std::regex_match(line, elasticRegex) || std::regex_match(line, inelasticRegex)) {
               

               
                    file.seekg(-static_cast<int>(line.length()) - 1, std::ios_base::cur);
                    break;
            

                }

                std::istringstream iss(line);
                std::string value;

                
                std::getline(iss, value, ',');
                
                double angle = std::stod(value);


                

                std::vector<double> sigmaRow;
                
                while (std::getline(iss, value, ',')) {
                
                    if (value.empty()) 
                        sigmaRow.push_back(9999999999);
                
                    else 
                        sigmaRow.push_back(std::stod(value));
                }

                if (isElastic) {
                    anglesElastic.push_back(angle);
                    sigmaElastic.push_back(sigmaRow);
                } else if (isInelastic) {
                    anglesInlastic.push_back(angle);
                    sigmaInelastic.push_back(sigmaRow);
                }
            }
        }
    }


    if (!elasticFound) {
        Qe->E.clear() ;
        Qe->E.push_back(0.) ;
    }


    if (!inelasticFound) {
        Qin->E.clear();
        Qin->E.push_back(0.);
    }


    if (elasticFound) {
        std::vector<std::pair<int, int>> invalidPositions;
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < sigmaElastic.size(); ++i) {
            for (int j = 0; j < sigmaElastic[0].size(); ++j) {
                if (sigmaElastic[i][j] == 9999999999) {
                    invalidPositions.emplace_back(i, j);
                }
            }
        }
        FixInvalidValues(Qe->E, anglesElastic, sigmaElastic, invalidPositions);
    }


    if (inelasticFound) {
        std::vector<std::pair<int, int>> invalidPositionsInelastic;
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < sigmaInelastic.size(); ++i) {
            for (int j = 0; j < sigmaInelastic[0].size(); ++j) {
                if (sigmaInelastic[i][j] == 9999999999) {
                    invalidPositionsInelastic.emplace_back(i, j);
                }
            }
        }
        FixInvalidValues(Qin->E, anglesInlastic, sigmaInelastic, invalidPositionsInelastic);
    }
}

void DcsLoader::FixInvalidValues( std::vector<double>& energies, 
    std::vector<double>& angles, std::vector<std::vector<double>>& sigma, 
        const std::vector<std::pair<int, int>>& invalidPositions) {

    for (const auto& pos : invalidPositions) {

        int row = pos.first;
        int col = pos.second;

        std::vector<double> validRow, validColumn, validAngles, validEnergies;

        
        for (int j = 0; j < sigma[row].size(); ++j) {
        
            if (std::find(invalidPositions.begin(), invalidPositions.end(), std::make_pair(row, j)) == invalidPositions.end()) {
        
                validRow.push_back(sigma[row][j]);
                validEnergies.push_back(energies[j]);
        
            }
        }

        
        for (int i = 0; i < sigma.size(); ++i) {
        
            if (std::find(invalidPositions.begin(), invalidPositions.end(), std::make_pair(i, col)) == invalidPositions.end()) {
        
                validColumn.push_back(sigma[i][col]);
                validAngles.push_back(angles[i]);
        
            }
        }

        try {
        
            double pcol = interpolateSpline(validEnergies, validRow, energies[col]);
            double prow = interpolateSpline(validAngles, validColumn, angles[row]);
            sigma[row][col] = std::max({pcol, prow});
        
        } catch (...) {
        
            sigma.erase(sigma.begin() + row);
            angles.erase(angles.begin() + row);
            --row;
        
        }
    }
}

void DcsLoader::Compute() {
    
    Init();

    IntegrateDifferentialCrossSections() ; 

    E = Qe->E;
    Q = Qe->Q;

    computed = true ; 

} ;

std::string DcsLoader::BuildFileName ( const std::string& name) { 

    return "DCS_"+customPrefix+name+".csv" ; 

};


void DcsLoader::Init() { if(!loaded) LoadData("Differential_Cross_Sections",Name); }


// ______________________ Implementation ChargeTransferCs _____________________

ChargeTransferCs::ChargeTransferCs( InteractionInterface* i, double A, double B ) : A(A), B(B) {
    
    Name = i->InteractionName() ;
    E = logspace ( log10(1.15e-3), log10(433), 50 ) ; 
    
    double mi = i->GetSp1()->getMass() ;
    double mj = i->GetSp2()->getMass() ;

    mu = ( mi * mj ) / ( mi + mj ) ; 
    mu *= ( 1. / amuKg ) ; 

};


void ChargeTransferCs::Compute() {
    
    Q.resize(E.size(), std::vector<double>(4));

    #pragma omp parallel for
    for ( int i = 0; i < E.size(); i++ ) {
        
        double gij = sqrt ( ( 8.*(E[i]/amuKg)*qe ) / ( std::numbers::pi * mu ) ) ; 
        
        for ( int j = 0; j < 5; j++ )
            Q[i][j] = (1./2.) * pow ( A - B*log(gij) , 2.) ; 

    }

    computed = true;

}
