 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include"TransportCrossSection.h"

std::string TransportCrossSectionCsv::BuildFileName(const std::string& filename ) const  {

    return "TCS_" + filename + ".csv";

}

void TransportCrossSectionCsv::PrepareHeader()  {

    header = "E [eV], Qe(1) [Å²], Qe(2) [Å²], Qe(3) [Å²], Qe(4) [Å²] " ;
    
}

void TransportCrossSectionCsv::PrepareInelasticHeader() {

    header = "E [eV], Qin(1) [Å²], Qin(2) [Å²], Qin(3) [Å²], Qin(4) [Å²] " ; 
    
}

void TransportCrossSectionCsv::PrepareMultiHeader(MultiCs* cscalc, int i) {

    header = "state " + std::to_string(i + 1) + ",E [eV], Q(1) [Å²], Q(2) [Å²], Q(3) [Å²], Q(4) [Å²]";

}

void TransportCrossSectionCsv::PrintMessage(const std::string& filename)  {

    std::cout << "Transport Cross Sections " << 
        tcs->GetIntInterface()->InteractionName() << " printed. \n";

}

void TransportCrossSectionCsv::PrepareData(const std::vector<double>& x, GasMixture* gasmix) {

    std::vector<double> e = tcs->TCScalculator->E ; 
    std::vector<std::vector<double>> q = tcs->TCScalculator->Q ; 

    /* In case of CsHolder the Elatic Cross Sections are stored in the E and Q of the whole class */
    data.resize(e.size(),std::vector<double>(5)) ; 
    for (int i = 0; i < e.size(); i++) {

        data[i][0] = e[i] ; 
        for (int j = 0; j < 4; j++)
            data[i][j+1] = q[i][j] ; 

    }
    
}

void TransportCrossSectionCsv::PrepareInelasticData ( CsHolder* tcsElIn ) {

    std::vector<double> e = tcsElIn->Qin->E ;
    std::vector<std::vector<double>> q = tcsElIn->Qin->Q ; 

    data.resize(e.size(),std::vector<double>(5)) ; 
    for (int i = 0; i < e.size(); i++) {

        data[i][0] = e[i] ; 
        for (int j = 0; j < 4; j++)
            data[i][j+1] = q[i][j] ; 

    }

}

void TransportCrossSectionCsv::PrepareData(MultiCs* cscalc, int i) {
    
    auto cs = (*cscalc)[i];
    if (!cs) return;

    const std::vector<double>& E = cs->E;
    const std::vector<std::vector<double>>& Q = cs->Q;

    data.clear();
    data.resize(E.size(), std::vector<double>(6, 0.0));  // state, E, Q1..Q4

    for (int row = 0; row < E.size(); ++row) {
        data[row][0] = 0.0;            // state N placeholder
        data[row][1] = E[row];         // energia
        for (int j = 0; j < 4; ++j) {
            data[row][2 + j] = Q[row][j];
        }
    }
}

void TransportCrossSectionCsv::Print ( const std::string& filename, const std::vector<double>& x, GasMixture* gasmix )  {

    
    /* if TCScalculator is null then Collision Integral is
    computed with methods not directly referring to Transport Cross Sections 
    this function will exit without doing anything */
    if (tcs->TCScalculator == nullptr) {
        /* std::cerr << " TCScalculator is null for "
            << tcs->GetIntInterface()->InteractionName() <<
            "Collision Integral handles without referring to Transport Cross Sections."<< 
            std::endl; */
        return;
    }
    
    if ( tcs->TCScalculator->Q.empty() )
        tcs->TCScalculator->Compute();

    std::string tempCustomFolder = customFolder;

    customFolder += "/TransportCrossSections_" + tempCustomFolder;

    // MultiCs case
    if (auto p = dynamic_cast<MultiCs*>(tcs->TCScalculator)) {

        for (int i = 0; i < p->Size(); ++i) {

            PrepareMultiHeader(p, i);

            PrepareData(p, i);
            
            if (i == 0)
            
                WriteData(filename);     // sovrascrive
            
            else
            
                AppendData(filename);    // aggiunge
        }

    } else if (auto p = dynamic_cast<CsHolder*>(tcs->TCScalculator)) {

        PrepareHeader();
        PrepareData(x, gasmix);
        WriteData(filename);

        PrepareInelasticHeader();
        PrepareInelasticData(p);
        AppendData(filename);

    } else {

        DataPrinter::Print(filename,x,gasmix) ; 

    }

    PrintMessage(filename);

    customFolder = tempCustomFolder;
}

TransportCrossSectionCsv::TransportCrossSectionCsv(TcsInterface* _tcs) : tcs(_tcs) { }
