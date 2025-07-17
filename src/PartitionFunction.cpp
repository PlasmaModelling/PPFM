#include"PartitionFunction.h"

std::string PartitionFunctionCsv::BuildFileName(const std::string& filename) const  {
    
    return "PF_" + filename + ".csv";

}

void PartitionFunctionCsv::PrepareHeader()  {

    header = "T [K], P [Pa], λD [m], Q [#]";

}

void PartitionFunctionCsv::PrepareData(const std::vector<double>& Ti, GasMixture* gasmix)  {

    double T0 = gasmix->getTemperature();
    double P = gasmix->getPressure();

    data.resize(Ti.size(), std::vector<double>(4));

    for (int i = 0; i < Ti.size(); ++i) {

        gasmix->setT(Ti[i]);

        double lambda = gasmix->Comp->getDebyeLength(Ti[i]);

        data[i][0] = Ti[i];
        data[i][1] = P;
        data[i][2] = lambda;

        pf->computePartitionFunction(Ti[i], P, lambda);
        data[i][3] = pf->getPf();

    }

    gasmix->setT(T0);
    gasmix->restartComposition();

}

void PartitionFunctionCsv::PrintMessage(const std::string& filename)  {

    std::cout << "Partition Function " << filename <<  " printed. \n" ;

}

PartitionFunctionCsv::PartitionFunctionCsv(PFinterface* _pf) : pf(_pf) {}

void PartitionFunctionCsv::Print ( const std::string& filename, const std::vector<double>& x, GasMixture* gasmix )  {
    
    std::string tempCustomFolder = customFolder ; 

    /*  PartitionFunctions will collect printed partition functions 
    from the reference Mixture in the user called folder */
    customFolder += "/PartitionFunctions_" + tempCustomFolder ;

    DataPrinter::Print(filename,x,gasmix) ; 

    customFolder = tempCustomFolder ; 

}
