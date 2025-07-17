#include "CollisionIntegral.h"

std::string CollisionIntegralCsv::BuildFileName(const std::string& filename ) const  {

    return "CI_" + filename + ".csv";

}

void CollisionIntegralCsv::PrepareHeader()  {

    header = "Tij* [#],"
    " Ω(1;1)[#], Ω(1;2)[#], Ω(1;3)[#], Ω(1;4)[#], Ω(1;5)[#], Ω(1;6)[#], Ω(1;7)[#],"
    " Ω(2;2)[#], Ω(2;3)[#], Ω(2;4)[#], Ω(2;5)[#], Ω(2;6)[#],"
    " Ω(3;3)[#], Ω(3;4)[#], Ω(3;5)[#],"
    " Ω(4;4)[#] ";

}

void CollisionIntegralCsv::PrepareData(const std::vector<double>& x, GasMixture* gasmix )  {

    data.resize(x.size(), std::vector<double>(17));

    for (size_t i = 0; i < x.size(); ++i) {

        gasmix->setT(x[i]);
        
        double lambda = gasmix->Comp->getDebyeLength(x[i]);
        
        double Te = x[i] * gasmix->theta->get();
        
        ci->ComputeCollisionIntegral(Te, x[i], lambda);

        data[i][0] = ci->TijStar ; 

        auto omega = ci->omega4th ;
        
        for (size_t j = 0; j < omega.size(); ++j)
            data[i][j+1] = omega[j];
    
    }

    // SETBACK
    gasmix->setT(x[0]);
    gasmix->restartComposition();
    
}

void CollisionIntegralCsv::PrintMessage(const std::string& filename)  {

    std::cout << "Collision Integral " << 
        ci->InteractionName() << " printed. \n";

}

void CollisionIntegralCsv::Print ( const std::string& filename, const std::vector<double>& x, GasMixture* gasmix )  {

    std::string tempCustomFolder = customFolder ; 

    /*  CollisionIntegrals_ folder will collect every printed collision integral 
    dor the reference mixture */
    customFolder += "/CollisionIntegrals_" + tempCustomFolder ;

    DataPrinter::Print(filename,x,gasmix) ; 

    customFolder = tempCustomFolder ; 

}

CollisionIntegralCsv::CollisionIntegralCsv(CInterface* _ci) : ci(_ci) { }