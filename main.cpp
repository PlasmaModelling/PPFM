#include "GasMixture.h"
#include "CiBox.h"
#include "PartitionFunction.h"
#include "Devoto.h"
#include "Thermodynamics.h"
#include "ZhangMurphyTP.h"
#include "Potential.h"
#include "CollisionIntegral.h"
#include "DataPrinter.h"

int main() {

    // Output folder
    std::string folder = "Hg/CollisionIntegrals";

    std::vector<double> T = arange(300., 50100., 100.);

    auto Ci = new CollisionIntegralCsv(new CollisionIntegral(new Argon, new Mercury), folder) ;
    Ci->solver->MultiPot({
            new ArHgX0(),
            new Morse3Param(0.008319, 1.09, 4.47) , 
            new Morse3Param(0.044919, 1.54, 3.38) 
        },{
            1.,
            3.,
            6.
        }
    );
    
    // Ci->Print("Ar-Hg", T, nullptr) ;
    
    Ci->solver->Pot(new ArHgX0());
    // Ci->Print("Ar-Hg.n6+Morse", T, nullptr) ;

    Ci->solver->Pot(new Morse3Param(0.008319, 1.09, 4.47) );
    // Ci->Print("Ar-Hg.Morse2", T, nullptr) ;

    Ci->solver->Pot(new Morse3Param(0.044919, 1.54, 3.38) );
    // Ci->Print("Ar-Hg.Morse3", T, nullptr) ;

    std::vector<double> Tstar = {0.004,0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1., 2., 4., 10., 20., 40., 100., 200.};
    for (size_t i = 0; i < Tstar.size(); i++)
        Tstar[i] *= 100. ; 
    
    Ci->solver->Pot(new Morse3Param(0.00861733,1.,4.693147));
    Ci->Print("test_on_Morse",Tstar,nullptr);
  
} 
 