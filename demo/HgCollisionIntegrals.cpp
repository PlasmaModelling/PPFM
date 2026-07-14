
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
    
    Ci->Print("Ar-Hg", T, nullptr) ;
    
    Ci->solver->Pot(new ArHgX0());
    Ci->Print("Ar-Hg.n6+Morse", T, nullptr) ;

    Ci->solver->Pot(new Morse3Param(0.008319, 1.09, 4.47) );
    Ci->Print("Ar-Hg.Morse2", T, nullptr) ;

    Ci->solver->Pot(new Morse3Param(0.044919, 1.54, 3.38) );
    Ci->Print("Ar-Hg.Morse3", T, nullptr) ;

    Ci->solver = new CollisionIntegral(new Krypton, new Mercury);

    Ci->solver->MultiPot(
        {
            new Morse3Param(0.022069, 1.40341, 4.07) , 
            new Morse3Param(0.077948, 1.51702, 3.35) ,
            new Morse3Param(0.012993, 1.01510, 4.58)
        },{
            1.,
            3.,
            6.
        }
    );
    
    //Ci->Print("Kr-Hg", T, nullptr) ;

    Ci->solver->Pot(new Morse3Param(0.022069, 1.40341, 4.07));
    //Ci->Print("Kr-Hg.Morse1", T, nullptr) ;

    Ci->solver->Pot(new Morse3Param(0.077948, 1.51702, 3.35));
    //Ci->Print("Kr-Hg.Morse2", T, nullptr) ;

    Ci->solver->Pot(new Morse3Param(0.012993, 1.01510, 4.58));
    //Ci->Print("Kr-Hg.Morse3", T, nullptr) ;

    Ci->solver = new CollisionIntegral(new Xenon, new Mercury);
    
    Ci->solver->MultiPot(
        {
            new Morse3Param(0.031492, 1.2456, 4.25) , 
            new Morse3Param(0.171208, 1.58133, 3.15) ,
            new Morse3Param(0.023259, 0.769038, 4.47)
        },{
            1.,
            3.,
            6.
        }
    );

    //Ci->Print("Xe-Hg", T, nullptr) ;
    Ci->solver->Pot(new Morse3Param(0.031492, 1.2456, 4.25));
    //Ci->Print("Xe-Hg.Morse1", T, nullptr) ;
    Ci->solver->Pot(new Morse3Param(0.171208, 1.58133, 3.15));
    //Ci->Print("Xe-Hg.Morse2", T, nullptr) ;
    Ci->solver->Pot(new Morse3Param(0.023259, 0.769038, 4.47));
    //Ci->Print("Xe-Hg.Morse3", T, nullptr) ;

    Ci->solver = new CollisionIntegral(new Mercury, new Mercury);
    
    Ci->solver->Pot(new HulburtHirschfelderUnreduced(328.,3.743,18.4,0.28,1.19e-2,1.7e-4));
    //Ci->Print("Hg-Hg.1", T, nullptr);

    Ci->solver->Pot(new HulburtHirschfelderUnreduced(365.,3.69,19.7,0.29,1.23e-2,1.7e-4));
    //Ci->Print("Hg-Hg.2", T, nullptr);

    Ci->solver->Pot(new Morse3Param(0.0452538594096887,1.25757,3.69));
    //Ci->Print("Hg-Hg.3", T, nullptr);

    Ci->solver = new CollisionIntegral(new Argon, new MercuryI);

    Ci->solver->Pot(new HulburtHirschfelderUnreduced(1650.,2.91,92.9,1.49,0.059,9.37e-4));
    //Ci->Print("Ar-Hg+", T, nullptr); // done

    Ci->solver = new CollisionIntegral(new Krypton, new MercuryI);

    Ci->solver->Pot(new HulburtHirschfelderUnreduced(2867.,2.89,92.,0.8,0.0338,2.85e-4));
    //Ci->Print("Kr-Hg+", T, nullptr); // done 

    Ci->solver = new CollisionIntegral(new Xenon, new MercuryI);
    
    Ci->solver->Pot(new HulburtHirschfelderUnreduced(5237.,2.95,100.7,0.44,0.0243,1.13e-4));
    //Ci->Print("Xe-Hg+", T, nullptr) ; //done 

    Ci->solver = new CollisionIntegral(new ArgonI, new Mercury);
    
    Ci->solver->Pot( new Capitelli( new ArgonI, new Mercury, 7.35501441681418, 0.2091111880553, 3.55758061480179 ));
    //Ci->Print("Hg-Ar+", T, nullptr) ; // done 

    Ci->solver = new CollisionIntegral(new KryptonI, new Mercury);
    
    Ci->solver->Pot( new Capitelli( new KryptonI, new Mercury, 7.2343314858538, 0.193340662665182, 3.71521288508569 ));
    //Ci->Print("Hg-Kr+", T, nullptr) ; // done 

    Ci->solver = new CollisionIntegral(new XenonI, new Mercury);
    
    Ci->solver->Pot( new Capitelli( new XenonI, new Mercury, 7.10327407957034, 0.178989158577641, 3.91064847125162 ));
    //Ci->Print("Hg-Xe+", T, nullptr) ; // done 

    Ci->solver = new CollisionIntegral(new MercuryI, new Mercury);
    
    Ci->solver->Pot( new Capitelli( new MercuryI, new Mercury, 7.09776681990272, 0.178424934962441, 3.91955775867222 ));
    //Ci->Print("Hg-Hg+.el", T, nullptr) ; // elastic done 
    Ci->solver->ChargeTransfer(10.13337,0.382136);
    //Ci->Print("Hg-Hg+", T, nullptr) ; // done

    Ci->solver = new CollisionIntegral(new Mercury, new Electron);
    Ci->solver->LoadElastic();
    //Ci->Print("Hg-e", T, nullptr) ; // done

} 
 