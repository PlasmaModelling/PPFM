/* 

    MAIN FILE USED IN THE VALIDATION PROCESS 
    DESCRIBED IN THE PURE ARGON VALIDATION FOR THE +
    ZHANG, MURPHY MODULE DESCRIBED ON 
    
    "PPFM (Plasma Properties For Many): An Object Oriented C++ 
    Library for Computing Thermodynamic and Transport Properties
    of Plasmas Under Different Operating Conditions" 
    
    by
    
    A.Vagnoni, E.Ghedini, M.Gherardi

*/

#include "GasMixture.h"
#include "CiBox.h"
#include "PartitionFunction.h"
#include "Devoto.h"
#include "Thermodynamics.h"
#include "ZhangMurphyTP.h"
#include "Potential.h"
#include "CollisionIntegral.h"

int main() {

    // Output folder
    std::string folder = "PureAr4Species";

    // GasMixture definition
    GasMixture* mix = new GasMixture ( 
        
        300.            , 101325.       , 
        new Argon       , new ArgonI    , new ArgonII           ,
        new Electron 
        
    );

    // Non equilibrium parameters to loop on
    std::vector<double> NEparam = { 1. , 2. , 3., 5., 7., 10. };

    // Set LTE T range to compute on 
    const std::vector<double> T = arange ( 300., 30100., 100. );

    // Editable partition function box via Composition accessor
    auto pp = mix->getCompositionObj()->getPfBox();

    (*pp)[0]->setAbInitio();
    (*pp)[1]->setAbInitio();
    (*pp)[2]->setAbInitio();

    mix->getCompositionObj()->getPfBox()->info();

    CiBox cibox (mix) ; 

    // Ar - Ar
    cibox[0]->Pot( new HFDTCS2_ArAr() ) ;   

    // Ar - Ar+ 
    cibox[1]->Load(false);
    auto arari = cibox[1]->GetIntInterface();

    // Elastic + Inelastic collision
    cibox[1]->TCScalculator = new CsHolder (

        // Elastic integration of the potential 
        new MultiCs (
            arari, {
                new AvrgChiIntegrator ( arari, new Morse3Param(1.34,1.69,2.43) ),
                new AvrgChiIntegrator ( arari, new Morse2Param(369.,    2.031) ),
                new AvrgChiIntegrator ( arari, new Morse3Param(0.21,1.63,3.08) ),
                new ThresholdCs (
                    arari, {
                        new AvrgChiIntegrator ( arari, new Morse2Param(2.68e+6,5.889) ) ,
                        new AvrgChiIntegrator ( arari,new Morse2Param(29100,4.154) )
                        },
                        {
                            10.
                        }
                    ),
                new AvrgChiIntegrator ( arari, new Morse3Param(0.10,1.79,3.16) ),
                new AvrgChiIntegrator ( arari, new Morse2Param(1.65e+4,3.88) ),
                
            },
            {
                1.,
                1.,
                1.,
                1.,
                1.,
                1.
            }
        ),
        // Two different Charge Transfer Cross Sections for the two different states
        new MultiCs ( 
            arari, { 

                new ChargeTransferCs ( arari, 26.39, 1.12 ) , 
                new ChargeTransferCs ( arari, 18.96, 0.83 ) 
            
            },
            // States degeneracies for Charge Transfers            
            {
                2.,
                4.
            }
        )
    );

    // Pure Argon polarizability
    double alphaAr = 1.641100; // Ang^3
      
    // Ar - Ar+2 
    cibox[2]->Pot( new Polarization ( new Argon, new ArgonII, alphaAr) );

    // Ar - e-
    cibox[3]->Load(false);    
    cibox[3]->LoadElastic(); 

    cibox.info() ; 

    ZhangTpCsv zhangmurphy ( &cibox , folder ) ; 

    std::vector<double> Th = T ;

    for (size_t i = 0; i < NEparam.size(); i++) {

        Th = T ; 
        for (size_t j = 0; j < T.size(); j++)
            Th[j] /= NEparam[i] ;  
        
        // int to write filenames
        int str = NEparam[i]; 

        // Set the NEparam
        mix->theta->set(NEparam[i]) ; 
        mix->setT(Th[0]) ; 
        mix->restartComposition() ;

        // Computing and printing data
        zhangmurphy.Print ( "ZM2013_Theta" + std::to_string(str), Th, mix );
        
    }
    
    return 0 ;

}
