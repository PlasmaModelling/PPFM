/* 

    MAIN FILE USED IN THE VALIDATION PROCESS 
    OF ARGON, KRIPTON, XENON AT DIFFERENT PRESSURES.

    REFERENCE ARTICLE:
    "Thermodynamic properties and transport
    coefficients of arc lamp plasmas: argon,
    krypton and xenon"

    BY: Anthony B Murphy and Eugene Tam
    DOI: 10.1088/0022-3727/47/29/295202

*/

#include "Potential.h"
#include "TransportCrossSection.h"
#include "CollisionIntegral.h"
#include "DataPrinter.h"

int main() {

    // Output folder
    std::string folder = "Extended_Argon_ArgonIon_interaction";

    auto CollIntArArIon = new CollisionIntegral(new Argon, new ArgonI) ; 

    CollIntArArIon->Load( false ) ; // Force to not use the loader 
    auto arariInterface = CollIntArArIon->GetIntInterface() ;

    CollIntArArIon->TCScalculator = new CsHolder (

        // Elastic integration of the potential 
        new MultiCs (

            arariInterface, {
            
                new AvrgChiIntegrator ( arariInterface, new Morse3Param(1.34,1.69,2.43) ),
                new AvrgChiIntegrator ( arariInterface, new Morse2Param(369.,    2.031) ),
                new AvrgChiIntegrator ( arariInterface, new Morse3Param(0.21,1.63,3.08) ),
                new ThresholdCs (
            
                    arariInterface, {
            
                        new AvrgChiIntegrator ( arariInterface, new Morse2Param(2.68e+6,5.889) ) ,
                        new AvrgChiIntegrator ( arariInterface, new Morse2Param(29100,4.154) )
            
                    },
                    {
                            10.
                    }
        
                ),
        
                new AvrgChiIntegrator ( arariInterface, new Morse3Param(0.10,1.79,3.16) ),
                new AvrgChiIntegrator ( arariInterface, new Morse2Param(1.65e+4,3.88) ),

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
            arariInterface, { 

                new ChargeTransferCs ( arariInterface, 26.39, 1.12 ) , 
                new ChargeTransferCs ( arariInterface, 18.96, 0.83 ) 
            
            },
            // States degeneracies for Charge Transfers            
            {
                2.,
                4.
            }
        )
    );

    std::vector<double> T = arange( 300. , 200100., 100. ) ; // K°

    auto CIprinter = new CollisionIntegralCsv ( CollIntArArIon , folder ) ;
    CIprinter->Print( "Ar_ArPlus_Omega" , T , nullptr) ;

    auto TCSprinter = new TransportCrossSectionCsv ( CollIntArArIon->GetTcsInterface() , folder ) ;
    TCSprinter->Print( "Ar_ArPlus_Q" , T, nullptr ) ;
    
    return 0;
}
