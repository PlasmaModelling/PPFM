/* 

    MAIN FILE USED IN THE VALIDATION PROCESS 
    OF KRIPTON, XENON AT DIFFERENT PRESSURES.

    REFERENCE ARTICLE:
    " Thermodynamic properties and transport
    coefficients of arc lamp plasmas: argon,
    krypton and xenon "

    BY : Anthony B Murphy and Eugene Tam
    DOI : 10.1088/0022-3727/47/29/295202

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
    std::string folder = "KRIPTON_XENON" ;

    double atm = 101325; // Pa

    std::vector<double> pressures = { atm, 2.*atm, 5.*atm, 10.*atm, 20.*atm, 50.*atm, 100.*atm } ; 

    std::vector<std::string> p_strings = { "1 atm","2 atm","5 atm","10 atm","20 atm","50 atm","100 atm" } ; 

    std::vector<double> T = arange ( 300., 30100., 100. ) ;

    // GasMixture definition
    GasMixture* krypton = new GasMixture ( 
        
        300.            , atm       , 
        new Krypton, new KryptonI, new KryptonII, 
        new KryptonIII, new KryptonIV, new Electron
    
    );

    GasMixture* xenon = new GasMixture (

        300.            , atm       ,
        new Xenon, new XenonI, new XenonII,
        new XenonIII, new XenonIV, new Electron

    );

    // setting the partition function calculations from electronic configurations

    krypton->Comp->Qbox->AllAbInitio() ; 

    xenon->Comp->Qbox->AllAbInitio() ; 

    // Thermodynamic modules

    auto thermKr = new ThermodynamicsCsv ( folder ) ; 

    auto thermXe = new ThermodynamicsCsv ( folder ) ; 

    // polarizabilities 

    double alphaKr = 2.498 ; // Ang^3

    double alphaXe = 4.005 ; // Ang^3

    // editable collision integrals boxes

    auto kryptonCI = CiBox ( krypton ) ;  

    auto xenonCI = CiBox ( xenon ) ; 

    // editing Krypton collision integrals

    kryptonCI.info() ; 

    kryptonCI[0]->Pot ( new HFD_B ( 44.44924, 1.1066257, -0.048669, 128.3040, 3947.999, 170000., 7.579694, 1.208 , true ) ) ; 

    kryptonCI[1]->MultiPot(
        {
            new Morse3Param ( 1.147,	1.532,	2.688 ) ,
            new Morse3Param ( 0.198,	1.945,	3.307 ) ,
            new Morse3Param ( 0.0487,	1.93,	4.084 ) ,
            new Morse3Param ( 0.01299,	1.065,	5.585 ) ,
            new Morse3Param ( 0.135,	1.333,	3.726 ) ,
            new Morse3Param ( 0.0262,	1.211,	4.628 ) 
            
        },
        {
            1.,1.,1.,1.,1.,1.
        }
    );

    kryptonCI[1]->ChargeTransfer ( 26.1, 1.13 ) ; 

    kryptonCI[2]->Pot ( new Polarization ( new Krypton, new KryptonII, alphaKr ) ) ; 
    kryptonCI[3]->Pot ( new Polarization ( new Krypton, new KryptonIII, alphaKr ) ) ; 
    kryptonCI[4]->Pot ( new Polarization ( new Krypton, new KryptonIV, alphaKr ) ) ; 
    kryptonCI[5]->LoadElastic() ; 

    // editing Xenon collision integrals  

    xenonCI.info() ; 

    xenonCI[0]->Pot ( new HFD_B ( 48.72733, 0.9127, -0.049061, 283.900, 11214., 619600., 8.249788, 1.114, true ) );

    xenonCI[1]->MultiPot(
        {
            new Morse3Param ( 0.98 ,	1.368 ,	3.114 ) ,
            new Morse3Param ( 0.236 ,	1.328 ,	3.695 ) ,
            new Morse3Param ( 0.074 ,	1.073 ,	4.395 ) ,
            new Morse3Param ( 0.02 ,	0.72 ,	5.774 ) ,
            new Morse3Param ( 0.199 ,	1.238 ,	3.983 ) ,
            new Morse3Param ( 0.046 ,	0.935 ,	4.773 ) 
        },
        {
            1.,1.,1.,1.,1.,1.
        }
    );
    xenonCI[1]->TCScalculator = new CsHolder ( 
        
        xenonCI[1]->TCScalculator,
        
        new ThresholdCs(xenonCI[1]->GetIntInterface(),
            {
                new ChargeTransferCs ( xenonCI[1]->GetIntInterface(), 78.3, 13.6 ) ,
                new ChargeTransferCs ( xenonCI[1]->GetIntInterface(), 45.7, 8.9 ) 
            },
            {
                10.
            }
        )
    );

    xenonCI[2]->Pot ( new Polarization ( new Xenon, new XenonII, alphaXe ) ) ; 
    xenonCI[3]->Pot ( new Polarization ( new Xenon, new XenonIII, alphaXe ) ) ; 
    xenonCI[4]->Pot ( new Polarization ( new Xenon, new XenonIV, alphaXe ) ) ; 
    xenonCI[5]->LoadElastic() ; 
    
    // Transport modules initialization

    auto DevKr = new DevotoTpCsv ( &kryptonCI, folder ) ;

    auto DevXe = new DevotoTpCsv ( &xenonCI, folder ) ;

    // print cycling on pressures
    for (size_t i = 0; i < pressures.size(); i++) {

        double p = pressures[i] ; 

        krypton->setP(p) ; 
        xenon->setP(p) ; 

        // Composition 
        auto compKr = new CompositionCsv ( krypton, krypton, folder ) ;
        compKr->Qbox->AllAbInitio() ;
        auto compXe = new CompositionCsv ( xenon, xenon, folder ) ;
        compXe->Qbox->AllAbInitio() ;

        compKr->Print ( "Kr Composition " + p_strings[i], T, krypton ) ; 
        compXe->Print ( "Xe Composition " + p_strings[i], T, xenon   ) ; 

        // Thermodynamic properties 
        thermKr->Print ( "Kr Thermodynamics " + p_strings[i], T, krypton ) ; 
        thermXe->Print ( "Xe Thermodynamics " + p_strings[i], T, xenon   ) ; 

        // Transport properties
        DevKr->Print ( "Kr Transport " + p_strings[i], T, krypton ) ; 
        DevXe->Print ( "Kr Transport " + p_strings[i], T, xenon   ) ; 

    }
    
    return 0;

}
