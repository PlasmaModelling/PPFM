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
    std::string folder = "ARGON_KRIPTON_XENON";

    double atm = 100000.; // Pa

    std::vector<double> pressures = { atm, 2.*atm, 5.*atm, 10.*atm, 20.*atm, 50.*atm, 100.*atm }; 
    std::vector<std::string> p_strings = { "1 atm","2 atm","5 atm","10 atm","20 atm","50 atm","100 atm" }; 
    std::vector<double> T = arange(300., 30100., 100.);

    // GasMixtures

    auto argon = new GasMixture ( 300., atm,
        
        new Argon, new ArgonI, new ArgonII,
        new ArgonIII, new ArgonIV, new Electron
    
    );
    
    auto krypton = new GasMixture ( 300., atm,
        
        new Krypton, new KryptonI, new KryptonII,
        new KryptonIII, new KryptonIV, new Electron
    
    );

    auto xenon = new GasMixture(300., atm,
        new Xenon, new XenonI, new XenonII,
        new XenonIII, new XenonIV, new Electron
    );

    // Partition functions setup
    auto argonPF = argon->getCompositionObj()->getPfBox();
    auto kryptonPF = krypton->getCompositionObj()->getPfBox();
    auto xenonPF = xenon->getCompositionObj()->getPfBox();

    argonPF->AllAbInitio();
    kryptonPF->AllAbInitio(); 
    xenonPF->AllAbInitio(); 

    argon->getCompositionObj()->getPfBox()->info();
    krypton->getCompositionObj()->getPfBox()->info();
    xenon->getCompositionObj()->getPfBox()->info();

    argon->setCompositionSolver(new GTSahaDHcorrection(argon,argon,argonPF));
    krypton->setCompositionSolver(new GTSahaDHcorrection(krypton,krypton,kryptonPF));
    xenon->setCompositionSolver(new GTSahaDHcorrection(xenon,xenon,xenonPF));

    // Polarizabilities NEW
    double alphaAr = 1.6423; // Å³
    double alphaKr = 2.4865; // Å³
    double alphaXe = 4.0484; // Å³

    // Polarizabilities OLD
    // double alphaAr = 1.62; // Å³
    // double alphaKr = 2.46; // Å³
    // double alphaXe = 3.99; // Å³

    // Collision Integrals
    auto argonCI = CiBox(argon);
    auto kryptonCI = CiBox(krypton);
    auto xenonCI   = CiBox(xenon);

    // Argon interactions
    
    // Ar - Ar
    argonCI[0]->Pot( new HFDTCS2_ArAr() ) ;   
    /* argonCI[0]->Pot( new Capitelli ( new Argon, new Argon, 8.1193, 0.0116273, 3.79383 )); */

    // Ar - Ar+ 
    argonCI[1]->Load(false);
    auto arari = argonCI[1]->GetIntInterface();

    // Elastic + Inelastic collision
    argonCI[1]->TCScalculator = new CsHolder (

        // Elastic integration of the potential 
        /* new AdaptChiIntegrator ( arari, new Capitelli ( new Argon, new ArgonI, 8.1193, 0.0236949, 4.70232 )), */
            
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
            { 1./6.,1./6.,1./6.,1./6.,1./6.,1./6. }
        ), 
       
        // Two different Charge Transfer Cross Sections for the two different states
        new MultiCs ( 
            arari, { 

                new ChargeTransferCsWithEnergy ( arari, 8.921, 0.3960 ) , 
                new ChargeTransferCsWithEnergy ( arari, 6.189, 0.2934 ) 
            
            },
            // States degeneracies for Charge Transfers 
            { 1./3., 2./3. }
        )
    );

    // WEIGHT THE ENERGIES AS IN THE ARTICLE FOR THE AR - AR+ ELASTIC COLLISION
    auto ariar = dynamic_cast<MultiCs*>(dynamic_cast<CsHolder*>(argonCI[1]->TCScalculator)->Qe);
    ariar->setEnergyWeights( { 0.,0.17752,0.,0.,0.,0.17752 } );

    // Ar - Ar+2 
    argonCI[2]->Pot( new Polarization ( new Argon, new ArgonII, alphaAr) );
    // Ar - Ar+3
    argonCI[3]->Pot(new Polarization(new Argon, new ArgonIII, alphaAr));
    // Ar - Ar+4
    argonCI[4]->Pot(new Polarization(new Argon, new ArgonIV, alphaAr));
    // Ar - e-
    argonCI[5]->Load(false);
    argonCI[5]->LoadElastic();
    argonCI.info();

    // Krypton interactions
    // Kr - Kr
     kryptonCI[0]->Pot(new HFD_B(6.9735391e+4,8.38802216,-2.79611543,1.06136003,0.56845577,0.42605480,4.011,1.2080,201.3));
    // kryptonCI[0]->Load(true);
     /* kryptonCI[0]->Pot ( new Capitelli (new Krypton, new Krypton, 7.84589, 0.0176239, 4.02566 )) ;  */

    // Kr - Kr+
    kryptonCI[1]->Load(false);
    kryptonCI[1]->TCScalculator = new CsHolder (

   /*      new AdaptChiIntegrator ( kryptonCI[1]->GetIntInterface(), new Capitelli ( new Krypton, new KryptonI, 7.84589, 0.0273176, 5.13097 ) 
        ), */
        
        new MultiCs( kryptonCI[1]->GetIntInterface(),
            {
                new AvrgChiIntegrator ( kryptonCI[1]->GetIntInterface(), new Morse3Param( 1.1358192415424 , 1.53051818584054  , 2.68822023138724) ),
                new AvrgChiIntegrator ( kryptonCI[1]->GetIntInterface(), new Morse3Param( 0.193415349504  , 1.48727664654195  , 3.30735756814375) ),
                new AvrgChiIntegrator ( kryptonCI[1]->GetIntInterface(), new Morse3Param( 0.0469900111936 , 1.12449104378486  , 4.08524806817116) ),
                new AvrgChiIntegrator ( kryptonCI[1]->GetIntInterface(), new Morse3Param( 0.01239841984   , 0.726959187763507 , 5.58811134713568) ),
                new AvrgChiIntegrator ( kryptonCI[1]->GetIntInterface(), new Morse3Param( 0.132043171296  , 1.34117458994349  , 3.73069933686615) ),
                new AvrgChiIntegrator ( kryptonCI[1]->GetIntInterface(), new Morse3Param( 0.0250448080768 , 0.923072806176961 , 4.63030059540125) )
            }, 
            {
                1./6.,1./6.,1./6.,1./6.,1./6.,1./6. 
            } 
        ), 
       
        new ChargeTransferCsWithEnergy ( kryptonCI[1]->GetIntInterface(), 8.6483, 0.3995 )
    );

    // WEIGHT THE ENERGIES AS IN THE ARTICLE FOR THE KR - KR+ ELASTIC COLLISION
    auto krikr = dynamic_cast<MultiCs*>(dynamic_cast<CsHolder*>(kryptonCI[1]->TCScalculator)->Qe);
    krikr->setEnergyWeights( { 0.,0.,0.,0.,0.6659,0.6659 } );

    // Kr - Kr+2
    kryptonCI[2]->Pot(new Polarization(new Krypton, new KryptonII, alphaKr)); 
    // Kr - Kr+3
    kryptonCI[3]->Pot(new Polarization(new Krypton, new KryptonIII, alphaKr)); 
    // Kr - Kr+4
    kryptonCI[4]->Pot(new Polarization(new Krypton, new KryptonIV, alphaKr)); 
    // Kr - e-
    kryptonCI[5]->LoadElastic(); 
    kryptonCI.info();

    // Xenon interactions
    // Xe - Xe
    xenonCI[0]->Pot(new HFD_B(5.44087277e+4,7.52958289,-3.3390428,1.00555220,0.58359858,0.47378306,4.3656,1.114,282.8));
    /* xenonCI[0]->Pot ( new Capitelli(new Xenon, new Xenon, 7.56933, 0.0254689, 4.31627 )); */


    // Xe - Xe+
    xenonCI[1]->Load(false);
    xenonCI[1]->TCScalculator = new CsHolder(
        /* new AdaptChiIntegrator ( xenonCI[1]->GetIntInterface(), new Capitelli ( new Xenon, new XenonI, 7.56933, 0.0334937, 5.64973 ) 
        ), */
            
            new MultiCs( xenonCI[1]->GetIntInterface(),
            {
                new AvrgChiIntegrator ( xenonCI[1]->GetIntInterface(), new Morse3Param ( 0.9722840838528 , 1.36752358722894   , 3.11420788616415 ) ) ,
                new AvrgChiIntegrator ( xenonCI[1]->GetIntInterface(), new Morse3Param ( 0.2325943561984 , 1.32789439985817   , 3.69471528652475 ) ) ,
                new AvrgChiIntegrator ( xenonCI[1]->GetIntInterface(), new Morse3Param ( 0.0721588034688 , 1.07313869362443   , 4.39534591376032 ) ) ,
                new AvrgChiIntegrator ( xenonCI[1]->GetIntInterface(), new Morse3Param ( 0.0195895033472 , 0.720170444102808  , 5.77438172537354 ) ) ,
                new AvrgChiIntegrator ( xenonCI[1]->GetIntInterface(), new Morse3Param ( 0.1955230808768 , 1.29151629657865   , 3.98258768925598 ) ) ,
                new AvrgChiIntegrator ( xenonCI[1]->GetIntInterface(), new Morse3Param ( 0.0445103272256 , 0.935076315085618  , 4.77264926513416 ) )
            },
            {
                1./6.,1./6.,1./6.,1./6.,1./6.,1./6. 
            }
        ), 
       
        new ChargeTransferCsWithEnergy(xenonCI[1]->GetIntInterface(), 9.716, 0.4204)
    );

    // WEIGHT THE ENERGIES AS IN THE ARTICLE FOR THE XE - XE+ ELASTIC COLLISION
    auto xeixe = dynamic_cast<MultiCs*>(dynamic_cast<CsHolder*>(xenonCI[1]->TCScalculator)->Qe);
    xeixe->setEnergyWeights( { 0.,0.,0.,0.,1.3066,1.3066 } );
    
    // Xe - Xe+2
    xenonCI[2]->Pot(new Polarization(new Xenon, new XenonII, alphaXe)); 
    // Xe - Xe+3
    xenonCI[3]->Pot(new Polarization(new Xenon, new XenonIII, alphaXe)); 
    // Xe - Xe+4
    xenonCI[4]->Pot(new Polarization(new Xenon, new XenonIV, alphaXe)); 
    
    std::vector<double> T2 = concatenate({300.,500.,750.},arange(1000.,30001.,1000.));
    
    // Xe - e-
    xenonCI[5]->Load(false);
    xenonCI[5]->LoadElastic();
    xenonCI.info();

    // Thermodynamic solvers + printers
    auto thermArSolver = new ThermodynamicsDHcorrected();
    auto thermKrSolver = new ThermodynamicsDHcorrected();
    auto thermXeSolver = new ThermodynamicsDHcorrected();
    auto thermAr = new ThermodynamicsCsv(thermArSolver, folder +"/UPDATED_SPECIESGASEOUS");
    auto thermKr = new ThermodynamicsCsv(thermKrSolver, folder +"/UPDATED_SPECIESGASEOUS"); 
    auto thermXe = new ThermodynamicsCsv(thermXeSolver, folder +"/UPDATED_SPECIESGASEOUS"); 

    // Transport printers
    auto devAr = new DevotoTpCsv( new DevotoLteLambdaR(&argonCI  , thermAr->solver), folder +"/UPDATED_SPECIESGASEOUS");
    auto devKr = new DevotoTpCsv( new DevotoLteLambdaR(&kryptonCI, thermKr->solver), folder +"/UPDATED_SPECIESGASEOUS");
    auto devXe = new DevotoTpCsv( new DevotoLteLambdaR(&xenonCI  , thermXe->solver), folder +"/UPDATED_SPECIESGASEOUS");

    devAr->solver->setOrders( {2,2,3,2,2,3} );
    devKr->solver->setOrders( {2,2,3,2,2,3} );
    devXe->solver->setOrders( {2,2,3,2,2,3} );

    // Loop on pressures
    for (size_t i = 0; i < pressures.size(); i++) {

        double p = pressures[i]; 

        argon->setP(p);
        krypton->setP(p);
        xenon->setP(p);

        // Composition printers
        auto compAr = new CompositionCsv(argon->getCompositionObj(), folder +"/UPDATED_SPECIESGASEOUS");
        auto compKr = new CompositionCsv(krypton->getCompositionObj(), folder +"/UPDATED_SPECIESGASEOUS");
        auto compXe = new CompositionCsv(xenon->getCompositionObj(),   folder +"/UPDATED_SPECIESGASEOUS");

        // Equilibrium compositions
         compAr->Print("Ar Composition " + p_strings[i], T, argon);
         compKr->Print("Kr Composition " + p_strings[i], T, krypton); 
         compXe->Print("Xe Composition " + p_strings[i], T, xenon); 

        // Thermodynamic properties
         thermAr->Print("Ar Thermodynamics " + p_strings[i], T, argon);
         thermKr->Print("Kr Thermodynamics " + p_strings[i], T, krypton);
         thermXe->Print("Xe Thermodynamics " + p_strings[i], T, xenon);

        // Transport properties
         devAr->Print("Ar Transport " + p_strings[i], T, argon);
         devKr->Print("Kr Transport " + p_strings[i], T, krypton); 
         devXe->Print("Xe Transport " + p_strings[i], T, xenon); 

    }

    return 0;

}