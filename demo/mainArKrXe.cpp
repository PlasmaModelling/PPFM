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

    double atm = 101325.; // Pa

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
    krypton->getCompositionObj()->getPfBox()->AllAbInitio(); 
    xenon->getCompositionObj()->getPfBox()->AllAbInitio(); 

    // Thermodynamic solvers + printers
    auto thermArSolver = new Thermodynamics();
    auto thermKrSolver = new Thermodynamics();
    auto thermXeSolver = new Thermodynamics();

    auto thermAr = new ThermodynamicsCsv(thermArSolver, folder);
    auto thermKr = new ThermodynamicsCsv(thermKrSolver, folder); 
    auto thermXe = new ThermodynamicsCsv(thermXeSolver, folder); 

    // Polarizabilities
    double alphaAr = 1.6411; // Å³
    double alphaKr = 2.4844; // 2.498; // Å³
    double alphaXe = 4.004 ; // 4.005; // Å³

    // Collision Integrals
    auto argonCI = CiBox(argon);
    auto kryptonCI = CiBox(krypton);
    auto xenonCI   = CiBox(xenon);

    // Argon interactions
    argonCI.info();

    argonCI[3]->Pot(new Polarization(new Argon, new ArgonIII, alphaAr));
    argonCI[4]->Pot(new Polarization(new Argon, new ArgonIV, alphaAr));

    // Krypton interactions
    kryptonCI.info();

    kryptonCI[0]->Pot(new HFD_B(44.44924, 1.1066257, -0.048669, 128.3040, 3947.999, 170000., 7.579694, 1.208, true));

    kryptonCI[1]->MultiPot(
        {
            new Morse3Param(1.147, 1.532, 2.688),
            new Morse3Param(0.198, 1.945, 3.307),
            new Morse3Param(0.0487, 1.93, 4.084),
            new Morse3Param(0.01299, 1.065, 5.585),
            new Morse3Param(0.135, 1.333, 3.726),
            new Morse3Param(0.0262, 1.211, 4.628)
        },
        { 1., 1., 1., 1., 1., 1. }
    );

    kryptonCI[1]->ChargeTransfer(26.1, 1.13); 
    
    kryptonCI[2]->Pot(new Polarization(new Krypton, new KryptonII, alphaKr)); 
    kryptonCI[3]->Pot(new Polarization(new Krypton, new KryptonIII, alphaKr)); 
    kryptonCI[4]->Pot(new Polarization(new Krypton, new KryptonIV, alphaKr)); 
    kryptonCI[5]->LoadElastic(); 

    // Xenon interactions
    xenonCI.info();

    xenonCI[0]->Pot(new HFD_B(48.72733, 0.9127, -0.049061, 283.900, 11214., 619600., 8.249788, 1.114, true));

    xenonCI[1]->MultiPot(
        {
            new Morse3Param(0.98, 1.368, 3.114),
            new Morse3Param(0.236, 1.328, 3.695),
            new Morse3Param(0.074, 1.073, 4.395),
            new Morse3Param(0.02, 0.72, 5.774),
            new Morse3Param(0.199, 1.238, 3.983),
            new Morse3Param(0.046, 0.935, 4.773)
        },
        { 1., 1., 1., 1., 1., 1. }
    );

    xenonCI[1]->TCScalculator = new CsHolder(
        xenonCI[1]->TCScalculator,
        new ThresholdCs(xenonCI[1]->GetIntInterface(),
            {
                new ChargeTransferCs(xenonCI[1]->GetIntInterface(), 78.3, 13.6),
                new ChargeTransferCs(xenonCI[1]->GetIntInterface(), 45.7, 8.9)
            },
            { 10. }
        )
    );

    xenonCI[2]->Pot(new Polarization(new Xenon, new XenonII, alphaXe)); 
    xenonCI[3]->Pot(new Polarization(new Xenon, new XenonIII, alphaXe)); 
    xenonCI[4]->Pot(new Polarization(new Xenon, new XenonIV, alphaXe)); 
    xenonCI[5]->LoadElastic(); 

    // Transport printers
    auto devAr = new DevotoTpCsv( new DevotoTP(&argonCI), folder);
    auto devKr = new DevotoTpCsv( new DevotoTP(&kryptonCI), folder);
    auto devXe = new DevotoTpCsv( new DevotoTP(&xenonCI)  , folder);

    // Loop on pressures
    for (size_t i = 0; i < pressures.size(); i++) {

        double p = pressures[i]; 

        argon->setP(p);
        krypton->setP(p);
        xenon->setP(p);

        // Composition printers
        auto compAr = new CompositionCsv(argon->getCompositionObj(), folder);
        auto compKr = new CompositionCsv(krypton->getCompositionObj(), folder);
        auto compXe = new CompositionCsv(xenon->getCompositionObj(),   folder);

        compAr->Print("Ar Composition " + p_strings[i], T, argon);
        compKr->Print("Kr Composition " + p_strings[i], T, krypton); 
        compXe->Print("Xe Composition " + p_strings[i], T, xenon); 

        // Thermodynamic properties
        thermAr->Print("Ar Thermodynamics " + p_strings[i], T, argon);
        thermKr->Print("Kr Thermodynamics " + p_strings[i], T, krypton);
        thermXe->Print("Xe Thermodynamics " + p_strings[i], T, xenon);

        // Transport properties
        // devAr->Print("Ar Transport " + p_strings[i], T, argon);
        // devKr->Print("Kr Transport " + p_strings[i], T, krypton); 
        // devXe->Print("Xe Transport " + p_strings[i], T, xenon); 
    }

    return 0;
}
