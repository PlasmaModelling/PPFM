/* 

    LIGHT REGRESSION TEST ON ArH2 PLASMA.
    FOR TESTING PURPOSES ONLY, NOT FOR PRODUCTION.

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

#include <string>
#include <vector>

#include "csv_compare.hpp"
#include <iostream>

int main() {

    // Output folder used by CTest
    std::string folder = "Test/ArH2_light";

    // GasMixture definition
    GasMixture* mix = new GasMixture ( 
        5000.       , 101325.,
        new Argon   , new ArgonI,
        new Hydrogen, new HydrogenI,
        new Electron 
    );

    // Single-temperature light regression case
    std::vector<double> T = { 5000. };

    // Editable partition function box
    auto editableQbox = mix->getCompositionObj()->getPfBox();

    // Regression target 1: edited Ar partition function
    (*editableQbox)[0]->setAbInitio();

    // Editable Collision Integrals Box
    CiBox cibox(mix);

    // Regression target 2: edited Ar-Ar interaction
    cibox[0]->TCScalculator = new AdaptChiIntegrator(
        cibox[0]->GetIntInterface(), 
        new HFDTCS2_ArAr()
    ); 

    // Output modules initialization
    Thermodynamics thSolver;
    ThermodynamicsCsv th(&thSolver, folder + "/Thermodynamics");

    DevotoTpCsv Dev(new DevotoLteLambdaR(&cibox), folder + "/Devoto"); 
    ZhangTpCsv  zm(new ZhangMurphyTP(&cibox), folder + "/Zhang"); 

    // Print edited partition functions and collision integrals
    editableQbox->PrintPartitionFunctions(T, mix, folder + "/Partition Functions");
    cibox.PrintCollisionIntegrals(T, mix, folder + "/Collision Integrals");
     
    // Compute and print thermodynamic and transport data
    th.Print("", T, mix); 
    Dev.Print("Dev1966", T, mix);
    zm.Print("ZM2013", T, mix);

    bool ok = true;

    const std::string projectRoot = PPFM_PROJECT_ROOT;

    const std::string computed =
        projectRoot + "/out/Test/ArH2_light";

    const std::string reference =
        projectRoot + "/tests/references/ArH2_light";

    ok &= compareCsvFiles(
        computed + "/Thermodynamics/TH_.csv",
        reference + "/Thermodynamics/TH_.csv"
    );

    ok &= compareCsvFiles(
        computed + "/Devoto/TP_Dev1966.csv",
        reference + "/Devoto/TP_Dev1966.csv"
    );

    ok &= compareCsvFiles(
        computed + "/Zhang/TP_ZM2013.csv",
        reference + "/Zhang/TP_ZM2013.csv"
    );

    ok &= compareCsvFiles(
        computed + "/Partition Functions/PF_Ar.csv",
        reference + "/Partition Functions/PF_Ar.csv"
    );

    ok &= compareCsvFiles(
        computed + "/Collision Integrals/CI_Ar_Ar.csv",
        reference + "/Collision Integrals/CI_Ar_Ar.csv"
    );

    delete mix;

    return ok ? 0 : 1;

}