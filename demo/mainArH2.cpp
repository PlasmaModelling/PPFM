/* 

    MAIN FILE USED IN THE VALIDATION PROCESS 
    DESCRIBED IN THE ILLUSTRATIVE EXAMPLE OF 
    
    "PPFM (Plasma Properties For Many): An Object Oriented C++ 
    Library for Computing Thermodynamic and Transport Properties
    of Plasmas Under Different Operating Conditions" 
    
    by
    
    A.Vagnoni, E.Ghedini, M.Gherardi

    DOI : 

    The scheme for transport cross sections and 
    collision integrals computation described in the tables of

    "Two-temperature thermodynamic and transport properties
    of argon–hydrogen and nitrogen–hydrogen plasmas"
    
    by 

    V Colombo, E Ghedini and P Sanibondi

    DOI : doi:10.1088/0022-3727/42/5/055213

    Has been faithfully reproduced.

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
std::string folder = "ArH2_NLTE";

// GasMixture definition
GasMixture* mix = new GasMixture ( 

    500.            , 101325.       , 
    new Argon       , new ArgonI    , new ArgonII           ,
    new Hydrogen    , new HydrogenI , new MolecularHydrogen ,
    new Electron 

);

// Non equilibrium parameters to loop on
std::vector<double> NEparam = { 1. , 2. , 3. };

// Set temperature range to compute on 
std::vector<double> T = arange ( 500., 30100., 100. );

// Setting editable partition function box
auto editableQbox = mix->getCompositionObj()->getPfBox();

// Get info on partition functions for editing
editableQbox->info(); 

// Set calculation from electronic configuration of H partition function
(*editableQbox)[3]->setAbInitio();

// Set desired mixture molar fraction
mix->setMoleFractions ( { 0.5 } , new MolecularHydrogen ); 

// Initialize editable Collision Integrals Box
CiBox cibox (mix);

// Get info on collision integrals for editing
cibox.info();

// Species polarizabities
double alphaAr = 1.641100; // Ang^3
double alphaH  = 0.666793; // Ang^3
double alphaH2 = 0.804225; // Ang^3


// Ar - Ar 
cibox[0]->TCScalculator = new AvrgChiIntegrator (

    // Manually setting an AvrgChiIntegrator for performances
    cibox[0]->GetIntInterface() , 
    new HFDTCS2_ArAr() 

); 

// Same result is obtained doing : 
// cibox[0].Pot(new HFDTCS2_ArAr())
// as algorithm is selected recognizing potentials 

// Ar - Ar+ 
cibox[1]->Load(false);
auto arari = cibox[1]->GetIntInterface();

// CsHolder for Elastic + Inelastic collision
cibox[1]->TCScalculator = new CsHolder (

    // Elastic integration of the potential 
    new AdaptChiIntegrator ( arari, 

        new HulburtHirschfelderUnreduced ( 

            // eps , re       , we      , weXe, Be   , Alphae
            0.01136, 3.759e-10, 3.068e-3, 256., 5.973, 0.39 

        ) 
    ) ,

    // Two different Charge Transfer CS for the inelastic part
    new MultiCs ( 
        arari, 

        // Cross - sections calculators
        { 

            new ChargeTransferCs ( arari, 26.39, 1.12 ), 
            new ChargeTransferCs ( arari, 18.96, 0.83 ) 

        } , 

        // State degeneracies
        {
            2.,
            4.
        }
    ) 

);

// Ar - Ar+2 
cibox[2]->Pot( new Polarization ( new Argon, new ArgonII, alphaAr ) );

// Ar - H 
cibox[3]->Pot( new Morse2Param ( 138.3 , 2.783 ) );

// Ar - H+
cibox[4]->Pot( new Morse3Param ( 4.1 , 1.792 , 1.33 ) );

// Ar - H2 
cibox[5]->Pot( new Morse3Param ( 7.019e-3 , 1.895 , 3.615 ) );

// Ar - e-
cibox[6]->Load(false);
InteractionInterface* Are = cibox[6]->GetIntInterface();

// Different Cross Sections for different energy ranges 
cibox[6]->TCScalculator = new ThresholdCs ( Are, 
    
    {
        // Quantum phase shifts from literature
        new PhaseShiftsLoader ( "Bell1984" , Are ), 
        new PhaseShiftsLoader ( "Gibson1996" , Are ),

        // Differential cross sections
        new DcsLoader(Are) 
    },
    
    {
        1., // Bell until 1.eV
        10. // Gibson until 10.eV, 
            // Dcs after.
    }

);

// Ar+ - H 
// Multi-state interaction with different potentials
cibox[9]->MultiPot ( 
    { 
        //                    Re,      D,    B0,   gamma,  lambda
        new Morse5Param ( 1.3229, 2.1806, 1.509, -0.4788, 0.13659 ),
        new Morse5Param ( 1.3229, 2.2873, 2.118, 0.21770, 0.21200 ),
        new Morse5Param ( 1.3229, 1.8190, 1.613, -0.5157, 0.14108 ),
        new Morse5Param ( 1.3229, 4.0051, 1.080, -0.4395, 0.15825 )
    },
    {
        // statistical degeneracies gk of the states
        2.,
        3.,
        6.,
        1.
    }
); 

// Ar+-H2 
cibox[11]->Pot( new Polarization ( new ArgonII , new MolecularHydrogen, alphaH2 ) );

// elastic + inelastic
cibox[11]->LoadInelastic();

// Ar+2 - H
cibox[14]->Pot ( new Polarization ( new ArgonII, new Hydrogen, alphaH ) ); 

// Ar+2 - H2
cibox[16]->Pot ( new Polarization ( new MolecularHydrogen, new ArgonII, alphaH2 ) ); 

// H - H+
cibox[19]->Load(false);

// getting the interface for initialization
InteractionInterface* i = cibox[19]->GetIntInterface();

// Multi state, with threshold, for elastic and charge transfer for inelastic
cibox[19]->TCScalculator = new CsHolder (

    new MultiCs ( i, {

        new AvrgChiIntegrator( i, new Morse2Param ( 51.75 , 1.677 )),
        new ThresholdCs ( i, 
            { 
                // Potential BEFORE 10. eV 
                new AvrgChiIntegrator ( i, new PowerPot ( -282. , 5.8 )) ,

                // Potential AFTER 10. eV
                new AvrgChiIntegrator ( i, new PowerPot ( -18.76 , 3.47 )) 
            },
            {
                // Threshold energy
                10.
            } 
        )
    },

    // degeneracies for the two states
    {
        2.,
        2.
    }

    ), // elastic done here.
    
    // inelastic
    new ChargeTransferCs ( i, 28.69, 1.3 ) 

);

// H-e- Integration of Momentum Transfer Cross Section 
cibox[21]->Load(false); 
cibox[21]->LoadElastic();  

// H+-H2
cibox[23]->Pot( new Polarization ( new HydrogenI , new MolecularHydrogen , alphaH2 ) );
cibox[23]->LoadInelastic();

// H2-e- Integration of differential scattering cross section
cibox[26]->Load(false); 
cibox[26]->LoadDCS(); 

// Output modules initialization
Thermodynamics thSolver;
ThermodynamicsCsv th(&thSolver, folder);
DevotoTpCsv Dev ( &cibox, folder ); 
ZhangTpCsv  zm  ( &cibox, folder ); 

// Loop on NEparam
for (size_t i = 0; i < NEparam.size(); i++) {

    // NEparam set to Mixture object
    mix->theta->set(NEparam[i]);

    // int to write filenames
    int str = NEparam[i]; 

    // Computing and printing data
    th.Print ( "Theta" + std::to_string(str), T, mix ); 
    Dev.Print ( "Dev1966_Theta" + std::to_string(str), T, mix );
    zm.Print ( "ZM2013_Theta" + std::to_string(str), T, mix );

}

return 0;

} 
