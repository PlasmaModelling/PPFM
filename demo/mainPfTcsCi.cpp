/* 

    TEST MAIN FOR TRANSPORT CROSS SECTIONS AND COLLISION INTEGRALS CALCULATIONS.

    #Include JUST the modules you intend to use.
    Check src/header.h files for info on the implemented classes
 
    this demo is intended to demonstrate the potentiality within PPFM 
    to compute large volumes of data required for thermodynamic and transport properties calculations 
    
*/

#include "GasMixture.h"
#include "PfBox.h"
#include "CiBox.h"
#include "PartitionFunction.h"
#include "CollisionIntegral.h"

int main() {
    
    std::cout<<" Testing Partition Functions, Transport Cross Sections and Collision Integrals calculations " << std::endl;

    /* It is convenient to define a gas mixture object to initialize boxes that will store  
    ** ordered obects for partition functions and collision integrals combinations.
    ** but single calculations can be performed as well, 
    ** see PertitionFunction.h and CollisionIntegral.h for further details */
    auto mix = new GasMixture ( 
        
        300.         , 101325.     , 
        new Argon    , new ArgonI  , new ArgonII  , 
        new Carbon   , new CarbonI , new CarbonII , 
        new Oxygen   , new OxygenI , new OxygenII , 
        new Electron
    
    ) ;
    
    // Temperature range definition
    std::vector<double> T = arange(300., 45100., 100.) ;
    
    // Boxes initialization.
    PfBox pfbox ( mix ) ;
    CiBox ciBox ( mix ) ;
    
    std::cout<< mix->getN() << " species in the mixture, hence, " << mix->getN() << " Partition Functions." << std::endl;
    std::cout<< mix->getN() << " species in the mixture, hence, " << ciBox.InteractionsNumber() << " (N(N+1)/2)) binary interactions." << std::endl;
    std::cout<< "Only Non-Coulomb interactions have to be computed." << std::endl;

    std::cout << " Boxes info pre editing\n" << std::endl;
    pfbox.info() ;
    ciBox.info() ;

    // Setting Partition functions calculation from electronic configuration
    pfbox.AllAbInitio(); 
    
    std::cout << " PfBox info post editing\n" << std::endl;
    pfbox.info() ;

    // Output folder for CSV files
    std::string folder = "demo/Pf_And_CI/" ; 

    std::cout << " Partition Functions calculation\n" << std::endl;
    
    // Print Partition Functions to CSV files
    pfbox.PrintPartitionFunctions(T, mix,folder+"PartitionFunctions/") ;

    std::cout<<std::endl;

    /* This test Transport Cross Sections and Collision Integrals calculations will be performed 
    ** with the phenomenological potential developed by Pirani et.al. 
    ** parameters for the potential have been computed apart, following */

    /* As a good approximation let's consider polarization potential for electron scattering,
    ** Polarization potential depend on the polarizability of the neutral particle involved in the interaction, 
    ** Polarizabilities from DOI: 10.1140/epjd/e2009-00192-7 
    ** Laricchiuta et.al. "High Temperature Mars Atmosphere Part I: Transport cross sections" 
    ** by analizing cibox.info() we can see Ar-e-,C-e- and O-e- interactions are in positions
    ** 9, 33 and 48 */
    std::vector<size_t> PolarizationIndices = {9, 33, 48} ;
    std::vector<double> Polarizabilities = {
        1.640, // Ar 
        // 0.919, // Ar+
        // 0.391, // Ar+2
        1.76,  // C
        // 0.79,  // C+
        // 0.368, // C+2
        0.80 //,  // O
        // 0.279, // O+
        // 0.228  // O+2
    };

    // Counter to navigate through polarization interactions parameters
    size_t polCounter = 0 ; 

    /* Parameters for the Capitelli potential, computed, for the present work following the procecedure explained in 
    ** DOI: 10.1016/j.cplett.2007.07.097 
    ** A. Laricchiuta, G. Colonna, D. Bruno, R. Celiberto, C. Gorse, F. Pirani, M. Capitelli 
    ** "Classical transport collision integrals for a Lennard-Jones like phenomenological model potential" */
    std::vector<std::vector<double>> CapitelliParameters = {

        // Beta,        eps0,         Re  
        {8.12E+00,	1.13E-02,	3.79E+00 },  // Ar_Ar
        {7.60E+00,	2.42E-02,	4.67E+00 },  // Ar_Ar+
        {7.48E+00,	7.75E-02,	4.63E+00 },  // Ar_Ar+2
        {7.04E+00,	1.38E-04,	7.67E+00 },  // Ar_C
        {7.65E+00,	2.60E-02,	4.56E+00 },  // Ar_C+
        {8.64E+00,	7.70E-01,	2.61E+00 },  // Ar_C+2
        {7.26E+00,	1.74E-04,	6.83E+00 },  // Ar_O
        {7.32E+00,	7.87E-03,	5.94E+00 },  // Ar_O+
        {7.66E+00,	1.14E-01,	4.19E+00 },  // Ar_O+2
        {6.90E+00,	2.65E-03,	8.25E+00 },  // Ar+_C
        {7.06E+00,	1.71E-03,	7.61E+00 },  // Ar+_O
        {6.86E+00,	9.73E-03,	7.92E+00 },  // Ar+2_C
        {7.00E+00,	5.79E-03,	7.42E+00 },  // Ar+2_O
        {6.69E+00,	1.08E-05,	1.15E+01 },  // C_C
        {6.91E+00,	2.69E-03,	8.17E+00 },  // C_C+
        {7.15E+00,	3.11E-02,	5.92E+00 },  // C_C+2
        {6.78E+00,	8.40E-06,	1.10E+01 },  // C_O
        {6.80E+00,	1.19E-03,	9.68E+00 },  // C_O+
        {6.92E+00,	1.18E-02,	7.53E+00 },  // C_O+2
        {7.08E+00,	1.75E-03,	7.51E+00 },  // C+_O
        {7.43E+00,	2.34E-02,	5.23E+00 },  // C+2_O
        {6.90E+00,	7.91E-06,	1.03E+01 },  // O_O
        {6.93E+00,	6.89E-04,	9.18E+00 },  // O_O+
        {7.08E+00,	7.19E-03,	7.00E+00 }   // O_O+2

    };

    /* NOTE: This is just a test case with in situ definition of parameters into variables, 
    ** high automatization, ideally with arbitrarly big N number of chemical species, 
    ** can be reached in this fashion; consider, for example, the definition of maps 
    ** (either csv text files and others) with the parameters of the interactions, 
    ** and the numbering of the interactions themselves into the cibox class, 
    ** that can be automatically loaded into the program and navigated through looping and indices. */

    // Counter to navigate through Capitelli interactions parameters
    size_t capCounter = 0 ;

    // Another thing that come in handy in this case is the Interaction Interface for the i-th interaction.
    auto interaction = ciBox[0]->GetIntInterface();

    // Now we can loop through the interactions in CiBox and set the potentials.
    for (size_t i = 0; i < ciBox.InteractionsNumber(); i++) {

        // We catch the species cause we want to set potentials only for 
        // neutral-neutral, or neutral-ions interactions.
        interaction = ciBox[i]->GetIntInterface();
        auto specie1 = interaction->GetSp1() ;
        auto specie2 = interaction->GetSp2() ;

        if (specie1->getCharge()*specie2->getCharge() == 0) {

            if (i == PolarizationIndices[polCounter]) {

                // Set of the polarization potential
                ciBox[i]->Pot(new Polarization (specie1, specie2, Polarizabilities[polCounter])) ;
                polCounter++ ;
            
            } else {
            
                // Set of the Capitelli potential
                ciBox[i]->Pot(new Capitelli (specie1, specie2, CapitelliParameters[capCounter][0], 
                    CapitelliParameters[capCounter][1], CapitelliParameters[capCounter][2])) ;
                capCounter++ ;

            }
        
        } else {
            continue ;
        }
    }

    /* If we assigned every potential correctly, no more "Raw Loader" should appear in ciBox.info() */
    std::cout << " CiBox info post editing\n" << std::endl;
    ciBox.info() ;

    std::cout << " Collision Integrals calculation\n" << std::endl;
    /* TCS and CI do not depend from the state of the mixture like Partition Functions
    ** do, so, the second argument could be left empty */
    ciBox.PrintDeflectionAngles(T, nullptr, folder+"DeflectionAngles/") ;
    ciBox.PrintTransportCrossSections(T, nullptr, folder+"TransportCrossSections/") ;
    ciBox.PrintCollisionIntegrals(T, nullptr, folder+"CollisionIntegrals/") ;

    return 0 ;
    
} 