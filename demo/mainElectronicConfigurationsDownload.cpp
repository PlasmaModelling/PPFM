#include"Species.h"
#include"ParfCalculator.h"

/* 

    TEST MAIN FOR ELECTRONIC CONFIGURATION DATA DOWNLOAD FROM NIST ATOMIC SPECTRA DATABASE

    #Include JUST the modules you intend to use.
    Check src/header.h files for info on the implemented classes

*/

// NO_MIXTURE used, Download Electronic configuration of Xe+4 and Kr+4

int main() {

    std::vector<Species*> AllOfTheSpecies = {

        new KryptonIV,
        new XenonIV 

    } ; 

    // LINKS TO ONLINE RESOURCES IN THE NIST ATOMIC SPECTRA DATABASE.
    // IN THE ATOMIC SPECTRA DATABASE SELECT THE SPECTRUM OF THE DESIRED SPECIES, TAB-DELIMITED,
    // CHECK "LEVEL" AND "g" ONLY, UNCHECK EVERYTHING ELSE, AND SELECT "eV" AS ENERGY UNIT.
    std::vector<std::string> AllOfTheLinks { 
    
        "https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=Kr+V&units=1&format=3&output=0&page_size=15&multiplet_ordered=0&level_out=on&g_out=on&temp=&submit=Retrieve+Data",
        "https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=Xe+V&units=1&format=3&output=0&page_size=15&multiplet_ordered=0&level_out=on&g_out=on&temp=&submit=Retrieve+Data"
    
    } ; 

    // Constructing ElectronicAtomicPF objects, the overloaded constructor will handle the download and the csv file writing.
    for (size_t i = 0; i < AllOfTheSpecies.size(); i++)
        ElectronicAtomicPF ( AllOfTheSpecies[i], AllOfTheLinks[i] ) ;
    
    std::cout << "Electronic configurations downloaded from online resources for " << AllOfTheSpecies.size() << " species: " << std::endl ;
    for (size_t i = 0; i < AllOfTheSpecies.size(); i++)
        std::cout << " - " << AllOfTheSpecies[i]->getFormula() << " Electronic Configuration downloaded from\n\t"<< AllOfTheLinks[i] << "\n and saved in data/ElectronicConfigurations/" << AllOfTheSpecies[i]->getFormula() << "_ElConfig.csv" << std::endl ;
    
    return 0 ;

} 
