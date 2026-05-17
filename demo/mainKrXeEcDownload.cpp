#include"Species.h"
#include"ParfCalculator.h"

// NO_MIXTURE Download Electronic configuration of Xe+4 and Kr+4

int main() {

    std::vector<Species*> AllOfTheSpecies = {

        new KryptonIV,
        new XenonIV 

    } ; 

    std::vector<std::string> AllOfTheLinks { 
        "https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=Kr+V&units=1&format=3&output=0&page_size=15&multiplet_ordered=0&level_out=on&g_out=on&temp=&submit=Retrieve+Data",
        "https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=Xe+V&units=1&format=3&output=0&page_size=15&multiplet_ordered=0&level_out=on&g_out=on&temp=&submit=Retrieve+Data"
    } ; 

    for (size_t i = 0; i < AllOfTheSpecies.size(); i++)
    {
        ElectronicAtomicPF ( AllOfTheSpecies[i], AllOfTheLinks[i] ) ;
    }

    return 0 ;

} 
