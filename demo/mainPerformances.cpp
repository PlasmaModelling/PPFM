/* 

    TEST MAIN FOR TRANSPORT PROPERTIES RAW CALCULATION PERFORMANCES

    #Include JUST the modules you intend to use.
    Check src/header.h files for info on the implemented classes
 
    Test on pure Argon 4 species: Ar, Ar+, Ar+2, e- 
    LTE, T.range 300-3000K.
    Required collision integrals are loaded. 

*/

#include "GasMixture.h"
#include "Devoto.h"
#include "ZhangMurphyTP.h"
#include <chrono>
#include <iostream>

int main() {
    
    std::cout << "Test on pure Argon 4 species: Ar, Ar+, Ar+2, e-." << std::endl;

    // GasMixture definition
    auto mix = new GasMixture ( 
        
        300.            , 101325.       , 
        new Argon       , new ArgonI    , new ArgonII           , 
        new Electron 
        
    ) ;

    // Non-equilibrium parameter Te/Th
    mix->theta->set(1.) ; 
    
    // Transport init with default CiBoxes
    DevotoTP transp ( mix ) ;
    ZhangMurphyTP zmtransp ( mix ) ;
        
    std::cout << "LTE, T.range 300–3000 K." << std::endl;

    // Set Temperature range to compute on
    std::vector<double> T = arange (300., 3100., 100.) ;

    std::cout << "Required collision integrals are loaded." << std::endl;

    // Loop 
    for (size_t i = 0; i < T.size(); i++) {
        
        std::cout << "T = " << T[i] << std::endl; 

        // Setting temperature, composition is recomputed.
        mix->setT(T[i]) ;
        
        // Devoto computation
        auto start = std::chrono::high_resolution_clock::now();
        transp.computeTransport(mix) ;        
        auto end = std::chrono::high_resolution_clock::now() ; 
        std::chrono::duration<double> duration = end - start ; 
        std::cout << "Devoto elapsed in: " << duration.count() << "s" << std::endl ;

        // Zhang–Murphy computation
        start = std::chrono::high_resolution_clock::now();
        zmtransp.computeTransport(mix) ;
        end = std::chrono::high_resolution_clock::now(); 
        duration = end - start;
        std::cout << "Zhang  elapsed in: " << duration.count() << "s" << std::endl ;
    }

    return 0 ;
}
