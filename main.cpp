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

int main() {

    std::cout<<"Test on pure Argon 4 species: Ar, Ar+, Ar+2, e- . " << std::endl;

    // GasMixture definition
    auto mix = new GasMixture (

        300.            , 101325.       ,
        new Argon       , new ArgonI    , new ArgonII           ,
        new Electron

    ) ;

    // non-equilibrium parameter Te/Th
    mix->theta->set(1.) ;

    // Transport init with default CiBoxes
    DevotoTP transp ( mix ) ;
    ZhangMurphyTP zmtransp ( mix ) ;

    std::cout << "LTE, T.range 300-3000K." << std::endl;

    // set Temperature range to compute on
    std::vector<double> T = arange (300., 3100., 100.) ;

    std::cout << "Required collision integrals are loaded." << std::endl;

    // loop
    for (size_t i = 0; i < T.size(); i++) {

        auto start = std::chrono::high_resolution_clock::now();

        std::cout << "T = "<<T[i]<<std::endl ;

        // Setting temperature, composition is computed.
        mix->setT(T[i]) ;

        // compute
        transp.computeTransport(mix) ;

        auto end = std::chrono::high_resolution_clock::now() ;
        std::chrono::duration<double> duration = end - start ;
        std::cout <<"Devoto  elapsed in: "<< duration.count() << "s" << std::endl ;

        start = end ;

        zmtransp.computeTransport(mix) ;

        end = std::chrono::high_resolution_clock::now() ;
        std::cout <<"Zhang   elapsed in: "<< duration.count() << "s" << std::endl ;

    }

    return 0 ;
}
