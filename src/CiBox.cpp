 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "CiBox.h"
#include "CollisionIntegral.h"
#include "GasMixture.h"

CiBox::CiBox(Mixture* mix ) {
    
    for (int i = 0; i < mix->cispecies.size() ; i++) {
        std::string formula ;
        for (int j = i; j < mix->cispecies.size() ; j++) {

            /* Creates rigth CollisionIntegrals with native pointers from 
            variant mix->cispecies in order to 
            built the correct CollisionIntegral<T1,T2> templetized object */
            std::visit([this](auto&& sp1, auto&& sp2){
                HybridInterface* integral = new CollisionIntegral(sp1,sp2) ;
                integrals.push_back(integral) ;
            }, mix->cispecies[i], mix->cispecies[j]);
        }
    }
    
    themix = "Binary interactions for a mixture of :\n " ;
    for (int i = 0; i < mix->cispecies.size()-1 ; i++) 
        themix += (*mix)(i)->getFormula() +", "; 
    
    int i = mix->cispecies.size() - 1 ; 
    themix += (*mix)(i)->getFormula() ;

}

HybridInterface* CiBox::operator[](int i){
    if (i < 0 || i >= integrals.size()) 
        throw std::out_of_range("Index out of CiBox::integrals range") ;
    else
        return integrals[i] ;
}

std::vector<double> CiBox::operator()(int i){
    if (i < 0 || i >= integrals.size()) 
        throw std::out_of_range("Index out of CiBox::integrals range") ;
    else
        return integrals[i]->omega4th ;
}

void CiBox::computeCollisionIntegrals(double Te, double Th, double lambda) {
    
    bool exception_occurred = false;
    std::vector<std::string> error_messages(integrals.size());

    if (firstCall) {
        // Esecuzione SEQUENZIALE alla prima chiamata
        for (int i = 0; i < integrals.size(); i++) {
            try {
               
                integrals[i]->ComputeCollisionIntegral(Te, Th, lambda);
               
            } catch (const std::exception& e) {
                exception_occurred = true;
                error_messages[i] = "Error in the " + std::to_string(i) + "-th collision integral: "
                                    + integrals[i]->InteractionName() + " " + e.what() + "\n\n\t"
                                    + "Initialize in main:\n\t  ciboxname[" + std::to_string(i) 
                                    + "]->setPot( see Potential.h for more );\n" 
                                    + "\t otherwise, \n" 
                                    + "\t - try loading it from Raw files placed in data/CollisionIntegral \n"
                                    + "\t - prefer a different calculation method, \n"
                                    + "\t see CollisionIntegral.h for more. ";
            }
        }
        firstCall = false;  // Imposta flag per chiamate successive
    } else {
        // Esecuzione PARALLELA dalle chiamate successive
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < integrals.size(); i++) {
            try {
                integrals[i]->ComputeCollisionIntegral(Te, Th, lambda);
            } catch (const std::exception& e) {
                #pragma omp critical
                {
                    exception_occurred = true;
                    error_messages[i] = "Error in the " + std::to_string(i) + "-th collision integral: "
                                        + integrals[i]->InteractionName() + " " + e.what() + "\n\n\t"
                                        + "Initialize in main:\n\t  ciboxname[" + std::to_string(i) 
                                        + "]->setPot( see Potential.h for more );\n" 
                                        + "\t otherwise, \n" 
                                        + "\t - try loading it from Raw files placed in data/CollisionIntegral \n"
                                        + "\t - prefer a different calculation method, \n"
                                        + "\t see CollisionIntegral.h for more. ";
                }
            }
        }
    }

    // Stampa errori
    if (exception_occurred) {
        for (const auto& msg : error_messages) {
            if (!msg.empty()) {
                std::cerr << msg << std::endl;
            }
        }
        std::exit(EXIT_FAILURE);
    }
}


int CiBox::InteractionsNumber() {
    return integrals.size() ; 
}

void CiBox::info() {

    std::cout << std::left ; 
    std::cout << themix << std::endl;
    std::cout << "[i] |    Interaction    |             Calculator type             " << std::endl ; 
    std::cout << "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾" << std::endl ;
    for(int i = 0; i< integrals.size() ; i++ ) {
        std::cout << std::left ; 
        std::cout << std::setw(8) << i; integrals[i]->info() ; std::cout << std::endl ; 
    }
    std::cout << std::endl ;

}

void CiBox::loadAll( bool b ) {
    for (int i = 0; i < integrals.size(); i++) {
        integrals[i]->Load(b) ;
    }
}

void CiBox::PrintCollisionIntegrals ( const std::vector<double>& Ti, GasMixture* gasmix, const std::string& folder ) {

    for ( auto* ci : integrals ) {

        CollisionIntegralCsv writer ( ci ) ; 

        writer.customFolder = folder;

        writer.Print ( ci->InteractionName(), Ti, gasmix );

    }
}

void CiBox::PrintTransportCrossSection ( const std::vector<double>& Ti, GasMixture* gasmix, const std::string& folder ) {

    for ( TcsInterface* tcs : integrals ) {

        TransportCrossSectionCsv writer(tcs) ; 

        writer.customFolder = folder;

        writer.Print ( tcs->GetIntInterface()->InteractionName(), Ti, gasmix );

    }
}