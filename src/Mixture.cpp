 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "Mixture.h"

int Mixture::getN() { return N ; }; 

int Mixture::getM() { return M ; }; 

Species* Mixture::operator()(int i) {
    
    if (i < 0 || i >= cispecies.size()) {
        throw std::out_of_range("Index out of range");
    }
    return std::visit([](auto&& arg) -> Species* { return arg; }, cispecies[i]);
}

std::vector<double> Mixture::masses() {
    
    std::vector<double> ms(N) ;
    for (int i = 0; i < N; i++) {
        ms[i] = (*this)(i)->getMass() ;
    }
    return ms ; 
}

std::vector<double> Mixture::masses(double conversion) {
    
    std::vector<double> ms(N) ;
    for (int i = 0; i < N; i++) {
        ms[i] = (*this)(i)->getMass() * conversion ;
    }
    return ms ; 
}

void Mixture::order(){
    
    AcceptedSpecies tmp ;
    
    for (int i = 1; i < N; i++) {
        if ( (*this)(i)->getFormula() < (*this)(i-1)->getFormula() ) {
    
            tmp = cispecies[i-1] ; 
            cispecies[i-1] = cispecies[i] ;
            cispecies[i] = tmp ;
            i = 0 ;
    
        }
    }
    /* 
    for (int i = 0; i < N; i++)
    {
        if ( (*this)(i)->getFormula() == "e-" ) {
            tmp = cispecies[i];
            // Sposta tutte le specie di un posto indietro
            for (int j = i; j > 0; j--) {
                cispecies[j] = cispecies[j-1];
            }
            // Metti l'elettrone in testa
            cispecies[0] = tmp;
            break;
        }
        
    }
     */

}

Mixture::~Mixture(){}

Gas::Gas(double pressure, double temperature ) :
    P(pressure), T(temperature) {
        theta = new Theta() ;
    }

double Gas::getPressure() { 
    return P ; 
}

double Gas::getTemperature() { 
    return T ; 
}

Gas::~Gas(){}

void Gas::setT ( double temperature ) {
    T = temperature;
}

void Gas::setP( double pressure    ) {
    P = pressure ;
}

bool Mixture::isBase(Species* specie) {
    for (auto& baseSpecie : std::get<0>(molefractions)) {
        if (typeid(*specie) == typeid(*baseSpecie)) {
            return true;
        }
    }
    return false;
}