 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "GasMixture.h"

void GasMixture::setT (double temperature ) {
    T = temperature;
    Comp->compositionSolve ( this , this ) ;
};

void GasMixture::setP ( double pressure    ) {
    P = pressure ;
    Comp->compositionSolve ( this , this ) ;
};

void GasMixture::restartComposition() {
    
    PfBox* tempPFBOX = this->Comp->Qbox;
    this->Comp = new Composition (this,this) ;
    this->Comp->setPfBox(tempPFBOX) ;
    Comp->compositionSolve (this,this) ;

}
