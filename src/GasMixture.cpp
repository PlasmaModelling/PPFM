 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "GasMixture.h"

// Destructor
GasMixture::~GasMixture() {
    delete Comp;
}

// Composition solver setter
void GasMixture::setCompositionSolver(Composition* solver) {
    if (Comp) delete Comp;
    Comp = solver;
}

// Setter T/P with composition recomputation
void GasMixture::setT(double temperature) {
    T = temperature;
    Comp->CompositionSolve(this, this);
}

// Setter T/P with composition recomputation
void GasMixture::setP(double pressure) {
    P = pressure;
    Comp->CompositionSolve(this, this);
}

// Restart composition and recompute internal state
void GasMixture::restartComposition() {

    PfBox* tempPFBOX = Comp->Qbox;
    delete Comp; 
    Comp = new GodinTrepSahaSolver(this, this);
    Comp->setPfBox(tempPFBOX);
    Comp->CompositionSolve(this, this);

}