#include "GasMixture.h"
#include "CiBox.h"
#include "CollisionIntegral.h"
#include "Potential.h"

int main() {

    // GasMixture definition
    GasMixture* mix = new GasMixture ( 
        
        300.            , 101325.       , 
        new Argon       , new ArgonI    , new ArgonII           ,
        new Hydrogen    , new HydrogenI , new MolecularHydrogen ,
        new Electron 
        
    );

    // Editable Collision Integrals Box
    CiBox cibox (mix);

    // Ar-Ar editing
    cibox[0]->Pot(new Morse3Param(1.,1.,1.)); 

    // Get info 
    cibox.info();

    return 0;

} 
