 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef NON_EQUILIBRIUM_PARAMETER   
#define NON_EQUILIBRIUM_PARAMETER   

class Theta {
    
    double nonEquilibriumParameter ;

    public:
    
    Theta (){
        nonEquilibriumParameter = 1. ;
    }

    double get() {
        if ( nonEquilibriumParameter > 1. )
            return nonEquilibriumParameter ; 
        else
            return 1. ;
    }

    void set ( double value ) {
    
        nonEquilibriumParameter = value ; 
    }
    
    ~Theta() {} ;

};

#endif