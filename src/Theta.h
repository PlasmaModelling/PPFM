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