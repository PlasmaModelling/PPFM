#ifndef INTERACTION_H
#define INTERACTION_H

#include"Potential.h"

/** @brief Interface class for the interaction, passed as pointer to TcsCalcolators classes
 * to gather for informations on the Interaction */
class InteractionInterface {
    
    public:

    /// @return the specie with its base type.
    virtual Species* GetSp1() = 0 ;
    
    /// @return the specie with its base type.
    virtual Species* GetSp2() = 0 ;

    /// @return std::string "Formula_Formula"
    virtual std::string InteractionName() = 0 ;
};

/** @brief Abstract class that represent the interaction through 
 * chemical species involved in a particle binary collison. It serves as an abstract 
 * concept to compute the physical quantities of transport cross sections 
 * and collision integrals. 
 * @tparam T1 1-st chemical specie with its native type.
 * @tparam T2 2-nd chemical specie with its native type.*/
template<typename T1, typename T2>
class Interaction : public InteractionInterface {

protected:

    T2* sp2 ;  
    T1* sp1 ; 
    
public:
   
    Interaction() ;
    /// @brief Construct the object by initializing the chemical species passed to it.
    /// @param t1 Pointer to generic chemical specie
    /// @param t2 Pointer to generic chemical specie
    Interaction( T1* t1,T2* t2 ) ;

    ~Interaction() ;

    
    /// @return the specie with its base type.
    Species* GetSp1() { return sp1 ; } ;
    /// @return the specie with its base type.
    Species* GetSp2() { return sp2 ; } ;

    /// @return the specie with its native type.
    T1* getSp1() ;
    /// @return the specie with its native type.
    T2* getSp2() ;

    /// @brief InteractionInterface ralization
    std::string InteractionName() ;

};

// ________________________ Implementazione Interaction ________________________

template<typename T1, typename T2>
Interaction<T1,T2>::Interaction() {
    sp1 = new T1 ;
    sp2 = new T2 ;
}

template<typename T1, typename T2>
Interaction<T1,T2>::Interaction(T1* t1, T2* t2): sp1(t1), sp2(t2) {}

template<typename T1, typename T2>
Interaction<T1,T2>::~Interaction() {}

template<typename T1, typename T2>
T1* Interaction<T1,T2>::getSp1(){return sp1;}

template<typename T1, typename T2>
T2* Interaction<T1,T2>::getSp2(){return sp2;}

template< typename T1, typename T2 >
std::string Interaction<T1,T2>::InteractionName() {
    return sp1-> getFormula() +"_"+sp2-> getFormula() ; 
}

#endif