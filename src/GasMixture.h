#ifndef GASMIXTURE_H
#define GASMIXTURE_H

#include <stdexcept>

#include "Composition.h"
#include "Mixture.h"

/// @brief Models a gas mixture combining thermodynamic state and species composition.
/// @see class Gas for pressure and temperature state.
/// @see class Mixture for chemical species management.
class GasMixture : public Gas, public Mixture {

    public:

    /// @brief Pointer to the composition solver.
    Composition* Comp;

    /// @brief Constructs a GasMixture from a set of species.
    template<typename... SpeciesTypes>
    GasMixture(SpeciesTypes... newSpecies);

    /// @brief Constructs a GasMixture with specified temperature, pressure, and species.
    template<typename... SpeciesTypes>
    GasMixture(double Temp, double Press, SpeciesTypes... newSpecies);

    /// @brief Sets the temperature and recomputes the composition.
    void setT(double temperature) override;

    /// @brief Sets the pressure and recomputes the composition.
    void setP(double pressure) override;

    /**
     * @brief Sets the mole fractions for the given set of species.
     * @param doublevalues List of mole fractions.
     * @param species Pointers to species to which mole fractions refer.
     * @details If only M–2 species are provided, the last one is auto-completed. Species are
     * ordered lexicographically and checked for consistency with the mixture. The function 
     * ensures normalization and raises exceptions if mismatches occur.
     * @throws std::invalid_argument if input sizes mismatch or species are missing or 
     * mole fractions do not sum to 1. */
    template<typename... speciestype>
    void setMoleFractions(std::initializer_list<double> doublevalues, speciestype... species);

    /// @brief Resets the composition object and recomputes the internal state.
    void restartComposition() override;
};


//_________________________ Implementazione _________________________


template<typename... SpeciesTypes>
GasMixture::GasMixture(SpeciesTypes... newSpecies): Mixture(newSpecies...){
    Comp = new Composition ( this , this ) ;                                    
    /* Comp->compositionSolve ( this , this ) ; */
}

template<typename... SpeciesTypes>
GasMixture::GasMixture(double Temp, double Press ,SpeciesTypes... newSpecies): Mixture(newSpecies...), Gas(Press,Temp ){
    Comp = new Composition ( this , this ) ;                                    
    /* Comp->compositionSolve ( this , this ) ; */
}

template < typename... speciestypes >
void GasMixture::setMoleFractions ( std::initializer_list<double> values , speciestypes... speciespointers ) {
    
    std::vector < Species* > givenSpecies ; 
    std::vector < double   > givenValues  ;

    ( givenSpecies.push_back ( speciespointers ) , ... ) ;
    for ( auto value : values )
        givenValues.push_back ( value ) ;

    if ( givenSpecies.size() == M-2 ) {
        double sum = 0.0 ; 
        for (int i = 0 ; i < givenSpecies.size() ; i++)
            sum += givenValues[i] ; 
        
        std::vector<Species*> existent = std::get<0>(molefractions) ;
        for (int i = 0; i < existent.size(); i++) {
            bool found = false ; 
            for (int j = 0; j < givenSpecies.size(); j++) {
                found = existent[i]->getFormula() == givenSpecies[j]->getFormula() ;
                if ( auto* ptr =  dynamic_cast < PolyAtomicMolecule* > (givenSpecies[j]) ) {
                    for (int k = 0; k < ptr->numberOfCostituents(); k++) {
                        Element* elem = (*ptr)[k] ;
                        found = existent[i]->getFormula() == elem->getFormula() ;
                    }
                }
                if (found) {
                   break;
                }       
            }   
            if ( !found ) {
                givenSpecies.push_back(existent[i]) ;
                givenValues.push_back(1.-sum) ;
                break;
            }
        }
    }
    

    // check input sizes
    if ( givenSpecies.size() != givenValues.size() ) {
        throw std::invalid_argument(
            "Different species numbers or molefractions, they have to be equal, check the input"
        ) ;
    }
    
    // orders the same way of chemical species. 
    if (givenSpecies.size() >= 2) {   

        Species* tmp ;
        double tmpdbl;
        for (int i = 1; i < givenSpecies.size() ; i++) {

            if ( givenSpecies[i]->getFormula() < givenSpecies[i-1]->getFormula() ) {

                tmp = givenSpecies[i-1] ; 
                tmpdbl = givenValues[i-1];
                givenSpecies[i-1] = givenSpecies[i] ;
                givenValues[i-1] = givenValues[i] ;
                givenSpecies[i] = tmp ;
                givenValues[i] = tmpdbl ;
                i = 0 ;
            
            }
        }
    
        // check the givenSpecies are present within the mixture
        for ( int j = 0 ; j < givenSpecies.size() ; j++ ) {

            bool found = false ; 
            
            for (int i = 0 ; i < this->getN() ; i++) {
                
                if (givenSpecies[j]->getFormula() == (*this)(i)->getFormula() ) {
                    
                    found = true ; 
                    break ; 
                
                }        
            }
            if (!found) {

                throw std::invalid_argument(
                    "Error: specie lacking in the mixture, check input and Mixture construction"
                ) ;
            }
            
            // chek inserted values
            double total = 0.0 ;
            for (int i = 0 ; i < givenValues.size() ; i++ )
                total += givenValues[i] ; 
            double rounded = std::round(total) ; 
            if ( std::abs ( total - rounded ) > 1e-5 )
                throw std::invalid_argument(
                    "The sum of molefractions far exceeds 1.0 \n "
                    "\t check molefractions."
            ) ; 
        }    
    }

    // set mole fractions
    std::get<0>(molefractions).clear() ;
    std::get<1>(molefractions).clear() ;
    for ( int i = 0 ; i < givenSpecies.size() ; i++ ) {
        std::get<0>(molefractions).push_back(givenSpecies[i]) ;
        std::get<1>(molefractions).push_back(givenValues[i] ) ;
    }

    // compute composition with newer mole fractions
    restartComposition() ;
    Comp->compositionSolve(this,this) ;
}


#endif
