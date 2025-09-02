 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef COMPOSITION_H
#define COMPOSITION_H

#include "PfBox.h"
#include "DataPrinter.h"

#include <map>
#include <string>
#include <typeindex>

class Gas ;
class Species ; 

/**
 * @brief Implement a matrix based algorithm to compute equilibrium composition of a Gas Mixture. 
 * @see <a href="../../articles/A Robust and Efficient Method for the Computation of Equilibrium 
 * Composition in Gaseous Mixtures.pdf"> Composition algorithm </a> for the implemented algorithm
 * details. */
class Composition {

    friend class GasMixture ;
    friend class Thermodynamics ; 

    /// @brief Composition Matrix, fixed on given Mixture
    std::vector<std::vector<double>> C ; 

    protected: 
    
    /// @brief Store composition values in particle densities [#/m^3]
    std::vector<double> ni ;

    /** @brief Compute and store in member ni the particle densities of the Species in the 
     * Mixture, at state specified by Gas.
     * @param mix Object of class Mixture.
     * @param gas Object of class Gas. */
    void compositionSolve ( Mixture* mix , Gas* gas ) ;

    public:

    /// @brief Partition Functions container
    PfBox* Qbox;
    
    /// @brief Constructor, initialize starting composition 
    /// @param mix Mixture object
    /// @param gas Gas object
    Composition( Mixture* mix , Gas* gas) ;

    /// @brief Constructor, initialize starting composition 
    /// @param mix Mixture object
    /// @param gas Gas object
    Composition( Mixture* mix , Gas* gas, PfBox* qbox ) ;

    /**
     * @brief Get the Debye Length of the computed composition 
     * @param temperature [K°]
     * @return lambdaD(T) [m] */
    double getDebyeLength ( double temperature ) ;
    
    /**
     * @brief get the i-th composition value
     * @param i 
     * @return ni[i] [#/m^3 */
    double operator() ( int i ) ;
    
    /**
     * @brief get the computed composition values as a vector of double values 
     * @return std::vector<double> ni ; */
    std::vector<double> compositions() ;
    
    /**
     * @brief get the computed composition values as a vector of double values 
     * multiplied by a conversion factor. 
     * @return std::vector<double> ni ; */
    std::vector<double> compositions(double conversion) ;

    /**
     * @brief get the computed total composition value 
     * @return sum ni ; */
    double ntot () { double sum = 0.0; for (const auto& ns : ni) sum += ns; return sum; }

    /// @brief set starting composition to member ni
    /// @param n0 starting composition
    void setn0 ( std::vector<double> n0 ) { ni = n0 ; }

    /// @brief set for a desired PF container with desired methods
    /// @param specifiedQbox 
    /// @details use PfBox[i].setMethod() to set method on the i-th Partition Function
    /// @see class PartitionFunction for the available methods 
    void setPfBox( PfBox* specifiedQbox ) { Qbox = specifiedQbox ; }    

    private:
    
    /**
     * @brief Build the row of the Composition Matrix C inspecting the incoming species 
     * type and its characteristics.
     * @param specie Object of class Species 
     * @param colmap Map of elements and related colummn in the Composition matrix C 
     * @return std::vector<double> row of composition Matrix */
    std::vector<double> Crow ( Species* specie, const std::map<std::type_index , int>& colmap )  ;
    
    /**
     * @brief Build the Composition Matrix C via the use of the method Crow.
     * @param mixx Object of class Mixture
     * @return std::vector<std::vector<double>> C */
    std::vector<std::vector<double>> CompositionMatrix ( Mixture& mixx ) ;
    
    /** 
     * @brief Build the Conservation Matrix A 
     * @details this matrix implements the conservation laws of nuclei, i.e. the stechiometry,
     * electrical neutrality of the ionized Gas Mixture and conservation of particle density.
     * @param mixx Object of class Mixture.
     * @param gass Object of class Gas.
     * @param C Composition Matrix
     * @return std::vector<std::vector<double>> A */
    std::vector<std::vector<double>> ConservationMatrix ( Mixture& mixx, Gas& gass, const std::vector<std::vector<double>>& C ) ;
    
    /// select a chemical basis of the algorithm based on 
    void baseCalc ( std::vector<int>& b, std::vector<int>& bs, const std::vector<std::vector<double>>& C ) ;
} ; 

//______________________________ Implementation ______________________________

class CompositionCsv : public Composition, public DataPrinter {

    /** @brief A pointer to use the Mixture during Printing */    
    Mixture* mixptr ; 

    /** @brief Inherited by DataPrinter, returns a default filename */
    std::string BuildFileName(const std::string& filename ) const override ;

    /** @brief Inherited by DataPrinter, returns a default string header of printables */
    void PrepareHeader() override ;

    /** @brief Inherited by DataPrinter, prepare DataPrintable::data matrix to be printed to Csv */
    void PrepareData ( const std::vector<double>& x, GasMixture* gasmix ) override ; 

    /** @brief Inherited by DataPrinter, print a default message to terminal */
    void PrintMessage(const std::string& filename) override ;
    
    public:

    /** @brief Base constructor of the class, PfBox of class Composition will be the default one */
    CompositionCsv ( Mixture* mix , Gas* gas ) ;

    /** @brief Constructor of the class, custom PfBox can be passed by the user  */
    CompositionCsv ( Mixture* mix , Gas* gas , PfBox* qbox ) ;

    /** @brief Base constructor of the class, PfBox of class Composition will be the default one,
     * a specified folder name can be passed by the user 
     * it will be placed in the /out folder of PPFM */    
    CompositionCsv ( Mixture* mix , Gas* gas , const std::string& folder ) ;

    /** @brief Constructor of the class, custom PfBox can be passed by the user, a specified folder name can be passed by the user 
     * it will be placed in the /out folder of PPFM */    
    CompositionCsv ( Mixture* mix , Gas* gas , PfBox* qbox, const std::string& folder ) ;

};


#endif 