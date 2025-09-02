 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef CIBOX_H
#define CIBOX_H

#include<vector>
#include<string>

// forward declarations
class HybridInterface ;
class Mixture ;
class GasMixture ; 

/** 
 * @class CiBox
 * @brief A container to store and operate on different Collision Integrals objects.
 * @details The class creates and stores the CollisionIntegrals Objects for a given Mixture. \n
 * This class implements useful methods to compute a series of collision integrals 
 * and store the computed values. \n
 * It also provide some operators to easy-call for methods of stored objects 
 * or to return the computed values. \n
 * @see <a href="../../articles/TransportFormulas.pdf"> TransportFormulas.pdf </a>
 * eq.17 and related reference. */
class CiBox {

    bool firstCall = true ;  

    /// @brief Collects the Collision Integrals object of the given Mixture 
    std::vector<HybridInterface*> integrals ; ///< [Ang^2]
    
    std::string themix ; 

    public :
    
    /** 
     * @brief Constructor 
     * @param mix Mixture object of chemical species, 
     * @see class Mixture for further details */
    CiBox(Mixture* mix ) ;
    
    /**
     * @brief Operator that returns the i-th stored Collision Integral object
     * @param i-th component of class member integrals  
     * @return HybridInterface* - base pointer to CollisionIntegrals objects */
    HybridInterface* operator[](int i) ;

    /**
     * @brief Operator that returns the 16 values needed for the 4th order approximation 
     * of the Chapman-Enskog method
     * @param i-th component of class member integrals  
     * @details Needed values of the collision integrals are : \n
     * (1,1),..,(1,7) \n
     * (2,2),..,(2,6) \n
     * (3,3),..,(3,5) \n
     * (4,4)     
     * @return std::vector<double> of i = 1..16 CollisionIntegrals */
    std::vector<double> operator()(int i) ;

    /**
     * @brief Prints Collision Integrals to CSV files for each interaction in CiBox
     * @param Ti Vector of temperatures
     * @param gasmix Pointer to the gas mixture
     * @param folder Output folder for CSV files
     * @details For each HybridInterface in CiBox, a CollisionIntegralCsv is created 
     * and used to print the corresponding integrals at the given temperatures. */
    void PrintCollisionIntegrals ( const std::vector<double>& Ti, GasMixture* gasmix, const std::string& folder );

    /**
     * @brief Prints TransportCrossSections to CSV files for each interaction in CiBox
     * @param Ti Vector of temperatures
     * @param gasmix Pointer to the gas mixture
     * @param folder Output folder for CSV files
     * @details For each HybridInterface in CiBox, a CollisionIntegralCsv is created 
     * and used to print the corresponding integrals at the given temperatures. */
    void PrintTransportCrossSection ( const std::vector<double>& Ti, GasMixture* gasmix, const std::string& folder );

    /**
     * @brief Compute and store for each collision integral in the collection 
     * @param temperature of the GasMixture [K°] 
     * @param lambda DebyeLength from the computed Composition [m]
     * @see Class Gas and class Composition for further details */
    void computeCollisionIntegrals(double Te, double Th, double lambda) ;
    
    /**
      * @brief Returns the number of binary interactions for the given mixture.
      * @details There are N(N+1)/2 binary interactions 
      * without repetitions in a Mixture of N chemical species. */ 
    int InteractionsNumber() ;
    void info() ; 
    void loadAll( bool b ) ; 
} ;

//_________________________ Implementation ____________________________

#endif