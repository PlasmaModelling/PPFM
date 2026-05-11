// PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
// (University of Bologna, Italy)                   // 
// Licensed under CC BY 4.0.                        // 
// To view a copy of this license, visit:           // 
// https://creativecommons.org/licenses/by/4.0/     // 

#ifndef THERMODYNAMICS_H
#define THERMODYNAMICS_H

/**
 * @file Thermodynamics.h
 * @brief This file contains the Thermodynamics class hierarchy for computing thermodynamic properties of gas mixtures.
 * The Thermodynamics class computes properties such as density, specific heats, and speed of sound based on the composition and state of a GasMixture. 
 * The ThermodynamicsDHcorrected class applies Debye–Hückel corrections to the computed properties for charged species.
*/

#include <vector>
#include <stdexcept>
#include <string>

// Forward declarations
class GasMixture;

/**
 * @class Thermodynamics
 * @brief Class for computing thermodynamic properties of a gas mixture.
 * @details Iterative dT and dP are implemented for stability of properties 
 * finite differences and Partition Functions derivates. 
*/
class Thermodynamics {

    protected:

    /** @brief Vector of computed thermodynamic quantities.
    * @details Indexing:  
    * [0] ρ [kg/m³]  
    * [1] R [J/(kg·K)]  
    * [2] hₑ [J/kg]  
    * [3] hₕ [J/kg]  
    * [4] eₑ [J/kg]  
    * [5] eₕ [J/kg]  
    * [6] Cₚ [J/(kg·K)]  
    * [7] Cᵥ [J/(kg·K)]  
    * [8] γ [–]  
    * [9] a [m/s] */
    std::vector<double> Td;

    /// @brief Vector of hentalphy values for species in the mixture [J/kg] (per particle).
    std::vector<double> henthalpies ;

    /// @brief Vector of dH/dT values for species in the mixture [J/(kg·K)] (per particle).
    std::vector<double> dhdT ; 

    /// @brief saving final dRdP to correct speed of sound in DH corrections 
    double dRdP_final ; 

    public:
    
    /// @brief Constructs the object and initializes the thermodynamic vector.
    Thermodynamics() { Td.resize(10, 0.0); }

    /**
     * @brief Computes thermodynamic properties from a given gas mixture state.
     * @param gasmix Reference to the GasMixture object (provides composition, T, P). */
    void computeThermodynamics(GasMixture& gasmix);

    // ---------------- Getters ---------------- //

    double rho()   const { return Td[0]; }
    double R()     const { return Td[1]; }
    double he()    const { return Td[2]; }
    double hh()    const { return Td[3]; }
    double ee()    const { return Td[4]; }
    double eh()    const { return Td[5]; }
    double cp()    const { return Td[6]; }
    double cv()    const { return Td[7]; }
    double gamma() const { return Td[8]; }
    double a()     const { return Td[9]; }
    
    double h(int i) const{ return henthalpies.at(i); }
    double dHdT(int i) const{ return dhdT.at(i); }
};

/**
 * @class ThermodynamicsDHcorrected
 * @brief Computes thermodynamic properties with Debye–Hückel corrections. 
*/
class ThermodynamicsDHcorrected : public Thermodynamics {

    /** @brief Computes thermodynamic properties as parent class 
     * and just applies for DH corrections */
    void computeThermodynamics(GasMixture& gasmix);
    
};

#endif
