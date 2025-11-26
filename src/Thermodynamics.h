// PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
// (University of Bologna, Italy)                   // 
// Licensed under CC BY 4.0.                        // 
// To view a copy of this license, visit:           // 
// https://creativecommons.org/licenses/by/4.0/     // 

#ifndef THERMODYNAMICS_H
#define THERMODYNAMICS_H

#include <vector>
#include <stdexcept>
#include <string>

class GasMixture;

/// @brief Computes thermodynamic properties for a given gas mixture.
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

    /// @brief Vector of hentalphy values for species in the mixture [J/kg].
    std::vector<double> henthalpies;

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

};

class ThermodynamicsDHcorrected : public Thermodynamics {

    /** @brief Computes thermodynamic properties as parent class 
     * and just applies for DH corrections */
    void computeThermodynamics(GasMixture& gasmix);
    
};

#endif
