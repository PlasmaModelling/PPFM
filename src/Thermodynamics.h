 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef THERMODYNAMICS_H
#define THERMODYNAMICS_H

#include<vector>
#include<stdexcept>
#include<string>
#include"DataPrinter.h"

class GasMixture ;

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

    public:
    
    /// @brief Constructs the object and initializes the thermodynamic vector.
    Thermodynamics() { Td.resize(10, 0.0); }

    /**
     * @brief Computes thermodynamic properties from a given gas mixture state.
     * @param gasmix Reference to the GasMixture object (provides composition, T, P).
     * @details 
     * - Uses finite differences to estimate derivatives via shifted thermodynamic states.
     * - Partition functions are recomputed at T ± ΔT and P ± ΔP.
     * - Uses the Godin formulation (eqs. 41–43) to compute Cₚ, Cᵥ, γ, and a.
     * - Restores the mixture to its original state after computation.
     * @see class PfBox for partition function evaluation.
     * @see class GasMixture for thermodynamic context.  */
    void computeThermodynamics(GasMixture& gasmix);
};


//______________________________ Implementation ______________________________

/// @brief Outputs thermodynamic properties to CSV.
/// @see class Thermodynamics for property computation.
/// @see class DataPrinter for CSV interface.
class ThermodynamicsCsv : public Thermodynamics, public DataPrinter {

    /// @brief Builds the output filename with "TH_" prefix.
    std::string BuildFileName(const std::string& filename) const;

    /// @brief Prepares the CSV header row.
    void PrepareHeader();

    /**
     * @brief Computes thermodynamic properties over a temperature range.
     * @param x Vector of electron temperatures (Te/θ).
     * @param gasmix Pointer to the gas mixture.
     * @details For each temperature, the following properties are computed and stored:
     * ρ [kg/m³], Cp [J/kg·K], hₑ + hₕ [J/kg], γ, a [m/s]. Initial gas state is restored at the end.
     * @see Thermodynamics::computeThermodynamics */
    void PrepareData(const std::vector<double>& x, GasMixture* gasmix) override;

    /**
     * @brief Prints a message confirming successful output.
     * @param filename Name of the generated file. */
    void PrintMessage(const std::string& filename) override;

    public:
    
    /// @brief Default constructor.
    ThermodynamicsCsv();

    /// @brief Constructor with output folder specification.
    ThermodynamicsCsv(const std::string& folder);

};

#endif