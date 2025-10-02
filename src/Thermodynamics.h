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
#include "DataPrinter.h"

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
};

/** @brief Outputs thermodynamic properties to CSV.
 ** @see class Thermodynamics for property computation.
 ** @see class DataPrinter for CSV interface. */
class ThermodynamicsCsv : public DataPrinter {

    /// @brief Solver reference (non-owning, external lifetime managed).
    Thermodynamics* solver;

    /// @brief Build the output filename.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepare the CSV header row.
    void PrepareHeader() override;

    /// @brief Prepare the data matrix for printing.
    void PrepareData(const std::vector<double>& x, GasMixture* gasmix) override;

    /// @brief Print completion message.
    void PrintMessage(const std::string& filename) override;

public:
    
    /// @brief Constructor with explicit solver
    ThermodynamicsCsv(Thermodynamics* solver);

    /// @brief Constructor with explicit solver and custom folder
    ThermodynamicsCsv(Thermodynamics* solver, const std::string& folder);

};

#endif
