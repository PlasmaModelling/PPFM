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
#include <vector>
#include <cmath>

class Gas;
class Species;
class Mixture;

/**
 * @brief Abstract interface for chemical composition solvers.
 *
 * A Composition object is responsible for computing the equilibrium particle densities
 * of a Gas Mixture at a given state. Different algorithms (matrix-based, Gibbs minimization, etc.)
 * can be implemented as derived classes of this abstract interface. */
class Composition {

    friend class GasMixture;
    friend class Thermodynamics;
    friend class CompositionCsv;

protected:

    /// @brief Mixture reference (defines species and elements count)
    Mixture* mixptr;

    /// @brief Gas reference (defines state variables)
    Gas* gasptr;

    /// @brief Store composition values in particle densities [#/m^3]
    std::vector<double> ni;

    /// @brief Partition Functions container (always available to solvers)
    PfBox* Qbox;

    enum class DebyeModel {

        Rat2002Th,    ///< Rat (2002) with Th (default)
        Rat2002Te,    ///< Rat (2002) with Te
        Ghourui       ///< Ghourui formulation
        // add here future extensions, implement in setDebyeModel 
        // and related private function.

    };

    /// @brief Selected Debye length formulation, defaults to Rat2002Th
    DebyeModel debyeChoice = DebyeModel::Rat2002Th;

public:

    /// @brief Base constructor
    Composition(Mixture* mix, Gas* gas) ;

    /// @brief Constructor with qbox assignment
    Composition( Mixture* mix , Gas* gas, PfBox* qbox ) ;

    virtual ~Composition() = default;

    /** 
     * @brief Compute and store in member ni the particle densities of the 
     * Species in the Mixture, at state specified by Gas. */
    virtual void CompositionSolve(Mixture* mix, Gas* gas) = 0;

    /**
     * @brief Get the Debye Length of the computed composition.
     * @param temperature Electron or heavy-particle temperature [K].
     * @return Debye length [m]. */
    virtual double getDebyeLength(double temperature) ;

    /// @brief get the i-th composition value [#/m^3]
    virtual double operator()(int i) { return ni[i]; }

    /// @brief get the computed composition values as a vector of double values 
    virtual std::vector<double> compositions() { return ni; }

    /// @brief get the computed composition values multiplied by a conversion factor
    virtual std::vector<double> compositions(double conversion) ;

    /// @brief get the computed total composition value 
    virtual double ntot() ;

    /// @brief set starting composition to member ni
    void setn0(std::vector<double> n0) { ni = std::move(n0); }

    /// @brief Get access to the current PfBox (non-const, advanced usage).
    PfBox* getPfBox() { return Qbox; }

    /// @brief set for a desired PF container with desired methods
    void setPfBox(PfBox* specifiedQbox) { Qbox = specifiedQbox; }

    /// @brief Set which Debye model to use
    void setDebyeModel(const std::string& modelName);

private:

    /// @brief Debye length via Rat (2002)
    double Debye_Rat2002(double T);

    /// @brief Debye length via Ghourui formulation
    double Debye_Ghourui(double T);

};

/**
 * @brief Implementation of the Godin–Trépanier matrix-based solver for Saha equilibrium.
 *
 * This class implements the algorithm described in Godin & Trépanier for the computation
 * of equilibrium composition of ionized gas mixtures, under the assumptions of ideal gas
 * state equation and charge neutrality. The algorithm introduces a base/non-base split 
 * for numerical stability.
 *
 * @see <a href="../../articles/A Robust and Efficient Method for the Computation of Equilibrium 
 * Composition in Gaseous Mixtures.pdf">Reference article</a> */
class GodinTrepSahaSolver : public Composition {

    /// @brief Composition Matrix, fixed on given Mixture
    std::vector<std::vector<double>> C;

public:
    
    /// @brief Base constructor
    GodinTrepSahaSolver(Mixture* mix, Gas* gas);

    /// @brief Constructor with qbox assignment.
    GodinTrepSahaSolver(Mixture* mix, Gas* gas, PfBox* qbox);

    /** @brief Solve the composition with the specified algorithm in this class reference
    and store results in member ni. */
    void CompositionSolve(Mixture* mix, Gas* gas) override;

private:

    /// @brief Construct the composition matrix rows for a given specie.
    std::vector<double> Crow(Species* specie, const std::map<std::type_index, int>& colmap);

    /// @brief Build the composition matrix for a given mixture.
    std::vector<std::vector<double>> CompositionMatrix(Mixture& mixx);

    /// @brief Build the conservation matrix for a given mixture and gas state.
    std::vector<std::vector<double>> ConservationMatrix(Mixture& mixx, Gas& gass,
                                                        const std::vector<std::vector<double>>& C);
    
    /// @brief Identify base and non-base species indices as described in the algorithm.
    void baseCalc(std::vector<int>& b, std::vector<int>& bs,
        const std::vector<std::vector<double>>& C);

};

/**
 * @brief CSV printer for composition results.
 *
 * Provides formatted CSV output for any Composition solver,
 * iterating over a given temperature range. */
class CompositionCsv : public DataPrinter {

    /// @brief Solver reference (non-owning, external lifetime managed).
    Composition* solver;

    /// @brief Build the output filename.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepare the CSV header.
    void PrepareHeader() override; 

    /// @brief Prepare the data matrix for printing.
    void PrepareData(const std::vector<double>& x, GasMixture* gasmix) override;

    /// @brief Print completion message.
    void PrintMessage(const std::string& filename) override;

public:

    /// @brief Constructor with explicit solver
    CompositionCsv(Composition* solver);

    /// @brief Constructor with explicit solver and custom folder
    CompositionCsv(Composition* solver, const std::string& folder);
};

#endif
