// PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
// (University of Bologna, Italy)                   // 
// Licensed under CC BY 4.0.                        // 
// To view a copy of this license, visit:           // 
// https://creativecommons.org/licenses/by/4.0/     // 

#ifndef COMPOSITION_H
#define COMPOSITION_H

#include "PfBox.h"

#include <map>
#include <string>
#include <typeindex>
#include <vector>
#include <cmath>
/* USE OPTIONAL ONLY TO AVOID CompositionSolve DUPLICATION FOR FROZEN LAMBDA
** TRY TO REASON WITH OVERLOADS IN EVERY OTHER SECTION WITHIN PPFM */
#include <optional>

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
        // and related new private function.

    };

    /// @brief Selected Debye length formulation, defaults to Rat2002Th
    DebyeModel debyeChoice = DebyeModel::Ghourui;

    /** @brief Effective pressure used in the Saha equations. 
     ** Peff = gasmix.getPressure() in the ideal cases. 
     ** @details Peff has nothing to do with the fixed pressure of the 
     ** gas mixture, it is the pressure particles feel inside the plasma. 
     ** remember to update it in solving methods. */
    double Peff; 

    /** @brief Formation energies
     ** @details Null for electrons and neutral atoms
     ** equal to the ionization energy for positive ions,
     ** bond dissociation energy and electron affinity with negative sign for molecules,
     ** and anions, respectively. To be lowered for corrections to ionization energies. */
    std::vector <double> epsf;

    /// @brief Debye length via Rat (2002)
    double Debye_Rat2002(double T);

    /// @brief Debye length via Ghourui formulation
    double Debye_Ghourui(double T);

    /// @brief Initialize composition values 
    virtual void initializeComposition();

    /// @brief Initialize formationEnergies
    virtual void initializeFormationEnergies();

    /// @brief Compute Total Partition Functions
    std::vector<double> totalPartitionFunctions(double T, double P, double lambdaD);

    virtual void restart() = 0;

    virtual ~Composition() = default;

    public:

    /// @brief Base constructor
    Composition(Mixture* mix, Gas* gas) ;

    /// @brief Constructor with qbox assignment
    Composition(Mixture* mix , Gas* gas, PfBox* qbox) ;

    /** @brief Compute and store in ni the particle densities of the Mixture at the Gas state.
     ** USE @p lambda only in advanced composition classes. */
    virtual void CompositionSolve(Mixture* mix, Gas* gas, 
        std::optional<double> lambda = std::nullopt) = 0;

    /** @brief Get the Debye Length of the computed composition.
     ** @param temperature Electron or heavy-particle temperature [K].
     ** @return Debye length [m]. */
    virtual double getDebyeLength(double T) ;

    /// @brief get the i-th composition value [#/m^3]
    virtual double operator()(int i) {return ni[i];}

    /// @brief get the computed composition values as a vector of double values 
    virtual std::vector<double> compositions() {return ni;}

    /// @brief get the computed composition values multiplied by a conversion factor
    virtual std::vector<double> compositions(double conversion) ;

    /// @brief get the computed total composition value 
    virtual double ntot() ;

    /// @brief set starting composition to member ni
    void setn0(std::vector<double> n0) {ni = std::move(n0);}

    /// @brief Get access to the current PfBox (non-const, advanced usage).
    PfBox* getPfBox() {return Qbox;}

    /// @brief set for a desired PF container with desired methods
    void setPfBox(PfBox* specifiedQbox) {Qbox = specifiedQbox;}

    /// @brief Set which Debye model to use
    void setDebyeModel(const std::string& modelName);

};

/** @brief Implementation of the Godin–Trépanier matrix-based solver for Saha equilibrium.
 **
 ** This class implements the algorithm described in Godin & Trépanier for the computation
 ** of equilibrium composition of ionized gas mixtures, under the assumptions of ideal gas
 ** state equation and charge neutrality. 
 ** The algorithm combines Newton method with Mass Action Law
 **
 ** @see <a href="../../articles/A Robust and Efficient Method for the Computation of Equilibrium 
 ** Composition in Gaseous Mixtures.pdf">Reference article</a> */
class GodinTrepSahaSolver : public Composition {
    
    public:
    
    /// @brief Base constructor
    GodinTrepSahaSolver(Mixture* mix, Gas* gas);
    
    /// @brief Constructor with qbox assignment.
    GodinTrepSahaSolver(Mixture* mix, Gas* gas, PfBox* qbox);
    
    /** @brief Solve the composition with the specified algorithm in this class reference
     ** and store results in member ni. */
    void CompositionSolve(Mixture* mix, Gas* gas, std::optional<double> lambda) override;
    
    /** @brief Override of the CompositionSolve in which Debye Length is not 
    updated during Newton steps */
    void CompositionSolveLambdaFrozen( Mixture* mix, Gas* gas, double lambda );

    protected:

    /// @brief Number of species and elements in the mixture, respectively and L = N-M. 
    int N, M, L; 

    /// @brief Composition Matrix, fixed on given Mixture
    std::vector<std::vector<double>> C;

    /** @brief Matrixes required for the algorithm, 
     ** built in constructor only for performances */
    std::vector<std::vector<double>> B,Bs,Binv,v ; 

    /** @brief int required for the algorithm, 
     ** built in constructor only for performances */
    std::vector<int> b, bs;

    /** @brief The RHS of the conservation matrix as described in the algorithm, 
     ** VERY IMPORTANT, as it defines the ntot of the mixture */
    std::vector<double> A0;

    /** @brief Construct the composition matrix rows for a given specie 
     ** polimorphically casting its nature. See Species.h for semantic definitions 
     ** of chemical species */
    std::vector<double> Crow(Species* specie, const std::map<std::type_index, int>& colmap);

    /// @brief Build the composition matrix for a given mixture as specified in the reference.
    std::vector<std::vector<double>> CompositionMatrix(Mixture& mixx);

    /// @brief Build the conservation matrix for a given mixture and gas state.
    std::vector<std::vector<double>> ConservationMatrix(Mixture& mixx, Gas& gass,
        const std::vector<std::vector<double>>& C);
    
    /// @brief Identify base and non-base species indices as described in the reference.
    void baseCalc(std::vector<int>& b, std::vector<int>& bs,
        const std::vector<std::vector<double>>& C);

    virtual void restart() override;

};

/** @brief Extension of GodinTrepSahaSolver with Debye–Hückel corrections.
 ** @details Implements the iterative Debye–Hückel correction to the equilibrium composition
 ** by updating the effective pressure and the formation energies of charged species
 ** until self-consistency between Debye length and total density is reached.
 ** References:
 ** - Capitelli Fundamental Aspects of Chemical Plasma Physics, 2016,
 ** - Godin, A., & Trépanier, J.-Y., (2002), *J. Thermophys. Heat Transfer*, 16(1), 145–152. */
class GTSahaDHcorrection : public GodinTrepSahaSolver {

    public:

    /// @brief Base constructor.
    GTSahaDHcorrection(Mixture* mix, Gas* gas);

    /// @brief Constructor with PfBox assignment.
    GTSahaDHcorrection(Mixture* mix, Gas* gas, PfBox* qbox);

    /// @brief Solve the composition including Debye–Hückel corrections.
    void CompositionSolve(Mixture* mix, Gas* gas, std::optional<double> lambda) override;

    protected:

    /// @brief Compute Debye–Hückel correction to pressure [Pa].
    double pressureDHcorrected(double T, double lambdaD) ;

    /// @brief Apply DH correction to formation energies (species with z ≠ 0).
    std::vector<double> formationEnergyDHcorrected(double lambdaD);

    /// @brief Convergence tolerance for Debye iteration.
    double tol = 1e-6;

    /// @brief Maximum number of iterations for self-consistent DH correction.
    int maxIter = 500;

    virtual void restart() override;

};

#endif
