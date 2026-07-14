 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef DATAPRINTER_H
#define DATAPRINTER_H

/**
 * @file DataPrinter.h
 * @brief This file contains the DataPrinter 
 * class hierarchy for outputting computed data to CSV files.
 * Each printing class has a composition relation with the related solver.
*/

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

class GasMixture; 
class Composition;
class Thermodynamics;
class DevotoTP;
class ZhangMurphyTP;
class TcsInterface;
class CInterface;
class HybridInterface;
class PFinterface;
class MultiCs;
class CsHolder;
class CsCalculator;
class ThresholdCs;
class AbInitioTcsIntegration;

/**
 * @brief Abstract base class for printing data to CSV files.
 * @details Provides a generic interface for constructing file paths, preparing headers,
 * and formatting data for CSV output. Derived classes must implement specific printing logic.
*/
class DataPrinter {

    protected:

    /// @brief Data matrix to print in Csv file    
    std::vector<std::vector<double>> data ;

    /** @brief Header file to construct in the printable object, 
     * remember to use "," as separator for the printables. 
    */
    std::string header ;

    /** @brief Implement to return a Default filename with coherent prefixes and extensions.  
     * @param filename Reference string for the file name ex: species->getFormula().  
    */
    virtual std::string BuildFileName( const std::string& filename ) const = 0 ;

    /** @brief Return the complete absolute path ( default: PathTo/PPFMfolder/out/" directory ). 
     * Override to specify other path when deriving this class. 
     * * @param filename Reference string for the file name ex: species->getFormula(). */
    virtual std::string BuildFilePath( const std::string& filename ) const ;

    /** @brief Override to implement data calculations and storing in this class member 
    vector<vector<double>> data. 
    * @param x Independent variable 
    * @param gasmix Pointer to class GasMixture  */
    virtual void PrepareData(const std::vector<double>& x, GasMixture* gasmix) = 0 ;

    // Override to implement the header to the data file
    virtual void PrepareHeader() = 0 ;

    /// @brief FileWriting centralized method.
    void WriteData(const std::string& filename ) ; 
    
    /// @brief FileAppendingCentralizedMethod
    void AppendData(const std::string& filename) ;

    public:

    /// @brief customizable folder, just assign it implementing constructor in derived classes
    std::string customFolder ;
    
    /// @brief Centralized printing method, override to customize Printing
    virtual void Print ( const std::string& filename, const std::vector<double>& x, GasMixture* gasmix ) ;

    /// @brief Print a message of finishing the print, override to customize message. 
    virtual void PrintMessage(const std::string& filename) = 0 ; 

};


/**
 * @brief CSV printer for composition results.
 *
 * Provides formatted CSV output for any Composition solver,
 * iterating over a given temperature range. */
class CompositionCsv : public DataPrinter {

    /// @brief Build the output filename.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepare the CSV header.
    void PrepareHeader() override; 
    
    /// @brief Prepare the data matrix for printing.
    void PrepareData(const std::vector<double>& x, GasMixture* gasmix) override;
    
    /// @brief Print completion message.
    void PrintMessage(const std::string& filename) override;
    
    public:

    /// @brief Solver reference (non-owning, external lifetime managed).
    Composition* solver;

    /// @brief Constructor with explicit solver
    CompositionCsv(Composition* solver);

    /// @brief Constructor with explicit solver and custom folder
    CompositionCsv(Composition* solver, const std::string& folder);
};

/** @brief Outputs thermodynamic properties to CSV.
 ** @see class Thermodynamics for property computation.
 ** @see class DataPrinter for CSV interface. */
class ThermodynamicsCsv : public DataPrinter {
    
    /// @brief Build the output filename.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepare the CSV header row.
    void PrepareHeader() override;

    /// @brief Prepare the data matrix for printing.
    void PrepareData(const std::vector<double>& x, GasMixture* gasmix) override;
    
    /// @brief Print completion message.
    void PrintMessage(const std::string& filename) override;
    
    public:
    
    /// @brief Solver reference (non-owning, external lifetime managed).
    Thermodynamics* solver;
    
    /// @brief Constructor with explicit solver
    ThermodynamicsCsv(Thermodynamics* solver);

    /// @brief Constructor with explicit solver and custom folder
    ThermodynamicsCsv(Thermodynamics* solver, const std::string& folder);

};

/// @brief CSV printer for transport properties computed via Devoto’s method.
/// @see class DevotoTP for transport computation.
/// @see class DataPrinter for CSV handling.
class DevotoTpCsv : public DataPrinter {
    
    /// @brief Builds the output filename with "TP_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header row.
    void PrepareHeader() override;

    /// @brief Computes and stores transport properties over a temperature range.
    void PrepareData(const std::vector<double>& temperatureRange, GasMixture* gasmix) override;
    
    /// @brief Prints confirmation message after file generation.
    void PrintMessage(const std::string& filename) override;
    
    public:
    
    /// @brief Solver reference (non-owning).
    DevotoTP* solver;

    /// @brief Constructor with explicit solver.
    DevotoTpCsv(DevotoTP* solver);

    /// @brief Constructor with solver and custom folder.
    DevotoTpCsv(DevotoTP* solver, const std::string& folder);
};


/// @brief CSV printer for transport properties computed via Zhang–Murphy model.
/// @see class ZhangMurphyTP for transport computation.
/// @see class DataPrinter for CSV handling.
class ZhangTpCsv : public DataPrinter {

    /// @brief Builds the output filename with "TP_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header row.
    void PrepareHeader() override;

    /// @brief Computes and stores transport properties over a temperature range.
    void PrepareData(const std::vector<double>& temperatureRange, GasMixture* gasmix) override;

    /// @brief Prints confirmation message after file generation.
    void PrintMessage(const std::string& filename) override;
    
    public:

    /// @brief Solver reference (non-owning).
    ZhangMurphyTP* solver;

    /// @brief Constructor with explicit solver.
    ZhangTpCsv(ZhangMurphyTP* solver);

    /// @brief Constructor with solver and custom folder.
    ZhangTpCsv(ZhangMurphyTP* solver, const std::string& folder);
};
class DeflectionAngleCsv : public DataPrinter {

    /// @brief Small non-owning handle to a printable deflection-angle solver.
    struct PrintableSolver {

        /// @brief Solver able to compute chi(b,E).
        AbInitioTcsIntegration* solver {nullptr};

        /// @brief Threshold solver composed of different chi solvers over energy ranges.
        ThresholdCs* thresholdSolver {nullptr};

        /// @brief Filename suffix used for states, thresholds, elastic branches, etc.
        std::string suffix;

    };

    std::vector<double> impact_parameter;
    std::vector<double> energy;

    /// @brief Current solver used by PrepareData.
    AbInitioTcsIntegration* solver {nullptr};

    /// @brief Current threshold solver used by PrepareData.
    ThresholdCs* thresholdSolver {nullptr};

    /// @brief List of all printable solvers found inside the cross-section structure.
    std::vector<PrintableSolver> printableSolvers;

    /// @brief Recursively collects all printable deflection-angle solvers.
    void CollectPrintableSolvers(
        CsCalculator* calculator,
        const std::string& suffix = ""
    );

    /// @brief Selects the proper AbInitioTcsIntegration branch of a ThresholdCs at energy E.
    AbInitioTcsIntegration* SelectThresholdSolver(
        ThresholdCs* threshold,
        double E
    ) const;

    /// @brief Checks if all branches of a ThresholdCs are printable as deflection angles.
    bool IsPrintableThreshold(ThresholdCs* threshold) const;

    /// @brief Builds the output filename with "CHI_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header row.
    void PrepareHeader() override;

    /** 
     * @brief Computes and stores deflection angle data over energy and impact parameter ranges.
     *
     * The temperature vector and gas mixture pointer are not used here.
     */
    void PrepareData(
        const std::vector<double>& null,
        GasMixture* nullPtr
    ) override;

    /// @brief Prints confirmation message after file generation.
    void PrintMessage(const std::string& filename) override;

public:

    /// @brief Constructor with explicit TcsInterface.
    DeflectionAngleCsv(TcsInterface* tcs);

    /// @brief Constructor with TcsInterface and custom folder.
    DeflectionAngleCsv(TcsInterface* tcs, const std::string& folder);

    /// @brief Set custom energy and impact parameter ranges for deflection angle computation.
    void setParameters(std::vector<double> b, std::vector<double> E) {
        impact_parameter = b;
        energy = E;
    }

    /// @brief Prints all available deflection-angle solvers found recursively.
    void Print(
        const std::string& filename,
        const std::vector<double>& x,
        GasMixture* gasmix
    ) override;

    /// @brief Returns true if no printable deflection-angle solver was found.
    bool ToSkip() const;

};
class TransportCrossSectionCsv : public DataPrinter {

    struct PrintableTcs {

        CsCalculator* calculator {nullptr};

        ThresholdCs* threshold {nullptr};

        /// @brief Calculator to be computed before reading calculator.
        CsCalculator* parentCalculator {nullptr};

        std::string suffix;

        bool inelastic {false};

    };

    /// @brief Pointer to transport cross section interface.
    TcsInterface* solver {nullptr};

    CsCalculator* currentParentCalculator {nullptr};

    /// @brief Current calculator used by PrepareData.
    CsCalculator* currentCalculator {nullptr};

    /// @brief Current threshold calculator, used only to enrich the header.
    ThresholdCs* currentThreshold {nullptr};

    /// @brief True when preparing inelastic data.
    bool currentInelastic {false};

    /// @brief All printable TCS blocks found recursively.
    std::vector<PrintableTcs> printableTcs;

    /// @brief Builds the output filename with "TCS_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header.
    void PrepareHeader() override;

    /// @brief Prepares data for the current calculator.
    void PrepareData(const std::vector<double>& x, GasMixture* gasmix) override;

    /// @brief Prints a confirmation message after writing the file.
    void PrintMessage(const std::string& filename) override;

    /// @brief Recursively collects printable TCS calculators.
    void CollectPrintableTcs(
        CsCalculator* calculator,
        const std::string& suffix = "",
        bool inelastic = false
    );

public:

    /// @brief Constructor with solver.
    TransportCrossSectionCsv(TcsInterface* solver);

    /// @brief Constructor with solver and custom folder.
    TransportCrossSectionCsv(TcsInterface* solver, const std::string& folder);

    /// @brief Overrides default print to handle MultiCs, ThresholdCs and CsHolder recursively.
    void Print(
        const std::string& filename,
        const std::vector<double>& x,
        GasMixture* gasmix
    ) override;

};

/// @brief CSV printer for Collision Integrals.
/// @details For printing single collision integrals not referring to a GasMixture
/// call for Print with nullptr as GasMixture pointer (third argument).
/// @see class DataPrinter for CSV handling.
class CollisionIntegralCsv : public DataPrinter {

    /// @brief Builds the output filename with "CI_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header.
    void PrepareHeader() override;

    /// @brief Prepares data across temperatures for Ω(l,s) integrals.
    void PrepareData(const std::vector<double>& x, GasMixture* gasmix) override;

    /// @brief Prints a confirmation message after writing the file.
    void PrintMessage(const std::string& filename) override;

    public:
    
    /// @brief Pointer to collision integral solver.
    HybridInterface* solver;

    /// @brief Constructor with solver.
    CollisionIntegralCsv(HybridInterface* solver);

    /// @brief Constructor with solver and custom output folder.
    CollisionIntegralCsv(HybridInterface* solver, const std::string& folder);

};

/// @brief Outputs partition functions computed from PFinterface to CSV.
/// @see class DataPrinter for CSV handling.
/// @see class PFinterface for partition function computation.
class PartitionFunctionCsv : public DataPrinter {

    /// @brief Builds the output filename with "PF_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header row.
    void PrepareHeader() override;

    /// @brief Computes and stores partition function data.
    void PrepareData(const std::vector<double>& Ti, GasMixture* gasmix) override;

    /// @brief Prints a confirmation message after writing the file.
    void PrintMessage(const std::string& filename) override;

    public:
    
    /// @brief Pointer to the partition function solver.
    PFinterface* solver;

    /// @brief Constructor with solver pointer.
    PartitionFunctionCsv(PFinterface* solver);

    /// @brief Constructor with solver pointer and custom folder.
    PartitionFunctionCsv(PFinterface* solver, const std::string& folder);

    /// @brief Redirects print to a dedicated subfolder.
    void Print(const std::string& filename, const std::vector<double>& x, GasMixture* gasmix) override;
};

#endif