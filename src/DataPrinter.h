 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef DATAPRINTER_H
#define DATAPRINTER_H

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
class PFinterface;
class MultiCs;
class CsHolder;

class DataPrinter {

protected:

    /// @brief Data matrix to print in Csv file    
    std::vector<std::vector<double>> data ;

    /** @brief Header file to construct in the printable object, 
     remember to use "," as separator for the printables.  */
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

    virtual void PrintMessage(const std::string& filename) = 0 ; 

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

/// @brief CSV printer for transport properties computed via Devoto’s method.
/// @see class DevotoTP for transport computation.
/// @see class DataPrinter for CSV handling.
class DevotoTpCsv : public DataPrinter {

    /// @brief Solver reference (non-owning).
    DevotoTP* solver;

    /// @brief Builds the output filename with "TP_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header row.
    void PrepareHeader() override;

    /// @brief Computes and stores transport properties over a temperature range.
    void PrepareData(const std::vector<double>& temperatureRange, GasMixture* gasmix) override;

    /// @brief Prints confirmation message after file generation.
    void PrintMessage(const std::string& filename) override;

public:

    /// @brief Constructor with explicit solver.
    DevotoTpCsv(DevotoTP* solver);

    /// @brief Constructor with solver and custom folder.
    DevotoTpCsv(DevotoTP* solver, const std::string& folder);
};


/// @brief CSV printer for transport properties computed via Zhang–Murphy model.
/// @see class ZhangMurphyTP for transport computation.
/// @see class DataPrinter for CSV handling.
class ZhangTpCsv : public DataPrinter {

    /// @brief Solver reference (non-owning).
    ZhangMurphyTP* solver;

    /// @brief Builds the output filename with "TP_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header row.
    void PrepareHeader() override;

    /// @brief Computes and stores transport properties over a temperature range.
    void PrepareData(const std::vector<double>& temperatureRange, GasMixture* gasmix) override;

    /// @brief Prints confirmation message after file generation.
    void PrintMessage(const std::string& filename) override;

public:

    /// @brief Constructor with explicit solver.
    ZhangTpCsv(ZhangMurphyTP* solver);

    /// @brief Constructor with solver and custom folder.
    ZhangTpCsv(ZhangMurphyTP* solver, const std::string& folder);
};

/// @brief CSV printer for Transport Cross Sections.
/// @details Transport Cross Sections are independent from GasMixture, place yours 
/// as the third argument of Print or use nullptr.
/// @see class DataPrinter for CSV handling.
class TransportCrossSectionCsv : public DataPrinter {

    /// @brief Pointer to transport cross section solver.
    TcsInterface* solver;

    /// @brief Builds the output filename with "TCS_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header.
    void PrepareHeader() override;

    /// @brief Prepares the CSV header for inelastic transport cross sections.
    void PrepareInelasticHeader();

    /// @brief Prepares data for elastic cross sections.
    void PrepareData(const std::vector<double>& x, GasMixture* gasmix) override;

    /// @brief Prepares elastic data from a CsHolder.
    void PrepareElasticData(CsHolder* tcsElIn);

    /// @brief Prepares inelastic data from a CsHolder.
    void PrepareInelasticData(CsHolder* tcsElIn);

    /// @brief Prepares data for MultiCs objects.
    void PrepareMultiCsData(MultiCs* multiCs, size_t i);

    /// @brief Prints a confirmation message after writing the file.
    void PrintMessage(const std::string& filename) override;

public:

    /// @brief Constructor with solver.
    TransportCrossSectionCsv(TcsInterface* solver);

    /// @brief Constructor with solver and custom folder.
    TransportCrossSectionCsv(TcsInterface* solver, const std::string& folder);

    /// @brief Overrides default print to handle MultiCs and CsHolder cases.
    void Print(const std::string& filename, const std::vector<double>& x, GasMixture* gasmix) override;
};

/// @brief CSV printer for Collision Integrals.
/// @details For printing single collision integrals not referring to a GasMixture
/// call for Print with nullptr as GasMixture pointer (third argument).
/// @see class DataPrinter for CSV handling.
class CollisionIntegralCsv : public DataPrinter {

    /// @brief Pointer to collision integral solver.
    CInterface* solver;

    /// @brief Builds the output filename with "CI_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header.
    void PrepareHeader() override;

    /// @brief Prepares data across temperatures for Ω(l,s) integrals.
    void PrepareData(const std::vector<double>& x, GasMixture* gasmix) override;

    /// @brief Prints a confirmation message after writing the file.
    void PrintMessage(const std::string& filename) override;

public:

    /// @brief Constructor with solver.
    CollisionIntegralCsv(CInterface* solver);

    /// @brief Constructor with solver and custom output folder.
    CollisionIntegralCsv(CInterface* solver, const std::string& folder);

};

/// @brief Outputs partition functions computed from PFinterface to CSV.
/// @see class DataPrinter for CSV handling.
/// @see class PFinterface for partition function computation.
class PartitionFunctionCsv : public DataPrinter {

    /// @brief Pointer to the partition function solver.
    PFinterface* solver;

    /// @brief Builds the output filename with "PF_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header row.
    void PrepareHeader() override;

    /// @brief Computes and stores partition function data.
    void PrepareData(const std::vector<double>& Ti, GasMixture* gasmix) override;

    /// @brief Prints a confirmation message after writing the file.
    void PrintMessage(const std::string& filename) override;

public:

    /// @brief Constructor with solver pointer.
    PartitionFunctionCsv(PFinterface* solver);

    /// @brief Constructor with solver pointer and custom folder.
    PartitionFunctionCsv(PFinterface* solver, const std::string& folder);

    /// @brief Redirects print to a dedicated subfolder.
    void Print(const std::string& filename, const std::vector<double>& x, GasMixture* gasmix) override;
};


#endif