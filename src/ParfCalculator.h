#ifndef PARF_CALCULATOR_H
#define PARF_CALCULATOR_H

// forward declaration of class Species
#include"DataLoader.h"
#include"Species.h"
#include<vector>
#include<string>
#include <stdexcept>

/** @class QCalculator
 *  @brief Base class for partition function calculation methods
 *  @details This is an abstract class that defines the interface for computing 
 *  partition functions using different methods (e.g., polynomial, tabulated data or Ab-initio Calculations). */
class QCalculator {

    public:

    /** @brief The method name used for computing partition functions */
    std::string metodo;

    /** @brief Pure virtual method to compute the partition function
     *  @param T Temperature [K°]
     *  @param P Pressure [Pa]
     *  @param debye Debye length [m]
     *  @return Computed partition function value */
    virtual double compute(double T, double P, double debye) = 0;

};

/** @class PfPoly
 *  @brief Base class for polynomial partition function computation
 *  @details Implements the polynomial method to compute partition functions based on a 
 *  set of hard-coded polynomial coefficients. This class is designed to be extended by 
 *  specializations for specific chemical species. */
class PfPoly : public QCalculator {

    protected:

    /** @brief Coefficients of a power descending polynomial
     * on Temperature. 
     * @details coefficients[n]*T^n + ... + coefficients[2]*T^2 + 
     * coefficients[1]*T + coefficients[0] */
    std::vector<double> coefficients;

    void Init( std::vector<double> cfs ) { 

        metodo = "Polinomial" ; 
        coefficients = cfs ; 
    
    }

    public:
    
    // Default (fails)
    PfPoly( Species* specieToFail ) { throw std::invalid_argument( "No Available Partition Function for "+specieToFail->getFormula()+" try with a different method.\n" ); }

    // Overload hard-coded partition functions
    PfPoly(Argon* argon ) ;
    PfPoly(ArgonI* argonI ) ;
    PfPoly(ArgonII* argonII ) ;
    PfPoly(Electron* argonII ) ;
    PfPoly(Hydrogen* elettrone ) ;
    PfPoly(MolecularNitrogen* ptr ) ;
    PfPoly(Nitrogen* ptr ) ;
    PfPoly(NitrogenI* ptr ) ;
    PfPoly(NitrogenII* ptr ) ;
    PfPoly(NitrogenIII* ptr ) ;
    PfPoly(MolecularHydrogen* elettrone ) ;
    PfPoly(HydrogenI* elettrone ) ;
    PfPoly(HydrogenAnion* ptr ) ;

    /** @brief Computes Partition Function using a polynomial formula
     *  @param T Temperature [K°]
     *  @param P Pressure [Pa]
     *  @param debye Debye length [m]
     *  @return Computed partition function value */
    double compute(double T, double P, double debye) override ;

};

class PfTtable : public QCalculator, public DataLoader {
    
    /// @brief Temperatures red from the file 
    std::vector<double> temperatures;
    
    /// @brief Partition Functions data only depending on Temperature
    std::vector<double> QQ;

    /** @brief Reads partition function data from temperature-only files
     *  @param file The input file stream is initialized in this class' constructor */
    void ParseFile(std::ifstream& file) override;

    /// @brief Builds the filename based on species name and optional prefix
    std::string BuildFileName(const std::string& name) override;

    void Init() override {}

    public:

    /// @brief Constructor that loads partition function default name datafile
    PfTtable(Species* sp);

    /// @brief Constructor that loads partition function datafile with a custom prefix
    PfTtable(Species* sp, const std::string& prefix);

    /** @brief Computes Partition Function interpolating dedicated files 
     *  @details The interpolation on Temperature is performed via a logaritmic
     *  (positive) cubic spline and the interpolation on Pressure is performed 
     *  via linear interpolation 
     *  @param T Temperature [K°]
     *  @param P Pressure [Pa]
     *  @param debye Debye length [m]
     *  @return Computed partition function value */
    double compute(double T, double P, double debye) override;
};

class PfTPtable : public QCalculator, public DataLoader {
    
    /// @brief Temperatures red from the file     
    std::vector<double> temperatures;
    
    /// @brief Pressures red from the file 
    std::vector<double> pressures;
    
    /// @brief Partition Functions data depending on Temperatures and Pressures
    std::vector<std::vector<double>> QQ2D;

    /// @brief Parses the CSV file for temperature and pressure-dependent partition function data
    void ParseFile(std::ifstream& file) override;

    /// @brief Builds the filename based on species name and optional prefix
    std::string BuildFileName(const std::string& name) override;

    void Init() override {}

    public:
    
    /// @brief Constructor that loads partition function default name datafile
    PfTPtable(Species* sp);

    /// @brief Constructor that loads partition function datafile with a custom prefix
    PfTPtable(Species* sp, const std::string& prefix);

    /** @brief Computes Partition Function interpolating dedicated files 
     *  @details The interpolation on Temperature is performed via a logaritmic
     *  (positive) cubic spline and the interpolation on Pressure is performed 
     *  via linear interpolation 
     *  @param T Temperature [K°]
     *  @param P Pressure [Pa]
     *  @param debye Debye length [m]
     *  @return Computed partition function value */
    double compute(double T, double P, double debye) override;

};

class ElectronicAtomicPF : public QCalculator, public DataLoader {

    /** @brief Stores electronic energy levels: 
     ** first column for degeneracies[#] second for energies [J] */
    std::vector<std::vector<double>> Levels; 
    
    /// @brief Pointer to the species object
    Species* sp;

    /// @brief Downloads and processes data from the NIST Atomic Spectra Database
    void CsvDownloadFromNISTAtomicSpectraDatabase(const std::string& hyperref, const std::string& filename);

    /// @brief Parses the CSV file containing electronic configuration data
    void ParseFile(std::ifstream& file) override;

    /// @brief Builds the filename based on species name and optional prefix
    std::string BuildFileName(const std::string& name) override;

    void Init() override {}
    
    public:
    
    /** @brief Reads the electronic configuration from a datafile placed in 
     *  /data/Electronic_Configuration/ folder, memorize energy levels constructing the object
     *  @details the file has to be named after Species->getFormula()_ElConfig.csv, 
     *  @see ParfCalculator.cpp for more details 
     *  @example Ar_ElConfig.csv */    
    ElectronicAtomicPF ( Species* sp ) ;
    
    /// @brief @see ElectronicAtomicPF(Species* sp); but with customized prefix
    ElectronicAtomicPF ( const std::string& prefix, Species* sp ) ;

    /** @brief Overload constructor to download data from NIST Atomic Spectra Database.
     *  It automatically writes the csv file that should be edited by the user otherwise.
     *  @details Go to https://physics.nist.gov/PhysRefData/ASD/levels_form.html and ask for the
     *  spectrum of the desired species, use spectroscopic nomenclature.\n
     *  set Level Units: eV ,\n 
     *  Format Output: TAB-DELIMETED (text).\n
     *  Check Term ordered.\n 
     *  Uncheck everything but "Level" and "g"
     *  @param sp species
     *  @param hyperref from https://physics.nist.gov/PhysRefData/ASD/levels_form.html 
     *  @example link https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=Ar+I&submit=Retrieve+Data&units=1&format=2&output=0&page_size=15&multiplet_ordered=0&level_out=on&g_out=on&temp= */
    ElectronicAtomicPF ( Species* sp, const std::string& hyperref );

    /** @brief Computes the atomic partition function value for given temperature,
     *  pressure, and Debye length 
     *  using the usual summatory on energy levels as described in the reference: 
     *  Capitelli Fundamental Aspects of Chemical Plasma Physics (FACPP) */
    double compute(double T, double P, double debye) override ;
};

#endif
