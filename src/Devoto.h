#ifndef DEVOTO_H
#define DEVOTO_H

#include "Transport.h"
#include "DataPrinter.h"

/** @class DevotoTP
 ** @brief Implementation of the theory formalized by R.S Devoto in \n
 ** @see 1) R.S. Devoto "Transport Properties of Ionized Monoatomic Gases", \n 
 ** The Physics of fluid, 9,6,June(1966) \n
 ** and \n
 ** @see 2) R.S. Devoto "Simplified Expressions for the Transport Properties 
 ** of Ionized Monoatomic Gases", \n  The Physics of fluid, 10,10,October(1967) */
class DevotoTP : public Appendix {

    private : 
    
    /**  @brief Function to compute Electrons Thermal Conductivity as in eq.20 of
     ** this class reference 2.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation */
    double ThermalCondEl ( GasMixture* gasmix, int order ) override ;
    double TotalThermalCondEl ( GasMixture* gasmix, int order ) ; 

    /**  @brief Function to compute Electrons Thermal Conductivity as in eq.14 of
     ** this class reference 1.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation  */
    double ThermalCondHeavy ( GasMixture* gasmix, int order ) override ;
    double TotalThermalCondHeavy ( GasMixture* gasmix, int order ) ; 

    /** @brief Function to compute Viscosity as in eq.21 of
     ** this class reference 1.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation  */
    double Viscosity ( GasMixture* gasmix, int order ) override ;
    

    /** @brief Function to compute Electrical Conductivity as in eq.16 of
     ** * this class reference 2.
     ** @param gasmix GasMixture object 
     ** @param order Desired order of approximation 
    */
    double ElCond ( GasMixture* gasmix, int order ) override ;
    
    /** @brief Function to compute Heat Exchange of electrons-heavy collisions
     ** @param gasmix GasMixture object  */
    double Qeh ( GasMixture* gasmix ) override ; 
    
    /** @brief Function to compute Diffusion Coefficients as in eq.8 of
     ** this class reference 1.
     ** @param gasmix GasMixture object
     ** @param order Desired order of approximation
     ** @param i-th specie 
     ** @param j-th specie */
    double Dij ( GasMixture* gasmix , int order , int i , int j ) override ; 

    /** @brief Function to compute Thermal Diffusion Coefficients as in eq.9 of
     ** this class reference 1.
     ** @param gasmix GasMixture object
     ** @param order Desired order of approximation
     ** @param i-th specie  */
    double DiT ( GasMixture* gasmix , int order , int i ) override ;

    public : 

    DevotoTP (  CiBox* cbx  ) : Appendix ( cbx ) {} 
    DevotoTP ( GasMixture* mix ) : Appendix ( mix ) {} 
    
    /** @brief Function that computes and store the Transport Coefficients in 
     ** the base class Transport members
     ** @param gasmix GasMixture object */
    void computeTransport ( GasMixture* gasmix ) override ;

} ;

//______________________________ Implementation ______________________________

/// @brief Handles CSV output of transport properties using Devoto's method.
class DevotoTpCsv : public DevotoTP, public DataPrinter {

    private:
    
    /// @brief Builds the output filename with "TP_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header row.
    void PrepareHeader();

    /**
     * @brief Computes and stores transport properties over a temperature range.
     * @param temperatureRange Electron temperatures (in θ units).
     * @param gasmix Pointer to the gas mixture.
     * @details For each temperature, computes λₑ, λₕ, μ, σ using Devoto's method.
     * Restores initial gas state after the loop.
     * @see class DevotoTP for computeTransport. */
    void PrepareData(const std::vector<double>& temperatureRange, GasMixture* gasmix) override;

    /**
     * @brief Prints confirmation message after file generation.
     * @param filename Name of the written CSV file. */
    void PrintMessage(const std::string& filename) override;

    public:

    /// @brief Constructor with CiBox pointer.
    DevotoTpCsv(CiBox* cbx);

    /// @brief Constructor with CiBox and custom folder.
    DevotoTpCsv(CiBox* cbx, const std::string& folder);

};

#endif