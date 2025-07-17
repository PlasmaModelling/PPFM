#ifndef ZMTRANSPORT_H
#define ZMTRANSPORT_H

#include"ZMCoefficients.h"
#include"DataPrinter.h"

/** @class ZhangMurphyTP
 * @brief Implementation of the theory formalized by Zhang et. al. in \n
 * @see 1) Xiao-Ning Zhang, He-Ping Li, Anthony B. Murphy, Wei-Dong Xia, \n 
 * "A Numerical Model of Non-Equilibrium Thermal Plasmas. \n 
 * I. Transport Properties," Phys. Plasmas, vol. 20, p. 033508, 2013.
 * @see 2) Comparison of the transport properties of 
 * two-temperature argon plasmas calculated using different methods. */
class ZhangMurphyTP : public ZMCoefficients {

    private:
    
    std::vector<std::vector<double>> DiffTheta ;
    std::vector<double> Dtheta ;

    /** @brief Function to compute Electrons Thermal Conductivity as in eq.34 of
     * this class reference 1.
     * @param gasmix GasMixture object 
     * @param order Desired order of approximation
     * @attention To confront electron thermal conductivity 
     * with other theories it has to be divided by theta as stated in eq.14 
     * of this class reference 2 */
    double ThermalCondEl    ( GasMixture* gasmix, int order ) override ;
    
    /** @brief Function to compute heavy particle Thermal Conductivity as in eq.34 of
     * this class reference 1.
     * @param gasmix GasMixture object 
     * @param order Desired order of approximation */
    double ThermalCondHeavy ( GasMixture* gasmix, int order ) override ;
    
    /** @brief Function to compute Viscosity as in eq.A.18 of
     * this class reference 2.
     * @param gasmix GasMixture object 
     * @param order Desired order of approximation */
    double Viscosity        ( GasMixture* gasmix, int order ) override ;

    /** @brief Function to compute Electrical Consuctivity as in eq.42 of
     * this class reference 1.
     * @param gasmix GasMixture object 
     * @param order Desired order of approximation */
    double ElCond           ( GasMixture* gasmix, int order ) override ;

    double Qeh              ( GasMixture* gasmix )            override { return 0. ; } ; 

    /** @brief Function to compute Ordinary Diffusion Coefficients as in eq.20 of
     * this class reference 1.
     * @param gasmix GasMixture object 
     * @param order Desired order of approximation */
    double Dij ( GasMixture* gasmix , int order , int i , int j ) override ;

    /** @brief Function to compute Thermal Diffusion Coefficients as in eq.21 of
     * this class reference 2.
     * @param gasmix GasMixture object 
     * @param order Desired order of approximation */
    double DiT ( GasMixture* gasmix , int order , int i ) override ;

    /** @brief Function to compute n.e. Ordinary Diffusion Coefficients as in eq.22 of
     * this class reference 1.
     * @param gasmix GasMixture object 
     * @param order Desired order of approximation */
    double DijTheta ( GasMixture* gasmix , int order , int i , int j ) ;
    
    /** @brief Function to compute n.e. thermal Diffusion Coefficients as in eq.36 of
     * this class reference 1.
     * @param gasmix GasMixture object 
     * @param order Desired order of approximation */
    double DiTheta ( GasMixture* gasmix , int order , int i ) ;
    
    /** @brief Function to compute n.e. Electrons Thermal Conductivity as in eq.35 of
     * this class reference 1.
     * @param gasmix GasMixture object 
     * @param order Desired order of approximation
     * @attention To confront electron thermal conductivity 
     * with other theories it has to be divided by theta as stated in eq.14 
     * of this class reference 2 */
    double NeThermalCondEl    ( GasMixture* gasmix, int order ) ;
    
    /** @brief Function to compute n.e. Heavy particle Thermal Conductivity as in eq.35 of
     * this class reference 1.
     * @param gasmix GasMixture object 
     * @param order Desired order of approximation */
    double NeThermalCondHeavy ( GasMixture* gasmix, int order ) ;
    
    public:
    
    ZhangMurphyTP ( CiBox* cbx  ) : ZMCoefficients ( cbx ) {} ;
    ZhangMurphyTP ( GasMixture* mix ) : ZMCoefficients ( mix ) {} 

    /** @brief Function that computes and store the Transport Coefficients in 
     * the base class Transport members
     * @param gasmix GasMixture object */
    void computeTransport ( GasMixture* gasmix ) override ;

    
};

//______________________________ Implementation ______________________________

/// @brief Outputs transport properties computed via Zhang-Murphy model to CSV.
/// @see class ZhangMurphyTP for transport computation.
/// @see class DataPrinter for CSV handling.
class ZhangTpCsv : public ZhangMurphyTP, public DataPrinter {

    /// @brief Builds the output filename with "TP_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header row.
    void PrepareHeader();

    /**
     * @brief Computes and stores transport properties over a temperature range.
     * @param temperatureRange Vector of Te/θ values.
     * @param gasmix Pointer to the gas mixture.
     * @details For each temperature, computes the Zhang-Murphy transport properties:
     * λₑ, λₕ, μ, σ, λₑθ, λₕθ. Initial gas state is restored after the loop.
     * @see ZhangMurphyTP::computeTransport */
    void PrepareData(const std::vector<double>& temperatureRange, GasMixture* gasmix) override;

    /**
     * @brief Prints a confirmation message after writing the file.
     * @param filename Name of the CSV file written. */
    void PrintMessage(const std::string& filename) override;

    public:

    /// @brief Constructor with pointer to collision integrals.
    ZhangTpCsv(CiBox* cbx);

    /// @brief Constructor with CiBox and custom output folder.
    ZhangTpCsv(CiBox* cbx, const std::string& folder);

};

#endif