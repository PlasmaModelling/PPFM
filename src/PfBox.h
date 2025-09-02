 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef PFBOX_H 
#define PFBOX_H

class PFinterface ;
class GasMixture ;
class Mixture ;
class Species ; 

#include <string>
#include <vector>
/// @brief Manages and computes partition functions for all species in a mixture.
/// @see class PFinterface for partition function model.
/// @see class PartitionFunctionCsv for CSV output.
class PfBox {

    /// @brief Vector of pointers to partition function interfaces.
    std::vector<PFinterface*> partitionfunctions;

    public:

    /**
     * @brief Constructs PfBox and initializes one partition function per species.
     * @param mix Pointer to the species container (Mixture). */
    PfBox(Mixture* mix);

    /**
     * @brief Accesses the i-th PFinterface.
     * @param i Index of the partition function.
     * @return Pointer to the i-th PFinterface.
     * @throws std::out_of_range if index is invalid. */
    PFinterface* operator[](int i);

    /**
     * @brief Returns the value of the i-th partition function.
     * @param i Index of the partition function.
     * @return Q value of the partition function.
     * @throws std::out_of_range if index is invalid. */
    double operator()(int i);

    /**
     * @brief Prints all partition functions to CSV files.
     * @param Ti Vector of temperatures [K].
     * @param gasmix Pointer to the gas mixture.
     * @param folder Output folder for CSV files.
     * @details For each species, constructs a `PartitionFunctionCsv` writer and 
     * prints its Q values over the given temperature range. */
    void PrintPartitionFunctions(const std::vector<double>& Ti, GasMixture* gasmix, const std::string& folder);

    /**
     * @brief Computes all partition functions in parallel.
     * @param temperature Temperature [K].
     * @param pressure Pressure [Pa].
     * @param lambda Debye length [m].
     * @details Uses OpenMP to parallelize the computation. If any function throws, 
     * an error message is printed and execution is terminated.
     * @throws std::exit(EXIT_FAILURE) if any computation fails. */
    void computePartitionFunctions(double temperature, double pressure, double lambda);

    /// @brief Displays a summary of partition functions and associated methods.
    void info();

};


#endif