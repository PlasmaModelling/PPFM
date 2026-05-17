 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef PARTITION_FUNCTION_H
#define PARTITION_FUNCTION_H

/**
 * @file ParfCalculator.h
 * @brief This file contains the concrete and interface class
 * for partition function in the same fashion as collision integrals.
*/

#include "ParfCalculator.h"
#include "GasMixture.h"
#include <stdexcept>
#include <iostream>
#include <iomanip>

/** @class PFinterface
 *  @brief Interface for Partition Function calculations
 *  @details This class defines the methods required to compute partition functions
 *  for various systems, allowing multiple implementations (e.g., Polynomial, Tabulated). */
class PFinterface {

    public:

    /** @brief Interface method to compute the partition function
     *  @param temperature Temperature [K°]
     *  @param pressure Pressure [Pa]
     *  @param lambdaD Debye length [m] */
    virtual void computePartitionFunction(double temperature, double pressure, double lambdaD) = 0;

    /** @brief Displays information about the partition function specie and method */
    virtual void info() = 0;

    /** @brief Sets the method to compute partition functions
     *  @details set PFtable interpolation method
     *  @see ParfCalculator.h for more details */
    virtual void setTable() = 0;


    /** @brief Sets the method to compute partition functions
     *  @details set AbInitio methods for atomic or molecular partition functions
     *  @see ParfCalculator.h for more details */
    virtual void setAbInitio() = 0;

    /** @brief Sets the method to compute partition functions
     *  @details set AbInitio methods for atomic or molecular partition functions
     *  @see ParfCalculator.h for more details */
    virtual void setAbInitio ( std::string hyperref ) = 0;


    /** @brief Returns the computed partition function value
     *  @return Partition function value */
    virtual double getPf() = 0;

    ///** @brief Returns a pointer to the species associated with this partition function.
    virtual Species* getSp() = 0;

    virtual ~PFinterface() = default;

};

 /** @class PartitionFunction
 *  @brief Partition function calculation class implementing PFinterface methods
 *  @tparam T Template parameter representing the species type
 *  @details This class allows computation of partition functions using different methods
 *  like polynomial approximation, tabulated data from files or Ab-Initio, select the method 
 *  calling setMethod on a PfBox[i] object. */
template <typename T>
class PartitionFunction : public PFinterface {

    friend class PfBox ; 

    protected:
    
    /** @brief Pointer to the calculation method */
    QCalculator* calculator;

    /** @brief Computed partition function value */
    double Q;

    /** @brief Pointer to the species object */
    T* sp;

    /** @brief Concrete template class Partition Function Constructor */
    PartitionFunction(T* Sp) ;
    
    /** @brief Tries to nitialize a calculator among all the possibilities */
    void initCalculator() ; 
    
    /** @brief Computes the partition function using the selected method
     *  @param temperature Temperature [K°]
     *  @param pressure Pressure [Pa]
     *  @param lambdaD Debye length [m] */
    void computePartitionFunction(double temperature, double pressure, double lambdaD) override ;

    /** @brief Sets the method to compute partition functions
     *  @details set PFtable interpolation method
     *  @see ParfCalculator.h for more details */
    void setTable() override ; 

    /** @brief Sets the method to compute partition functions
     *  @details set AbInitio methods for atomic or molecular partition functions
     *  @see ParfCalculator.h for more details */
    void setAbInitio() override ;

    /** @brief Sets the method to compute partition functions
     *  @details set AbInitio methods for atomic or molecular partition functions
     *  @see ParfCalculator.h for more details */
    virtual void setAbInitio ( std::string hyperref ) override ;

    /** @brief Returns the computed partition function value
     *  @return Partition function value */
    double getPf() override ;

    /// @brief Return a pointer to the species associated with this partition function. 
    Species* getSp() { return this->sp ; }

    /** @brief Displays information about the partition function and method used */
    void info() override ;

    public:

    ~PartitionFunction() { delete calculator; delete sp; }

} ;

//________________________ Implementation _____________________________

template <typename T>
void PartitionFunction<T>::initCalculator() {

    sp = new T ;
    
    try
    {
        calculator = new PfPoly(sp) ; 
    }
    catch(const std::exception& e)
    {
        try
        {
            calculator = new PfTtable(sp) ;
        }
        catch(const std::exception& e)
        {
            try
            {
                calculator = new PfTPtable(sp) ;
            }
            catch(const std::exception& e)
            {
                try
                {
                    calculator = new ElectronicAtomicPF(sp) ; 
                }
                catch(const std::exception& e)
                {
                    std::cerr << "No available Partition Function, see ParfCalculator.h for more \n "<<
                        e.what() << '\n';
                }
            }
        }
    }
}

template <typename T>
PartitionFunction<T>::PartitionFunction(T* Sp) { initCalculator(); }

template <typename T>
void PartitionFunction<T>::computePartitionFunction( double temperature, 
    double pressure, double lambdaD ) {
    
    try {
    
        this->Q = calculator->compute(temperature, pressure, lambdaD);
    
    } catch (const std::exception& e) {
       
        throw;
        std::exit(EXIT_FAILURE);

    }
            
}

template <typename T>
void PartitionFunction<T>::setTable() {

    delete calculator;

    try 
    {

        calculator = new PfTtable(this->sp) ;
    
    } catch(const std::exception& e) {

        try {

        calculator = new PfTPtable(this->sp) ;
    
        } catch(const std::exception& e) {

            std::cerr <<"Error in creating PfTable:\n"<< e.what() <<"\nsee ParfCalculator.h for more. ";
        }
    }
}

template <typename T>
void PartitionFunction<T>::setAbInitio() {

    delete calculator;

    if(this->sp->getFormula() != "e-") {

        try {
            
            calculator = new ElectronicAtomicPF(this->sp) ;
            
        } catch(const std::exception& e) {
            
            std::cerr <<"No Table found:\n"<< e.what() << '\n' ;
            std::exit(EXIT_FAILURE); // o std::abort();
            
        }

    }
        
}

template <typename T>
void PartitionFunction<T>::setAbInitio(std::string hyperref) {

    delete calculator;

    try {

        calculator = new ElectronicAtomicPF( this->sp, hyperref) ;
    
    } catch(const std::exception& e) {

        std::cerr <<"No Table found:\n"<< e.what() << '\n' ;
    }

}

template <typename T>
double PartitionFunction<T>::getPf() {
    
    double q = this->Q ; 
    if (q < 0.) {
        
        throw std::invalid_argument(
            "Negative partition function value for "  +this->sp->getFormula()+  
            "!" +std::to_string(q)+ " \nCheck for polynomial coefficients \n"
            "implemented in HCParFun.h in the interested TemperatureRange or\n"
            "try with different method"
        );
    }

    return this->Q;
}

template <typename T>
void PartitionFunction<T>::info() {

    int w = 12;
    std::cout << std::left;
    std::cout << std::setw(w) << this->sp->getFormula() << std::setw(w)
              << this->calculator->metodo << std::endl;
}

#endif
