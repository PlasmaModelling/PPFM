 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef MIXTURE_H
#define MIXTURE_H

#include "AcceptedSpecies.h"
#include "Theta.h"

// forward declarations
class Species;
class Element;
class Electron;


/// @brief Abstract base class for gas mixtures with chemical species.
class Mixture {

    friend class CiBox;
    friend class PfBox;
    friend class Composition;

    protected:

    /// @brief Vector of accepted species (ions, atoms, electrons).
    std::vector<AcceptedSpecies> cispecies;

    /// @brief Tuple containing base species and their mole fractions.
    std::tuple<std::vector<Species*>, std::vector<double>> molefractions;

    /// @brief Number of species in the mixture.
    int N;

    /// @brief Number of elemental species (atoms or electrons).
    int M;

    public:

    /// @brief Constructs a Mixture from a variable number of AcceptedSpecies.
    template<typename... speciestypes>
    Mixture(speciestypes... newspecies);

    /// @brief Returns the number of species.
    int getN();

    /// @brief Returns the number of elemental species.
    int getM();

    /**
     * @brief Returns the i-th species as a base pointer.
     * @param i Index of the species.
     * @return Pointer to the species.
     * @throws std::out_of_range if i is out of bounds.
     */
    Species* operator()(int i);

    /**
     * @brief Returns a vector of species masses.
     * @return Vector of masses [amu].
     */
    std::vector<double> masses();

    /**
     * @brief Returns a vector of converted species masses.
     * @param conversion Multiplicative factor (e.g. amu to kg).
     * @return Vector of converted masses.
     */
    std::vector<double> masses(double conversion);

    /// @brief Virtual destructor.
    virtual ~Mixture() = 0;

    private:

    /// @brief Orders species by their chemical formula.
    void order();

    /**
     * @brief Checks if a species belongs to the base set.
     * @param specie Pointer to the species.
     * @return True if it matches a base species in molefractions.
     */
    bool isBase(Species* specie);

    /**
     * @brief Resets the composition of the mixture.
     * @details Must be implemented in derived classes to reset mole fractions or species states.
     */
    virtual void restartComposition() = 0;

};

/// @brief Basic thermodynamic state: temperature and pressure.
class Gas {

    protected:

    /// @brief Temperature [K].
    double T;

    /// @brief Pressure [Pa].
    double P;

    public:

    /// @brief Pointer to the temperature scaling object (θ).
    Theta* theta;

    /// @brief Constructs a gas object with pressure and temperature.
    Gas(double pressure, double temperature);

    /// @brief Returns the pressure [Pa].
    double getPressure();

    /// @brief Returns the temperature [K].
    double getTemperature();

    /// @brief Sets the temperature [K].
    virtual void setT(double temperature);

    /// @brief Sets the pressure [Pa].
    virtual void setP(double pressure);

    /// @brief Destructor.
    ~Gas();

};

//___________________________________ Implementation ___________________________________

template<typename... speciestypes>
Mixture::Mixture(speciestypes... newspecies) {
    
    (cispecies.push_back(newspecies), ... ) ;
    N = cispecies.size() ;
    M = 0 ;
    for (int i = 0; i < N; i++) {
        if (dynamic_cast<Element*>((*this)(i)) != nullptr || dynamic_cast<Electron*>((*this)(i)) != nullptr ){
            M++ ;
        } 
    }
    order() ; 
    for (int i = 0; i < N ; i++) {
        if (Element* ptr = dynamic_cast<Element*>((*this)(i))){
            std::get<0>(molefractions).push_back((*this)(i)) ;
            std::get<1>(molefractions).push_back((1./((double)M-1.) )) ;
        } 
    }

} 

#endif