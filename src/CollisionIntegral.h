 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

/** 
 * @file CollisionIntegral.h
 * @brief This file contains the declaration of the Collision Integral related classes.
 * @details COMPUTATION OF COLLISION INTEGRALS AS A FUNCTION OF T(K°)    \n
 * reference: \n
 * R. S. Devoto, "Transport Properties of Ionized Monatomic Gases",
 * Physics of Fluids, vol. 9, n. 6, pp. 1230–1240, June 1966, DOI: 10.1063/1.1761885.
 * Formula 11.
*/

#ifndef COLLISION_INTEGRAL_H
#define COLLISION_INTEGRAL_H


#include"TransportCrossSection.h"
#include"OmegaCalculator.h"


#include <regex>

/** 
 * @class CInterface
 * @brief Interface class only to operate on collision integrals.
 * @details As collision integrals, transport cross sections, deflection angles,
 * and potentials are different mathematical objects they're kept separated
 * in order to operate on them with different methods and algorithms. \n
*/
class CInterface {

    friend class Transport ; 
    friend class CiBox ; 
    friend class CollisionIntegralCsv ; 

    protected: 

    /** @brief Collision Integral Ω(T) needed for the 4th order 
     * Chapman-Enskog approximation to compute transport properties
     * @details Needed values of the collision integrals are : \n
     * (1,1),..,(1,7) \n
     * (2,2),..,(2,6) \n
     * (3,3),..,(3,5) \n
     * (4,4)
     * so omega4th will be 16 elements long
    */
    std::vector<double> omega4th;

    /// @brief Effective temperature of collision between species i and j.
    double TijStar;

    public: 

    /** @brief Pointer to the OmegaCalculator object
     * @details A pointer to the OmegaCalculator hierarchy, 
     * which classes are the actual objects implementing the algorithms. 
    */
    OmegaCalculator* Omega;
    
    /** @brief Interface function to compute the Collision Integral
     * @param temperature Temperature of collision [K°]
     * @param lambda Debye length [m] 
    */
    virtual void ComputeCollisionIntegral ( double Te, double Th, double lambda ) = 0 ;
    
    /** @brief Interface function to write the name of the interaction
     * @return std::string  
    */
    virtual std::string InteractionName() = 0 ;
    
    /** @brief Interface function to set the Interaction Potential 
     * @param pott Potential object
     * @details Assign the Potential to the Collision Integral of a Non Coulomb interaction.
    */
    virtual void Pot ( Potential* pott ) = 0;
    
    /** @brief Interface function to set the Interaction Potentials for multi-state interactions. 
     * @param potentials std::initializer_list<Potential*> of { Potential objects }
     * @param gs std::initializer_list<double> of { State degeneracies }
     * @see class Potential and derived for further details on implemented Potentials 
    */
    virtual void MultiPot ( std::initializer_list<Potential*> potentials ,
        std::initializer_list<double> gs ) = 0 ;
    
    /** @brief Interface function to set parameters for the Charge Exchange collision cross 
     ** section Q(g)=(A-Bln(g))^2 
     ** @param A [Ang]
     ** @param B [Ang] 
    */
    virtual void ChargeTransfer ( double A, double B ) = 0 ; 
    
    /// @brief Info on the collision integrals of the Mixture 
    virtual void info() = 0 ;
    
    /// @brief Interface function to call to load the collision integral from its raw file
    virtual void Load( bool b ) = 0 ;

} ;


/// @brief An hybrid class to represent a mixed interface for Collision Integrals and Transport Cross Sections.
class HybridInterface : public virtual CInterface, public virtual TcsInterface {};

/**
 * @class CollisionIntegral
 * @brief A concrete template class to represent the collision integral for a pair of species. 
 * @tparam T1 the first chemical specie of the interaction, must be derived from Species.
 * @tparam T2 the second chemical specie of the interaction, must be derived from Species.
 * @details CollisionIntegral inherits from TransportCrossSection to have access to the 
 * TCS calculator and from HybridInterface to have access to the Collision Integral interface. \n
 * The class is designed to compute the collision integrals needed for the Chapman-Enskog method,
 * and to store them in the omega4th member. \n
 * The class also provides methods to set the interaction potential, 
 * to set the parameters for charge exchange interactions, and to print info on the collision integrals.
*/
template<typename T1, typename T2>
class CollisionIntegral : public TransportCrossSection<T1,T2>, public HybridInterface {
        
    /// @brief initialize the object
    void InitOmega4th() { omega4th.resize(16) ; }

    /// @brief Initialize the right OmegaCalculator object, @see OmegaCalculator.h for more 
    void InitCalculator() ;

    /// @brief Default init with no loader overload
    void InitCalculator(bool b) ;

    public:

    /// @brief Constructor for the CollisionIntegral class
    /// @param t1 Pointer to the first chemical specie
    /// @param t2 Pointer to the second chemical specie
    CollisionIntegral(T1* t1,T2* t2) : TransportCrossSection<T1,T2>(t1,t2) {InitOmega4th(); InitCalculator();}

    /** @brief Compute collision integral and store it in the omega4th member of this class.
     * this will be the preferred method if the file is available. \n
     * @param temperatura  [K°]
     * @param lambda Debye Length [m] 
     * @details calls the computation on the OmegaCalculator object for the required collision integrals.
     * implement in case of further orders (l,s)
    */
    void ComputeCollisionIntegral(double Te, double Th, double lambda) override ;

    /// @brief Function returning "Specie1_Specie2" interaction name.
    std::string InteractionName() override ;
        
    /** @brief Interface function to set the Interaction Potential 
     * @param pott Potential object
     * @details Assign the Potential to the Collision Integral of a Non Coulomb interaction. 
    */
    void Pot ( Potential* pott ) ;
    
    /** @brief Interface function to set the Interaction Potentials for multiple interactions 
     * @details Initialize a MultiPotOmega and a MultiCs new objects to members Omega and TCScalculator 
     * @param potentials Initializer_list of { Potential objects }
     * @param gs Initializer_list of { interaction statistical degeneracies }
     * @see class Potential and derived for further details on implemented Potentials 
    */
    void MultiPot ( std::initializer_list<Potential*> potentials ,
        std::initializer_list<double> gs ) override ;

    /** @brief Interface function to set the degeneracies for multi state interactions, 
     * @details Initialize a MultiPotOmega and a nullptr MultiCs new objects to members Omega and TCScalculator.
     * CsCalculators of the MultiCs have to be specified or initialized by the user.   
     * @param gs Initializer_list of { interaction statistical degeneracies } */
    /* virtual void MultiState ( std::initializer_list<double> gs ) override ; */

    /** @brief Interface function to set parameters for the Charge Exchange collision cross 
     ** section Q(g)=(A-Bln(g))^2 
     ** @param A [Ang]
     ** @param B [Ang] */
    void ChargeTransfer ( double A, double B ) override ; 
    
    /// @brief Info on the Collision Integral
    void info() override ; 
    
    /// @brief Try LoaderOmega if argument is true, default if false
    void Load( bool b ) override { b ? InitCalculator() : InitCalculator(b) ;}

};

//________________________________ Implementation ________________________________

template<typename T1 ,typename T2 >
void CollisionIntegral<T1, T2>::InitCalculator() {

    if constexpr ( std::is_base_of_v<ChargedSpecies, T1> && std::is_base_of_v<ChargedSpecies, T2> ) {

        InitCalculator(true);
        return;

    } else {

        try {

            Omega = new LoaderOmega(this);

        } catch (const std::exception&) {

            InitCalculator(true);

        }
    }
}

template<typename T1 ,typename T2 >
void CollisionIntegral<T1, T2>::InitCalculator(bool) {

    if constexpr ( std::is_base_of_v<ChargedSpecies, T1> && std::is_base_of_v<ChargedSpecies, T2> )
        Omega = new CoulombOmega(new T1, new T2);
    else
        Omega = new NonCoulombOmega(this);
}

template<typename T1 ,typename T2 >
void CollisionIntegral<T1, T2>::Pot ( Potential* pott ) {
    
    Omega = new NonCoulombOmega( this, this, pott );

}

template<typename T1, typename T2>
void CollisionIntegral<T1, T2>::MultiPot(std::initializer_list<Potential*> potentials,
    std::initializer_list<double> gss) {

    if (potentials.size() != gss.size()) {

        throw std::runtime_error(
        
            "Lacking potentials in MultiPot setting \n"
            "  potentials and degeneracies"
            "must be the same number."
        
        );

    } else {

        std::vector<CsCalculator*> chis ; 
        std::vector<double> gs ; 
        
        for (auto chi : potentials ) {
            if (auto morsePot = dynamic_cast<Morse*>(chi)) 
                chis.push_back( new AvrgChiIntegrator ( this, chi ) ) ;
            else
                chis.push_back( new AdaptChiIntegrator ( this, chi ) ) ;        
        }
  
        for (auto g : gss)  
            gs.push_back(g) ; 
        
        TCScalculator = new MultiCs(this->GetIntInterface(),chis,gs) ; 
        
    }
}

template<typename T1 ,typename T2 >
void CollisionIntegral<T1, T2>::ChargeTransfer( double A, double B ) {

    Omega = new ChargeExchangeOmega( this, Omega, A, B ) ;

} 

template<typename T1 ,typename T2 >
void CollisionIntegral<T1, T2>::ComputeCollisionIntegral(double Te, double Th, double lambda) {
    
    try {

        double Ti, Tj;
        
        if constexpr (std::is_base_of_v<Electron, T1>)
            Ti = Te;
        else
            Ti = Th;

        if constexpr (std::is_base_of_v<Electron, T2>)
            Tj = Te;
        else
            Tj = Th;

        double mi = this->getSp1()->getMass();
        double mj = this->getSp2()->getMass();
        
        TijStar = pow((1 / (mi + mj)) * ((mi / Tj) + (mj / Ti)), -1);

        omega4th = {
            
            Omega->Compute(1, 1, TijStar, lambda, this->TCScalculator), 
            Omega->Compute(1, 2, TijStar, lambda, this->TCScalculator),
            Omega->Compute(1, 3, TijStar, lambda, this->TCScalculator), 
            Omega->Compute(1, 4, TijStar, lambda, this->TCScalculator),
            Omega->Compute(1, 5, TijStar, lambda, this->TCScalculator), 
            Omega->Compute(1, 6, TijStar, lambda, this->TCScalculator),
            Omega->Compute(1, 7, TijStar, lambda, this->TCScalculator), 
            
            Omega->Compute(2, 2, TijStar, lambda, this->TCScalculator), 
            Omega->Compute(2, 3, TijStar, lambda, this->TCScalculator), 
            Omega->Compute(2, 4, TijStar, lambda, this->TCScalculator), 
            Omega->Compute(2, 5, TijStar, lambda, this->TCScalculator), 
            Omega->Compute(2, 6, TijStar, lambda, this->TCScalculator),
            
            Omega->Compute(3, 3, TijStar, lambda, this->TCScalculator), 
            Omega->Compute(3, 4, TijStar, lambda, this->TCScalculator), 
            Omega->Compute(3, 5, TijStar, lambda, this->TCScalculator), 
            
            Omega->Compute(4, 4, TijStar, lambda, this->TCScalculator)
        };
    
    } catch (const std::exception& e) {
        throw;
    }
}


template<typename T1, typename T2>
std::string CollisionIntegral<T1,T2>::InteractionName(){

    return this->getSp1()->getFormula() + "_" + this->getSp2()->getFormula() ;

}


template<typename T1, typename T2>
void CollisionIntegral<T1,T2>::info() {
    
    double w = 12.;
    
    std::string loadd = "none";

    std::string interaction_name = this->InteractionName();
    
    std::regex pattern("([^_]+)_([^_]+)");
    std::smatch match;
    
    std::string stringa1, stringa2;

    if (std::regex_search(interaction_name, match, pattern) && match.size() >= 3) {
        stringa1 = match[1].str();
        stringa2 = match[2].str();
    } else {
  
        stringa1 = interaction_name;
        stringa2 = "";
    }
  
    std::cout << std::left;
    std::cout << std::setw(6) << stringa1  << "-" <<  std::right << std::setw(6) 
        << stringa2 << std::setw(2+((3./2.)*w)) ; Omega->info() ; 
}

#endif