/* 
COMPUTATION OF COLLISION INTEGRALS FUNCTION OF T(K°)    
reference: 
    R. S. Devoto, "Transport Properties of Ionized Monatomic Gases",
    Physics of Fluids, vol. 9, n. 6, pp. 1230–1240, June 1966, DOI: 10.1063/1.1761885.
    Formula 11.
*/

#ifndef COLLISION_INTEGRAL_H
#define COLLISION_INTEGRAL_H


#include"TransportCrossSection.h"
#include"OmegaCalculator.h"
#include"DataPrinter.h"

#include <regex>

/// @brief Interface class for CollisionIntegral objects. 
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
     * so omega4th will be 16 elements long, those needed for the approximation */
    std::vector<double> omega4th;

    /// @brief Effective temperature of collision
    double TijStar;

    public: 

    OmegaCalculator* Omega;
    
    /** @brief Interface function to compute the Collision Integral
     * @param temperature Temperature of collision [K°]
     * @param lambda Debye length [m] */
    virtual void ComputeCollisionIntegral ( double Te, double Th, double lambda ) = 0 ;
    
    /** @brief Interface function to write the name of the interaction
     * @return std::string  */
    virtual std::string InteractionName() = 0 ;
    
    /** @brief Interface function to set the Interaction Potential 
     * @param pott Potential object
     * @see class Potential and derived for further details on implemented Potentials */
    virtual void Pot ( Potential* pott ) = 0;
    
    /** @brief Interface function to set the Interaction Potentials for multiple interactions 
     * @param potentials Initializer_list of { Potential objects }
     * @param gs Initializer_list of { interaction statistical degeneracies }
     * @see class Potential and derived for further details on implemented Potentials */
    virtual void MultiPot ( std::initializer_list<Potential*> potentials ,
        std::initializer_list<double> gs ) = 0 ;
    
    /** @brief Interface function to set the degenereacies for multi state interactions 
     * @param gs Initializer_list of { interaction statistical degeneracies } */
/*     virtual void MultiState ( std::initializer_list<double> gs ) = 0 ;*/

    /** @brief Interface function to set parameters for the Charge Exchange collision cross 
     ** section Q(g)=(A-Bln(g))^2 
     ** @param A [Ang]
     ** @param B [Ang] */
    virtual void ChargeTransfer ( double A, double B ) = 0 ; 
    
    /// @brief Info on the collision integrals of the Mixture 
    virtual void info() = 0 ;
    
    /// @brief Interface function to call to load the collision integral from its raw file
    virtual void Load( bool b ) = 0 ;

} ;

// CLass for the user to operate on both Transport Cross Sections and Collision Integrals
class HybridInterface : public virtual CInterface, public virtual TcsInterface {};

template<typename T1, typename T2>
class CollisionIntegral : public TransportCrossSection<T1,T2>, public HybridInterface {
        
    /// @brief initialize the object
    void InitOmega4th() { omega4th.resize(16) ; }

    /// @brief Initialize the right OmegaCalculator object, @see OmegaCalculator.h for more 
    void InitCalculator() ;

    /// @brief Default init with no loader overload
    void InitCalculator(bool b) ;

    public:

    CollisionIntegral(T1* t1,T2* t2) : TransportCrossSection<T1,T2>(t1,t2) {InitOmega4th(); InitCalculator();}

    /** @brief Compute collision integral and store it in the omega4th member of this class.
     * this will be the preferred method if the file is available. \n
     * @param temperatura  [K°]
     * @param lambda Debye Length [m] */
    void ComputeCollisionIntegral(double Te, double Th, double lambda) override ;

    /// @brief function returning "Formula1_Formula2"
    std::string InteractionName() override ;
        
    /** @brief Interface function to set the Interaction Potential 
     * @param pott Potential object
     * @see class Potential and derived for further details on implemented Potentials */
    void Pot ( Potential* pott ) ;
    
    /** @brief Interface function to set the Interaction Potentials for multiple interactions 
     * @details Initialize a MultiPotOmega and a MultiCs new objects to members Omega and TCScalculator 
     * @param potentials Initializer_list of { Potential objects }
     * @param gs Initializer_list of { interaction statistical degeneracies }
     * @see class Potential and derived for further details on implemented Potentials */
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

//________________________________ Printing class ________________________________

/// @brief Prints collision integrals to CSV files.
/// @see class DataPrinter for CSV output interface.
class CollisionIntegralCsv : public DataPrinter {

    friend class CiBox;

    /// @brief Pointer to collision integral interface.
    CInterface* ci;

    /// @brief Builds the output filename with "CI_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header for collision integrals.
    void PrepareHeader() override;

    /**
     * @brief Computes and stores collision integrals at different temperatures.
     * @param x Vector of reduced temperatures Tij*.
     * @param gasmix Pointer to the gas mixture.
     * @details For each temperature, computes the 4th order set of Ω(l,s) integrals,
     * using Debye length and θ to evaluate each term. Results are stored for CSV output.
     * @see class CInterface for ComputeCollisionIntegral. */
    void PrepareData(const std::vector<double>& x, GasMixture* gasmix) override;

    /**
     * @brief Prints a message confirming successful output.
     * @param filename Name of the written file. */
    void PrintMessage ( const std::string& filename ) override;

    /**
     * @brief Overrides default print to prepend CollisionIntegrals_ subfolder.
     * @param filename Base filename.
     * @param x Vector of reduced temperatures.
     * @param gasmix Pointer to the gas mixture.
     * @details Temporarily modifies the output folder to include a dedicated
     * subfolder for collision integrals, then calls the base print method.
     * @see class DataPrinter for Print method. */
    void Print ( const std::string& filename, const std::vector<double>& x, GasMixture* gasmix ) override;

    public:
    
    /// @brief Constructor from a HybridInterface pointer.
    CollisionIntegralCsv ( CInterface* _ci );

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

/* template<typename T1, typename T2>
void CollisionIntegral<T1, T2>::MultiState ( std::initializer_list<double> gss ) {
    
    std::vector<CsCalculator*> nullcalcs (gss.size(),nullptr) ; 
    Omega = new MultiPotOmega(this,this,gss) ;
    TCScalculator = new MultiCs(this,nullcalcs) ;
    
} */

template<typename T1 ,typename T2 >
void CollisionIntegral<T1, T2>::ChargeTransfer( double A, double B ) {

    Omega = new ChargeExchangeOmega( this, Omega, A, B ) ;

} 

template<typename T1 ,typename T2 >
void CollisionIntegral<T1, T2>::ComputeCollisionIntegral(double Te, double Th, double lambda) {
    
    try {

        double Ti, Tj;
        
        Ti = dynamic_cast<Electron*>(this->getSp1()) ? Te : Th;
        Tj = dynamic_cast<Electron*>(this->getSp2()) ? Te : Th;

        double mi = this->getSp1()->getMass();
        double mj = this->getSp2()->getMass();
        
        // Effective temperature of collision
        TijStar = pow((1 / (mi + mj)) * ((mi / Tj) + (mj / Ti)), -1);

        // direct calculation through the OmegaCalculator object. 
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