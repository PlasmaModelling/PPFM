 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

/* 
COMPUTATION OF TRANSPORT CROSS SECTION (Ang^2) FUNCTION OF E(eV) 
reference: 
    G. Colonna, A. Laricchiuta, 
    "General numerical algorithm for classical collision integral calculation", 
    Comput. Phys. Commun. 178 (2008) 809–816, DOI: 10.1016/j.cpc.2008.01.039.
*/

#ifndef TRANSPORT_CROSS_SECTION_H
#define TRANSPORT_CROSS_SECTION_H

#include"Interaction.h"
#include"TcsCalculator.h"
#include"GasMixture.h"


/// @brief Interface class for Transport Cross Section (TCS) objects.
/// @see CsCalculator for cross section evaluation.
class TcsInterface {

    public:
    
    /** @brief Pointer to the base class for cross section calculators.
     * @details Derived calculator objects (Elastic, Inelastic, DCS) can be assigned polymorphically.
     * This pointer is used internally for computing transport cross sections.
     * @see CsCalculator for details. */
    CsCalculator* TCScalculator ;

    /// @brief Loads elastic cross section data.
    virtual void LoadElastic() = 0;

    /// @brief Loads elastic cross section with a custom prefix.
    virtual void LoadElastic(const std::string& customPrefix) = 0;

    /// @brief Loads inelastic cross section data.
    virtual void LoadInelastic() = 0;

    /// @brief Loads inelastic cross section with a custom prefix.
    virtual void LoadInelastic(const std::string& customPrefix) = 0;

    /// @brief Loads differential cross section (DCS) data.
    virtual void LoadDCS() = 0;

    /// @brief Loads DCS with a custom prefix.
    virtual void LoadDCS(const std::string& customPrefix) = 0;

    /// @brief Downloads DCS from a hyperlinked source.
    virtual void DownloadDCS(const std::string& hyperref) = 0;

    /// @brief Returns the base interface to the current interaction.
    virtual InteractionInterface* GetIntInterface() = 0;

    /**
     * @brief Get the interface by cloning the current object polymorphically.
     * @details Realization in concrete class.
     * @return Pointer to a deep copy of the derived `TcsInterface` instance. */
    virtual TcsInterface* GetTcsInterface() = 0;

};

/**
 * @brief Template class for Transport Cross Section objects.
 * @tparam T1 Type of the first interacting species.
 * @tparam T2 Type of the second interacting species.
 * @see Interaction for basic species interaction structure. */
template <typename T1, typename T2>
class TransportCrossSection : public Interaction<T1, T2>, public virtual TcsInterface {

    protected:
    /**
     * @brief Constructs a TransportCrossSection for the given species pair.
     * @param t1 Pointer to the first species.
     * @param t2 Pointer to the second species.
     * @details Inherits base interaction data and enables cross section computation via CsCalculator. */
    TransportCrossSection(T1* t1, T2* t2);

    /// @brief Loads elastic cross section using the default path.
    void LoadElastic() override;

    /// @brief Loads elastic cross section with a given prefix.
    void LoadElastic(const std::string& customPrefix) override;

    /// @brief Loads inelastic cross section using default settings.
    void LoadInelastic() override;

    /// @brief Loads inelastic cross section with a given prefix.
    void LoadInelastic(const std::string& customPrefix) override;

    /// @brief Loads differential cross section data.
    void LoadDCS() override;

    /// @brief Loads DCS with a given prefix.
    void LoadDCS(const std::string& customPrefix) override;

    /// @brief Downloads DCS from an external reference.
    void DownloadDCS(const std::string& hyperref) override;

    
    public:

    /// @brief Returns the current object as a base interaction interface.
    InteractionInterface* GetIntInterface() override { return this; }
    
    /** @brief Returns a polymorphic clone of the current object.
     ** @return Deep copy of this TransportCrossSection. */
    TcsInterface* GetTcsInterface() override;
    

};

//__________________________________ Implementation ___________________________________

template <typename T1, typename T2>
void TransportCrossSection<T1, T2>::LoadElastic() {

    this->TCScalculator = new ElasticLoader( this ) ;

}; 

template <typename T1, typename T2>
void TransportCrossSection<T1, T2>::LoadElastic( const std::string& customPrefix ) {

    this->TCScalculator = new ElasticLoader( customPrefix , this );

}

template <typename T1, typename T2>
void TransportCrossSection<T1, T2>::LoadInelastic() {

    this->TCScalculator = new CsHolder(TCScalculator, new InelasticLoader(this));

}

template <typename T1, typename T2>
void TransportCrossSection<T1, T2>::LoadInelastic( const std::string& customPrefix ) {

    this->TCScalculator = new CsHolder(TCScalculator, 
        new InelasticLoader(customPrefix,this));

}

template <typename T1, typename T2>
void TransportCrossSection<T1, T2>::LoadDCS() {

    this->TCScalculator = new DcsLoader( this ) ;

}; 

template <typename T1, typename T2>
void TransportCrossSection<T1, T2>::LoadDCS( const std::string& customPrefix ) {

    this->TCScalculator = new DcsLoader( customPrefix, this ) ;

}; 

template <typename T1, typename T2>
void TransportCrossSection<T1, T2>::DownloadDCS( const std::string& hyperref ) {

    this->TCScalculator = new DcsLoader( this , hyperref ) ;

}; 

template <typename T1, typename T2>
TransportCrossSection<T1, T2>::TransportCrossSection(T1* t1, T2* t2) : Interaction<T1, T2>(t1, t2) {
    TCScalculator = nullptr ; 
}

template <typename T1, typename T2>
TcsInterface* TransportCrossSection<T1, T2>::GetTcsInterface() { return new TransportCrossSection<T1, T2>(*this); }

#endif
