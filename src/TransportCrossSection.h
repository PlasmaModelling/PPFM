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
     * @brief Clones the current object polymorphically.
     * @return Pointer to a deep copy of the derived `TcsInterface` instance. */
    virtual TcsInterface* clone() = 0;

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

    /// @brief Returns the current object as a base interaction interface.
    InteractionInterface* GetIntInterface() override { return this; }

    /**
     * @brief Returns a polymorphic clone of the current object.
     * @return Deep copy of this TransportCrossSection. */
    TcsInterface* clone() override;

};


/// @brief Prints collision integrals to CSV files.
/// @see class DataPrinter for CSV output interface.
class TransportCrossSectionCsv : public DataPrinter {

    friend class CiBox;

    /// @brief Pointer to transport cross section interface.
    TcsInterface* tcs ;

    /// @brief Builds the output filename with "TCS_" prefix.
    std::string BuildFileName(const std::string& filename) const override;

    /// @brief Prepares the CSV header for elastic transport cross sections.
    void PrepareHeader() override;

    /// @brief Prepares the CSV header for inelastic transport cross sections.
    void PrepareInelasticHeader() ; 

    void PrepareMultiHeader( MultiCs* cscalc , int i ) ; 

    /**
     * @brief Computes and stores elastic transport cross sections.
     * @param x Not-used.
     * @param gasmix Not used.
     * @details Transport Cross Section computation is static within 
     * energy ranges in eV, no data should be passed as they're incapsulated
     * in the tcs->TCScalculator object.
     * @see classes in file TcsCalculator.h. */
    void PrepareData(const std::vector<double>& x, GasMixture* gasmix) override ;

    /// @brief Overload function to use in case of a MultiCs.
    /// @param cscalc 
    void PrepareData( MultiCs* cscalc , int i ) ;

    /**
     * @brief Computes and stores inelastic transport cross sections.
     * @param tcsElIn Holder to extract Qin.
     * @details Transport Cross Section computation is static within 
     * energy ranges in eV, no external data should be passed as they're incapsulated
     * in the tcs->TCScalculator object.
     * @see classes in file TcsCalculator.h. */
    void PrepareInelasticData ( CsHolder* tcsElIn ) ;

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
    TransportCrossSectionCsv ( TcsInterface* _ci );

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
TcsInterface* TransportCrossSection<T1, T2>::clone() { return new TransportCrossSection<T1, T2>(*this); }

#endif
