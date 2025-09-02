 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef TCS_CALCULATOR_H
#define TCS_CALCULATOR_H 

#include"DataLoader.h"
#include<string>
#include<vector>
#include <stdexcept>

// forward declarations
class InteractionInterface ; 
class Potential ; 

/// @brief Abstract base class for transport cross section calculators. 
class CsCalculator {

    public:
    
    /// @brief Interaction label used to build file names (e.g. "Ar_e-").
    std::string Name;

    /// @brief Energy grid [eV] for cross section data.
    std::vector<double> E;

    /// @brief Matrix of Qᶩ(E), l = 1..4, in [Å²]. Size: E.size() x 4.
    std::vector<std::vector<double>> Q;

    /**
     * @brief Returns the vector Q^l(E) for a given order l = 1..4.
     * @param l Order of the momentum transfer (1–4).
     * @return Vector of cross sections for given l.
     * @throws std::invalid_argument if l is out of range.  */
    std::vector<double> operator()(int l);

    bool computed = false ; 

    /// @brief Computes the cross sections. Must be implemented in derived classes.
    virtual void Compute() = 0;

    protected:

    /**
     * @brief Constructs a CsCalculator with interaction metadata.
     * @param i Pointer to the interaction. */
    CsCalculator(InteractionInterface* i);

    CsCalculator() = default;
    virtual ~CsCalculator() = default;
};

/// @brief Combines elastic and inelastic cross section calculators.
class CsHolder : public CsCalculator {

    public:

    /// @brief Pointer to elastic cross section calculator.
    CsCalculator* Qe;

    /// @brief Pointer to inelastic cross section calculator.
    CsCalculator* Qin;

    /// @brief Constructs from two calculators.
    CsHolder(CsCalculator* Qe, CsCalculator* Qin);

    /// @brief Constructs from interaction interface (placeholder).
    CsHolder(InteractionInterface* i);

    /// @brief Default constructor.
    CsHolder() = default;

    ~CsHolder() = default;

    /// @brief Calls Compute() on both Qe and Qin.
    void Compute() override;

};

/// @brief Container for multiple cross section calculators (e.g., multi-channel).
class MultiCs : public CsCalculator {

    public:
    
    /// @brief Vector for multiple calculators.
    std::vector<CsCalculator*> Qs;

    /// @brief State degeneracies "g" vector 
    std::vector<double> statesG ; 

    /// @brief Construct MultiCs object, calculators must be passed by the user.
    MultiCs ( InteractionInterface* i, std::vector<CsCalculator*> c, std::vector<double> gs );

    /// @brief Access calculator by index.
    CsCalculator*& operator[](int i);

/*     /// @brief Access Q^l of s-th calculator.
    std::vector<double> operator()(int s, int l);
 */
    /// @brief Returns number of calculators.
    int Size();

    /// @brief Computes all sub-calculators.
    void Compute() override;

    protected:

    /**
     * @brief USE ONLY on logics of multpiple cross sections 
     * without state degeneracies
     * @example ThresholdCs */
    MultiCs ( InteractionInterface* i, std::vector<CsCalculator*> c) ; 

    ~MultiCs() = default;

};

class ThresholdCs : public MultiCs {

    protected:

    /// @brief Limiting energies for Elim.size()+1 calculators
    std::vector<double> Elim ; 

    ~ThresholdCs(){}

    public:

    /// @brief Construct a new object, n calculators and n-1 limiting energies must be passed by user
    ThresholdCs( InteractionInterface* i, std::vector<CsCalculator*> c, std::vector<double> Elim ) ;

    /// @brief Computes sub-calculators and builds correct CsCalculator::E and CsCalculator::Q
    virtual void Compute() override ;

};

/// @brief Abstract class to implement default CSV parser for TCS_ datafiles.
class TcsParser {

    protected:
    /**
     * @brief Parses a default TCS_ CSV file.
     * @param file Input file stream (already opened).
     * @param E Output vector of energy values [eV].
     * @param Q Output matrix of cross sections (4 identical columns for each energy).
     * @details Assumes the CSV contains two columns:  
     * - First column: energy in eV  
     * - Second column: transport cross section Qᵐ(E) in Å²  
     * Each row of Q is initialized with four identical values Q^l(E) = Qᵐ(E) for l = 1..4.
     * This is a placeholder for generalizing to multipolar expansions.
     * @throws std::runtime_error if data is malformed or inconsistent.  */
    void TcsParse(std::ifstream& file, std::vector<double>& E, std::vector<std::vector<double>>& Q);

};
/// @brief Loads elastic transport cross section data from CSV files.
/// @see class CsCalculator for interface to cross section computation.
/// @see class DataLoader for file handling.
/// @see class TcsParser for CSV parsing.
class ElasticLoader : public CsCalculator, public DataLoader, public TcsParser {

    public:

    /**
     * @brief Constructs an ElasticLoader with default path and interaction pointer.
     * @param i Pointer to the interaction.  */
    ElasticLoader(InteractionInterface* i);

    /// @brief Constructs an ElasticLoader with a custom prefix and interaction.
    ElasticLoader(const std::string& prefix, InteractionInterface* interaction);

    /// @brief Default constructor.
    ElasticLoader(){};

    /// @brief Destructor.
    ~ElasticLoader(){};

    /// @brief Loads data if not already loaded (calls Init).
    void Compute() override;

    protected:
    
    /// @brief Initializes the loader by triggering data load.
    void Init() override;

    /// @brief Builds the expected file name from the interaction name.
    std::string BuildFileName(const std::string& name) override;

    /**
     * @brief Parses the loaded CSV file.
     * @param file Stream of the file to be parsed.
     * @details Calls TcsParse from the base class.
     * @throws std::runtime_error if parsing fails.  */
    void ParseFile(std::ifstream& file) override;
};

/// @brief Loads inelastic transport cross section data from CSV files.
/// @see class CsCalculator for interface to cross section computation.
/// @see class DataLoader for file handling.
/// @see class TcsParser for CSV parsing.
class InelasticLoader : public CsCalculator, public DataLoader, public TcsParser {

    public:
    
    /**
     * @brief Constructs an InelasticLoader with default path and interaction pointer.
     * @param i Pointer to the interaction.  */
    InelasticLoader(InteractionInterface* i);

    /// @brief Constructs an InelasticLoader with a custom prefix and interaction.
    InelasticLoader(const std::string& prefix, InteractionInterface* interaction);

    /// @brief Default constructor.
    InelasticLoader(){};

    /// @brief Destructor.
    ~InelasticLoader(){};

    /// @brief Loads data if not already loaded (calls Init).
    void Compute() override;

    protected:
    
    /// @brief Initializes the loader by triggering data load.
    void Init() override;

    /// @brief Builds the expected file name from the interaction name.
    std::string BuildFileName(const std::string& name) override;

    /**
     * @brief Parses the loaded CSV file.
     * @param file Stream of the file to be parsed.
     * @details Calls TcsParse from the base class.
     * @throws std::runtime_error if parsing fails. */
    void ParseFile(std::ifstream& file) override;

};
/// @brief Loads phase shift data and computes transport cross sections from it.
/// @see class CsCalculator for base interface.
/// @see class DataLoader for file parsing and loading.
class PhaseShiftsLoader : public CsCalculator, public DataLoader {

    /// @brief Matrix of phase shifts ηₗ(E) read from file.
    std::vector<std::vector<double>> etaL;

    /// @brief Pointer to the associated interaction interface.
    InteractionInterface* inT;

    public:
    /**
     * @brief Constructs the loader and initializes from file.
     * @param i Pointer to the interaction.  */
    PhaseShiftsLoader(InteractionInterface* i);

    /// @brief Constructs the loader with custom prefix and interaction pointer.
    PhaseShiftsLoader(const std::string& prefix, InteractionInterface* i);

    /// @brief Default constructor.
    PhaseShiftsLoader(){};

    /// @brief Destructor.
    ~PhaseShiftsLoader(){};

    /**
     * @brief Triggers the loading and computation of cross sections.
     * @details Loads the phase shifts and computes multipole TCS integrals
     * Q^(l) for l = 1, 2, 3, 4 using quantum-mechanical formulae. */
    void Compute() override;

    protected:
    
    /// @brief Loads the CSV file if not already loaded.
    void Init() override;

    /**
     * @brief Computes transport cross sections from phase shifts.
     * @details Uses the quantum formula:
     * \f[
     * Q^{(l)}(E) = \sum_{\ell} \mathcal{F}_\ell^{(l)} \sin^2(\eta_\ell - \eta_{\ell'})
     * \f]
     * where l = 1, 2, 3, 4. Contributions are summed and scaled to Å².
     * Assumes `etaL` is padded to allow access to ηₗ₊₁...ηₗ₊₄. */
    void ComputeFromPhaseShifts();

    /// @brief Builds the filename using prefix and interaction name.
    std::string BuildFileName(const std::string& name);

    /**
     * @brief Parses the phase shift CSV file.
     * @param file Input file stream.
     * @details First column is E [eV], the rest are ηₗ(E). If < 4 columns are present,
     * values are extended by repeating the last column. Missing rows are skipped.
     * @throws std::runtime_error if file is malformed or inconsistent. */
    void ParseFile(std::ifstream& file) override;

};

/// @brief Analytic charge-exchange cross section model based on a log-velocity fit.
/// @see class CsCalculator for interface to energy-dependent cross section computation.
class ChargeTransferCs : public CsCalculator {

    /// @brief Fit parameter A in the analytic cross section expression.
    double A;

    /// @brief Fit parameter B in the analytic cross section expression.
    double B;

    /// @brief Reduced mass of the system [amu].
    double mu;

    public:
    /**
     * @brief Constructs the charge-exchange cross section object.
     * @param i Pointer to the interaction.
     * @param A Coefficient A in σ = ½ (A - B log(g))².
     * @param B Coefficient B in σ = ½ (A - B log(g))².
     * @details Initializes the energy grid (log-spaced from 1.15 meV to 433 eV),
     * computes reduced mass in atomic mass units. */
    ChargeTransferCs(InteractionInterface* i, double A, double B);

    /**
     * @brief Computes the cross section on the energy grid.
     * @details The cross section is computed as:
     * \f[
     * \sigma(E) = \frac{1}{2} \left( A - B \log g_{ij} \right)^2
     * \f]
     * where g₍ᵢⱼ₎ is the relative velocity at energy E.
     * Units: σ in Å² if A, B are consistent with log10(cm/s). */
    void Compute() override;
};

/** @brief Interface functions for Ab Initio classes to perform integration of deflection angles and transport 
 * cross sections from the model interaction potential for the interacting particles */
class AbInitioTcsIntegration : public CsCalculator {

    protected:

    /// @brief Interaction Potential
    Potential* pot ; 
    
    /// @brief Initialize energy range, override it to customize
    virtual void InitE() { E = logspace ( log10(1.15e-3), log10(433), 50 ); }

    /// @brief Set the model potential passed by the user
    void SetPot( Potential* pot ){ pot = pot ; }  
    
    /// @brief Get the stored model potential 
    Potential* GetPot() { return pot; }

    /// @brief Construct a new object
    AbInitioTcsIntegration ( InteractionInterface* i, Potential* Pot ) : CsCalculator(i) { pot = Pot; InitE(); }

    /// @brief Computes the deflection angle chi from model interaction potential
    virtual double deflectionAngle ( double b , double E ) = 0 ;

};

/** @brief Computes the deflection angle integrating the interaction potential 
 * through the algorithm presented in
 * @see G. Colonna, A. Laricchiuta, 
 * "General numerical algorithm for classical collision integral calculation", 
 * Comput. Phys. Commun. 178 (2008) 809–816, DOI: 10.1016/j.cpc.2008.01.039. */
class AdaptDeflAngle : public AbInitioTcsIntegration {

    /** @brief Formula 1 of the reference DOI: 10.1016/j.cpc.2008.01.039.
     * @param r interparticle distance [Ang] 
     * @param b Impact parameter [Ang]
     * @param E Energy of collision [eV]
     * @return double Chi, DeflectionAngle */
    double IntegrandKi( double r, double b, double E ) ;

    /// @brief Adapt integration step. 
    std::tuple<double, bool> AdaptStep( double dx1, double F_bar_x1, double F_x1 ) ;
    
    /// @brief formula 1 integral extremes. 
    double AnalyzeIntegrandKiMin( double b, double E ) ;
    
    /// @brief formula 1 integral extremes.
    double AnalyzeIntegrandKiMax( double b )  ;
    
    /// @brief Fractal adaptive integration of formula 1. 
    double FractalKi( double rmin,  double rmax, double b, double E ) ;

    protected:

    /// @brief pi - ki 
    double theta ; 

    /// @brief Construct a new DeflectionAngle object
    AdaptDeflAngle(InteractionInterface* i, Potential* pot) : AbInitioTcsIntegration(i,pot) {}

    /// @brief Computes the deflection angle chi from model interaction potential
    double deflectionAngle ( double b , double E ) override ;

};

/** @brief Computes the deflection angle integrating the interaction potential 
 * through the algorithm presented in
 * @see J.A.Barker, W.Fock and F. Smith 
 * "Calculation of Gas Transport Properties and The Interaction of Argon Atoms", 
 * Phys. Fluids 7, 897–903 (1964) , DOI: 10.1063/1.1711301 */
class AvrgDeflAngle : public AbInitioTcsIntegration {

    protected:

    /// @brief reference variables defined in the reference 
    double r0{1.}, e0{1.} ;

    /// @brief Adimensional variables of the algorithm defined in the reference 
    std::vector<double> ws, v, Gs ; 

    /// @brief  number of points for ws and impact parameter
    int N {2000} ;

    /// @brief Initialize the energy range 
    virtual void InitE() override { E = logspace ( log10(1.1501404511032503e-3), log10(433), 150 ); }

    /// @brief Construct a new AvrgDefAngle object 
    AvrgDeflAngle( InteractionInterface* i, Potential* pot ) ;
    
    /// @brief Computes the deflection angle from model interaction potential 
    double deflectionAngle ( double Bsi, double Gsi ) override ;

};

/** @brief Compute Elastic Cross Section integrating the DeflectionAngle 
 ** through the algorithm presented in
 ** @see G. Colonna, A. Laricchiuta, 
 ** "General numerical algorithm for classical collision integral calculation", 
 ** Comput. Phys. Commun. 178 (2008) 809–816, DOI: 10.1016/j.cpc.2008.01.039. */
class AdaptChiIntegrator : public AdaptDeflAngle {

    /// @brief Integrand for formula 2 of the reference DOI: 10.1016/j.cpc.2008.01.039. 
    std::vector<double> IntegrandQ ( double b, double E ) ; 
    
    /// @brief Fractal integration of formula 2 in DOI: 10.1016/j.cpc.2008.01.039. 
    std::vector<double> FractalQ( double bmin, double bmax, double Ei ) ;

    public:

    /// @brief Construct a new AdaptChiIntegrator object by calling its base constructor
    AdaptChiIntegrator( InteractionInterface* i , Potential* pot ) : AdaptDeflAngle(i,pot) {} 

    /** @brief Compute and populates Q, this method represent the start for the threefold 
     * integration to collision integral with the algorithm described in the reference to 
     * this class hierarchy */
    void Compute() override ; 

};

/** @brief Compute Elastic Cross Section integrating the DeflectionAngle 
 ** through the algorithm presented in
* @see J.A.Barker, W.Fock and F. Smith 
 * "Calculation of Gas Transport Properties and The Interaction of Argon Atoms", 
 * Phys. Fluids 7, 897–903 (1964) , DOI: 10.1063/1.1711301 */
class AvrgChiIntegrator : public AvrgDeflAngle {

    /// @brief Numerical determination of asymptotic reduced impact parameter 
    std::tuple < double, double > AsymptoticBmax ( double Bmax0 , double Gsi ) ;

    /// @brief Integrate the deflection angle to transport cross sections  
    std::vector<double> IntegrateChi ( std::tuple<double, double> B , double Gsi ) ;

    public:

    /// @brief Construct a new AvrgChiIntegrator object 
    AvrgChiIntegrator( InteractionInterface* i, Potential* pot ) : AvrgDeflAngle(i,pot) {}

    /** @brief Compute and populates Q, this method represent the start for the threefold 
     * integration to collision integral with the algorithm described in the reference to 
     * this class hierarchy */
    void Compute() override ; 

};

/// @brief Loads and integrates differential cross section (DCS) data.
/// @see class CsHolder for composite structure of Qe and Qin calculators.
/// @see class DataLoader for file management and loading.
class DcsLoader : public CsHolder, public DataLoader {

    /// @brief Scattering angles for elastic DCS [deg].
    std::vector<double> anglesElastic;

    /// @brief σ(θ, E) matrix for elastic scattering.
    std::vector<std::vector<double>> sigmaElastic;

    /// @brief Scattering angles for inelastic DCS [deg].
    std::vector<double> anglesInlastic;

    /// @brief σ(θ, E) matrix for inelastic scattering.
    std::vector<std::vector<double>> sigmaInelastic;

    /**
     * @brief Initializes data loading from file or URL.
     * @details Called before computation. Loads CSV if not already loaded. */
    void Init() override;

    /**
     * @brief Integrates DCS to compute multipolar transport cross sections Q^(l).
     * @details Uses trapezoidal integration over θ for each Qᵉ and Qⁱ, where:
     * \f[
     * Q^{(l)}(E) = \int_0^\pi \sigma(\theta, E) (1 - \cos \theta)^{l} \sin \theta d\theta
     * \f] */
    void IntegrateDifferentialCrossSections();

    /**
     * @brief Downloads DCS data from a URL and saves it locally.
     * @param hyperref Direct link to raw CSV data.
     * @throws std::runtime_error if CURL fails or download is empty. */
    void DownloadAndSaveData(const std::string& hyperref);

    /**
     * @brief Interpolates or removes invalid σ values marked by sentinel (9999999999).
     * @param energies Vector of energy values.
     * @param angles Vector of angles.
     * @param sigma 2D matrix of DCS values.
     * @param invalidPositions List of invalid indices. */
    void FixInvalidValues ( std::vector<double>& energies,
        std::vector<double>& angles,
            std::vector<std::vector<double>>& sigma,
                const std::vector<std::pair<int, int>>& invalidPositions );

    public:
    
    /// @brief Constructs the loader and initializes from default DCS file.
    DcsLoader(InteractionInterface* interaction);

    /// @brief Constructs the loader with a custom prefix and file name.
    DcsLoader(const std::string& prefix, InteractionInterface* interaction);

    /// @brief Constructs the loader by downloading data from URL.
    DcsLoader(InteractionInterface* interaction, const std::string& hyperref);

    /**
     * @brief Computes integrated transport cross sections.
     * @details Loads and parses DCS data if needed, then computes Qᵉ and Qⁱ. */
    void Compute() override;

    protected:
    
    /// @brief Builds expected DCS file name.
    std::string BuildFileName(const std::string& name) override;

    /**
     * @brief Parses the DCS CSV format.
     * @param file Input stream of the CSV file.
     * @details Recognizes "Elastic" / "Inelastic" sections, extracts energy/angle/sigma,
     * handles missing values with placeholder and post-processing. */
    void ParseFile(std::ifstream& file) override;
    
};

#endif
