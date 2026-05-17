 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "DataPrinter.h"
#include "GasMixture.h"
#include "Composition.h"
#include "Thermodynamics.h"
#include "TcsCalculator.h"
#include "Devoto.h"
#include "ZhangMurphyTP.h"
#include "PartitionFunction.h"
#include "CollisionIntegral.h"

// Costruisce il path completo (di default: nella directory "out/")
std::string DataPrinter::BuildFilePath( const std::string& filename ) const {

    std::filesystem::path base = std::filesystem::path(PPFM_PROJECT_ROOT) / "out";

    if (!customFolder.empty())
        base /= customFolder;

    std::filesystem::create_directories ( base ); 

    return ( base / BuildFileName ( filename ) ).string() ;
}

void DataPrinter::WriteData(const std::string& filename ) {

    std::ofstream file(BuildFilePath(filename));
    
    if (!file.is_open()) {
    
        std::cerr << "Error: Unable to open file '" << filename << "'." << std::endl;
        return;
    
    }

    file << header << std::endl;

    int N = data.size();

    std::vector<std::string> lines( N );

    for (size_t i = 0; i < data.size(); ++i) {
    
        std::ostringstream oss;
        const auto& row = data[i];
    
        for (size_t j = 0; j < row.size(); ++j) {
            oss << row[j];
            if (j != row.size() - 1)
                oss << ",";
    
        }
    
        lines[i] = oss.str();
    
    }

    for (const auto& line : lines) {
    
        file << line << std::endl;
    
    }

    file.close();

}

void DataPrinter::AppendData(const std::string& filename) {

    std::ofstream file(BuildFilePath(filename), std::ios::app);

    if (!file.is_open()) {

        std::cerr << "Error: Unable to open file '" << filename << "' for appending." << std::endl;
        return;

    }

    // Scrive l'header prima dei nuovi dati
    file << header << std::endl;

    // Scrive i dati
    for (const auto& row : data) {
        for (size_t j = 0; j < row.size(); ++j) {
            file << row[j];
            if (j != row.size() - 1) {
                file << ",";
            }
        }
        file << std::endl;
    }

    file.close();
}

// Metodo finale che esegue la stampa su file
void DataPrinter::Print ( const std::string& filename, const std::vector<double>& x, GasMixture* gasmix ) {
    
    PrepareHeader();

    PrepareData(x, gasmix);
    
    if (data.empty() || header.empty()) {

        std::cerr << "Error: data or header not prepared." << std::endl;
        return;
    
    }

    WriteData ( filename ) ; 

    PrintMessage(filename) ; 

}

// __________________________ CompositionCsv __________________________ //

CompositionCsv::CompositionCsv(Composition* solver) : solver(solver) {}

CompositionCsv::CompositionCsv(Composition* solver, const std::string& folder)
    : solver(solver) {
        customFolder = folder;
}

std::string CompositionCsv::BuildFileName(const std::string& filename) const {
    return "Comp_" + filename + ".csv";
}

void CompositionCsv::PrepareHeader() {

    header = "Te [K],";
    for (int i = 0; i < solver->compositions().size(); i++) {
        header += "n_{" + (*solver->mixptr)(i)->getFormula() + "}[#/m^3],";
    }
    header += "n_{tot}[#/m^3]";
}

void CompositionCsv::PrepareData(const std::vector<double>& temperatureRange, GasMixture* mix) {

    double T0 = mix->getTemperature();
    double theta = mix->theta->get();

    std::vector<std::vector<double>> ns(temperatureRange.size(),
        std::vector<double>(mix->getN()));

    data.clear() ; 
    data.resize(temperatureRange.size(), std::vector<double>(mix->getN() + 2));

    solver = mix->getCompositionObj();

    for (int i = 0; i < temperatureRange.size(); i++) {
        
        mix->setT(temperatureRange[i]);

        ns[i] = solver->compositions();

        data[i][0] = temperatureRange[i] * theta;
        for (int j = 1; j < ns[0].size() + 1; j++)
            data[i][j] = ns[i][j - 1];
        
        data[i][ns[0].size() + 1] = solver->ntot();
    }

    mix->setT(T0);
    mix->restartComposition();
}

void CompositionCsv::PrintMessage(const std::string& filename) {
    std::cout << "Composition " << filename << " printed." << std::endl;
}


// ---------------- ThermodynamicsCsv ---------------- //

std::string ThermodynamicsCsv::BuildFileName(const std::string& filename) const {
    return "TH_" + filename + ".csv";  
}

void ThermodynamicsCsv::PrepareHeader() {
    header = "Te [K], ρ [kg/m³], Cₚ [J/(kg·K)], hₑ + hₕ [J/kg], γ [–], a [m/s]";
}
    
void ThermodynamicsCsv::PrintMessage(const std::string& filename) {
    std::cout << "Thermodynamic properties " << filename << " printed." << std::endl;
}

ThermodynamicsCsv::ThermodynamicsCsv(Thermodynamics* solver) : solver(solver) {}

ThermodynamicsCsv::ThermodynamicsCsv(Thermodynamics* solver, const std::string& folder) : solver(solver) {
    customFolder = folder;
}

void ThermodynamicsCsv::PrepareData(const std::vector<double>& temperatureRange, GasMixture* gasmix) {
        
    double T0 = gasmix->getTemperature();  
    double theta = gasmix->theta->get();

    data.resize(temperatureRange.size(), std::vector<double>(6));

    for (int i = 0; i < temperatureRange.size(); i++) {
        gasmix->setT(temperatureRange[i]);        
        solver->computeThermodynamics(*gasmix);
    
        data[i][0] = temperatureRange[i] * theta;
        data[i][1] = solver->rho();
        data[i][2] = solver->cp();
        data[i][3] = solver->he() + solver->hh();
        data[i][4] = solver->gamma();
        data[i][5] = solver->a();
    }

    gasmix->setT(T0); 
    gasmix->restartComposition(); 
    solver->computeThermodynamics(*gasmix);
} 

// __________________________ DevotoCsv __________________________ //


std::string DevotoTpCsv::BuildFileName(const std::string& filename) const {
    return "TP_" + filename + ".csv";
}

void DevotoTpCsv::PrepareHeader() {
    
    if (dynamic_cast<DevotoLteLambdaR*>(solver))
        header = "Th [K], λₑ [W/(m·K)], λₕ [W/(m·K)], λᵣ [W/(m·K)], λint [W/(m·K)], μ [Pa·s], σ [S/m], HeatExchange[J]";    
    else    
        header = "Th [K], λₑ [W/(m·K)], λₕ [W/(m·K)], μ [Pa·s], σ [S/m], HeatExchange[J]";

    }

void DevotoTpCsv::PrintMessage(const std::string& filename) {
    std::cout << "Devoto transport properties " << filename << " printed." << std::endl;
}

DevotoTpCsv::DevotoTpCsv(DevotoTP* solver) : solver(solver) {}

DevotoTpCsv::DevotoTpCsv(DevotoTP* solver, const std::string& folder)
    : solver(solver) {
    customFolder = folder;
}

void DevotoTpCsv::PrepareData(const std::vector<double>& temperatureRange, GasMixture* gasmix) {

    double T0 = gasmix->getTemperature();
    data.resize(temperatureRange.size(), std::vector<double>(solver->Tp.size()+1));

    for (size_t i = 0; i < temperatureRange.size(); ++i) {
        
        gasmix->setT(temperatureRange[i]);
        solver->computeTransport(gasmix);
        // solver->lightComputeTransport(gasmix);
        
        data[i] = concatenate({temperatureRange[i]}, solver->Tp) ;

    }

    gasmix->setT(T0);
    gasmix->restartComposition();

}

// __________________________ ZhangTpCsv __________________________ //


std::string ZhangTpCsv::BuildFileName(const std::string& filename) const {
    return "TP_" + filename + ".csv";
}

void ZhangTpCsv::PrepareHeader() {
    header = "Th [K], kₑ [W/(m·K)], kₕ [W/(m·K)], μ [Pa·s], σ [S/m], kₑθ [W/m], kₕθ [W/m]";
}

void ZhangTpCsv::PrintMessage(const std::string& filename) {
    std::cout << "Zhang–Murphy transport properties " << filename << " printed." << std::endl;
}

ZhangTpCsv::ZhangTpCsv(ZhangMurphyTP* solver) : solver(solver) {}

ZhangTpCsv::ZhangTpCsv(ZhangMurphyTP* solver, const std::string& folder)
    : solver(solver) {
    customFolder = folder;
}

void ZhangTpCsv::PrepareData(const std::vector<double>& temperatureRange, GasMixture* gasmix) {

    double T0 = gasmix->getTemperature();
    data.resize(temperatureRange.size(), std::vector<double>(7));

    for (size_t i = 0; i < temperatureRange.size(); ++i) {
        gasmix->setT(temperatureRange[i]);
        solver->computeTransport(gasmix);

        data[i][0] = temperatureRange[i];
        for (int j = 0; j < 6; ++j)
            data[i][j + 1] = solver->Tp[j];
    }

    gasmix->setT(T0);
    gasmix->restartComposition();
    solver->computeTransport(gasmix);
}
// ___________________________ DeflectionAngleCsv __________________________ //

std::string DeflectionAngleCsv::BuildFileName(const std::string& filename) const {

    return "CHI_" + filename + ".csv";

}


void DeflectionAngleCsv::PrepareHeader() {

    header = "1st column: b [A] 1st row: E [eV] matrix: chi(b,E) [rad]";

    if (thresholdSolver != nullptr) {

        header += " | Threshold energies [eV]:";

        for (size_t i = 0; i < thresholdSolver->E.size(); ++i)
            header += " " + std::to_string(thresholdSolver->E[i]);

    }

}


void DeflectionAngleCsv::PrintMessage(const std::string& filename) {

    std::cout << "Deflection angle " << filename << " printed." << std::endl;

}


bool DeflectionAngleCsv::IsPrintableThreshold(ThresholdCs* threshold) const {

    if (threshold == nullptr)
        return false;

    if (threshold->Qs.empty())
        return false;

    for (size_t i = 0; i < threshold->Qs.size(); ++i) {

        if (dynamic_cast<AbInitioTcsIntegration*>(threshold->Qs[i]) == nullptr)
            return false;

    }

    return true;

}


AbInitioTcsIntegration* DeflectionAngleCsv::SelectThresholdSolver(
    ThresholdCs* threshold,
    double E
) const {

    if (threshold == nullptr)
        return nullptr;

    if (threshold->Qs.empty())
        return nullptr;

    size_t index = 0;

    while (
        index < threshold->E.size() &&
        E > threshold->E[index]
    ) {
        ++index;
    }

    if (index >= threshold->Qs.size())
        return nullptr;

    return dynamic_cast<AbInitioTcsIntegration*>(threshold->Qs[index]);

}


void DeflectionAngleCsv::CollectPrintableSolvers(
    CsCalculator* calculator,
    const std::string& suffix
) {

    if (calculator == nullptr)
        return;


    if (auto chiSolver = dynamic_cast<AbInitioTcsIntegration*>(calculator)) {

        printableSolvers.push_back({chiSolver, nullptr, suffix});
        return;

    }


    // IMPORTANT: ThresholdCs must be checked before MultiCs,
    // because ThresholdCs derives from MultiCs.
    //
    // ThresholdCs is not printed branch-by-branch.
    // It is printed as a single composite solver using the proper branch
    // depending on the current energy.
    if (auto threshold = dynamic_cast<ThresholdCs*>(calculator)) {

        if (IsPrintableThreshold(threshold)) {

            printableSolvers.push_back({
                nullptr,
                threshold,
                suffix + "_threshold"
            });

        }

        return;

    }


    if (auto multi = dynamic_cast<MultiCs*>(calculator)) {

        for (size_t i = 0; i < multi->Qs.size(); ++i) {

            CollectPrintableSolvers(
                multi->Qs[i],
                suffix + "_state" + std::to_string(i + 1)
            );

        }

        return;

    }


    if (auto holder = dynamic_cast<CsHolder*>(calculator)) {

        CollectPrintableSolvers(
            holder->Qe,
            suffix + "_elastic"
        );

        return;

    }

}


DeflectionAngleCsv::DeflectionAngleCsv(TcsInterface* tcs) {

    impact_parameter = arange(0.1, 10.0, 0.2);
    energy = logspace(log10(1.15e-3), log10(433), 50);

    CollectPrintableSolvers(tcs->TCScalculator);

}


DeflectionAngleCsv::DeflectionAngleCsv(
    TcsInterface* tcs,
    const std::string& folder
) {

    impact_parameter = arange(0.1, 10.0, 0.2);
    energy = logspace(log10(1.15e-3), log10(433), 50);

    customFolder = folder;

    CollectPrintableSolvers(tcs->TCScalculator);

}


void DeflectionAngleCsv::PrepareData(const std::vector<double>&, GasMixture*) {

    data.clear();

    data.resize(
        impact_parameter.size() + 1,
        std::vector<double>(energy.size() + 1)
    );


    for (size_t i = 0; i < energy.size(); ++i)
        data[0][i + 1] = energy[i];


    for (size_t j = 0; j < impact_parameter.size(); ++j)
        data[j + 1][0] = impact_parameter[j];


    for (size_t i = 0; i < energy.size(); ++i) {

        AbInitioTcsIntegration* activeSolver = solver;

        if (thresholdSolver != nullptr)
            activeSolver = SelectThresholdSolver(thresholdSolver, energy[i]);

        if (activeSolver == nullptr)
            continue;


        for (size_t j = 0; j < impact_parameter.size(); ++j) {

            data[j + 1][i + 1] =
                activeSolver->deflectionAngle(
                    impact_parameter[j],
                    energy[i]
                );

        }

    }

}


bool DeflectionAngleCsv::ToSkip() const {

    return printableSolvers.empty();

}


void DeflectionAngleCsv::Print(
    const std::string& filename,
    const std::vector<double>& x,
    GasMixture* gasmix
) {

    if (printableSolvers.empty()) {

        std::cout << "Skipping " << filename
                  << ": no valid deflection-angle solver found."
                  << std::endl;

        return;

    }


    for (const auto& printable : printableSolvers) {

        solver = printable.solver;
        thresholdSolver = printable.thresholdSolver;

        if (solver == nullptr && thresholdSolver == nullptr)
            continue;

        std::string outputName = filename + printable.suffix;

        DataPrinter::Print(outputName, x, gasmix);

        solver = nullptr;
        thresholdSolver = nullptr;

    }

}
// ___________________________ TransportCrossSectionCsv __________________________ //

std::string TransportCrossSectionCsv::BuildFileName(const std::string& filename) const {

    return "TCS_" + filename + ".csv";

}


void TransportCrossSectionCsv::PrepareHeader() {

    if (currentInelastic)
        header = "E [eV], Qin(1) [A^2], Qin(2) [A^2], Qin(3) [A^2], Qin(4) [A^2]";
    else
        header = "E [eV], Qe(1) [A^2], Qe(2) [A^2], Qe(3) [A^2], Qe(4) [A^2]";


    if (currentThreshold != nullptr) {

        header += ", Threshold energies [eV]:";

        for (size_t i = 0; i < currentThreshold->E.size(); ++i)
            header += " " + std::to_string(currentThreshold->E[i]);

    }

}


void TransportCrossSectionCsv::PrintMessage(const std::string& filename) {

    std::cout << "Transport cross section " << filename << " printed." << std::endl;

}

void TransportCrossSectionCsv::PrepareData(const std::vector<double>&, GasMixture*) {

    data.clear();

    if (currentCalculator == nullptr)
        return;

    if (currentParentCalculator != nullptr)
        currentParentCalculator->Compute();
    else
        currentCalculator->Compute();

    auto& e = currentCalculator->E;
    auto& q = currentCalculator->Q;

    data.resize(e.size(), std::vector<double>(5));

    for (size_t i = 0; i < e.size(); ++i) {

        data[i][0] = e[i];

        for (size_t j = 0; j < 4; ++j)
            data[i][j + 1] = q[i][j];

    }

}

void TransportCrossSectionCsv::CollectPrintableTcs(
    CsCalculator* calculator,
    const std::string& suffix,
    bool inelastic
) {

    if (calculator == nullptr)
        return;


    if (auto threshold = dynamic_cast<ThresholdCs*>(calculator)) {

        printableTcs.push_back({
            threshold,
            threshold,
            nullptr,
            suffix + "_threshold",
            inelastic
        });

        return;

    }


    // IMPORTANT: DcsLoader must be checked before CsHolder,
    // because DcsLoader derives from CsHolder.
    //
    // The DcsLoader itself must be computed, because its Compute()
    // integrates the differential cross sections and fills Qe and Qin.
    if (auto dcs = dynamic_cast<DcsLoader*>(calculator)) {

        printableTcs.push_back({
            dcs->Qe,
            nullptr,
            dcs,
            suffix + "_elastic_DCS",
            false
        });

        printableTcs.push_back({
            dcs->Qin,
            nullptr,
            dcs,
            suffix + "_inelastic_DCS",
            true
        });

        return;

    }


    if (auto holder = dynamic_cast<CsHolder*>(calculator)) {

        CollectPrintableTcs(
            holder->Qe,
            suffix + "_elastic",
            false
        );

        CollectPrintableTcs(
            holder->Qin,
            suffix + "_inelastic",
            true
        );

        return;

    }


    if (auto multi = dynamic_cast<MultiCs*>(calculator)) {

        for (size_t i = 0; i < multi->Qs.size(); ++i) {

            CollectPrintableTcs(
                multi->Qs[i],
                suffix + "_state" + std::to_string(i + 1),
                inelastic
            );

        }

        return;

    }


    printableTcs.push_back({
        calculator,
        nullptr,
        nullptr,
        suffix,
        inelastic
    });

}

TransportCrossSectionCsv::TransportCrossSectionCsv(TcsInterface* solverr) {

    solver = solverr;

    if (solver != nullptr)
        CollectPrintableTcs(solver->TCScalculator);

}


TransportCrossSectionCsv::TransportCrossSectionCsv(
    TcsInterface* solverr,
    const std::string& folder
) {

    solver = solverr;
    customFolder = folder;

    if (solver != nullptr)
        CollectPrintableTcs(solver->TCScalculator);

}


void TransportCrossSectionCsv::Print(
    const std::string& filename,
    const std::vector<double>& x,
    GasMixture* gasmix
) {

    if (printableTcs.empty()) {

        std::cout << "Skipping: no transport cross section calculator found for "
                  << filename << "." << std::endl;

        return;

    }


    bool firstBlock = true;

    for (const auto& printable : printableTcs) {

        currentCalculator = printable.calculator;
        currentThreshold = printable.threshold;
        currentParentCalculator = printable.parentCalculator;
        currentInelastic = printable.inelastic;

        data.clear();
        header.clear();

        PrepareHeader();
        PrepareData(x, gasmix);

        if (!printable.suffix.empty())
            header += ", (" + printable.suffix.substr(1) + ")";

        if (data.empty() || header.empty()) {

            std::cerr << "Error: data or header not prepared for "
                    << filename + printable.suffix << "." << std::endl;

            continue;

        }

        if (firstBlock) {
            WriteData(filename);
            firstBlock = false;
        } else {
            AppendData(filename);
        }

        currentCalculator = nullptr;
        currentThreshold = nullptr;
        currentParentCalculator = nullptr;
        currentInelastic = false;

    }

    PrintMessage(filename);

}

// ___________________________ CollisionIntegralCsv __________________________ //

std::string CollisionIntegralCsv::BuildFileName(const std::string& filename) const {
    return "CI_" + filename + ".csv";
}

void CollisionIntegralCsv::PrepareHeader() {
    header =
        "Tij* [K],"
        " πΩ(1;1)[Å^2], πΩ(1;2)[Å^2], πΩ(1;3)[Å^2], πΩ(1;4)[Å^2], πΩ(1;5)[Å^2], πΩ(1;6)[Å^2], πΩ(1;7)[Å^2],"
        " πΩ(2;2)[Å^2], πΩ(2;3)[Å^2], πΩ(2;4)[Å^2], πΩ(2;5)[Å^2], πΩ(2;6)[Å^2],"
        " πΩ(3;3)[Å^2], πΩ(3;4)[Å^2], πΩ(3;5)[Å^2],"
        " πΩ(4;4)[Å^2]";
}

void CollisionIntegralCsv::PrepareData(const std::vector<double>& x, GasMixture* gasmix) {
    data.resize(x.size(), std::vector<double>(17));

    // if gasmix is provided, store its initial temperature
    double T0 = 0.0;
    bool hadGas = (gasmix != nullptr);
    if (hadGas) T0 = x.front();

    for (size_t i = 0; i < x.size(); ++i) {
        const double T = x[i];

        double lambda = std::numeric_limits<double>::infinity(); // default: no Debye, nullsafe
        double Te     = T;                                       // default: Te = T

        if (hadGas) {
            gasmix->setT(T);
            lambda = gasmix->getCompositionObj()->getDebyeLength(T);
            Te     = T * gasmix->theta->get();
        }

        solver->ComputeCollisionIntegral(Te, T, lambda);

        data[i][0] = solver->TijStar;

        const auto& omega = solver->omega4th;
        for (size_t j = 0; j < omega.size(); ++j)
            data[i][j + 1] = omega[j];
    }

    if (hadGas) {
        gasmix->setT(T0);
        gasmix->restartComposition();
    }
}

void CollisionIntegralCsv::PrintMessage(const std::string& filename) {
    std::cout << "Collision Integral "
              << solver->InteractionName()
              << " printed (" << filename << ")." << std::endl;
}

CollisionIntegralCsv::CollisionIntegralCsv(CInterface* solver)
    : solver(solver) {}

CollisionIntegralCsv::CollisionIntegralCsv(CInterface* solver, const std::string& folder)
    : solver(solver) {
    customFolder = folder;
}

// ___________________________ PartitionFunctionCsv __________________________ //

std::string PartitionFunctionCsv::BuildFileName(const std::string& filename) const {
    return "PF_" + filename + ".csv";
}

void PartitionFunctionCsv::PrepareHeader() {
    header = "T [K], P [Pa], λD [m], Q [#]";
}

void PartitionFunctionCsv::PrepareData(const std::vector<double>& Ti, GasMixture* gasmix) {
    double T0 = gasmix->getTemperature();
    double P  = gasmix->getPressure();

    data.resize(Ti.size(), std::vector<double>(4));

    for (size_t i = 0; i < Ti.size(); ++i) {
        gasmix->setT(Ti[i]);
        double lambda = gasmix->getCompositionObj()->getDebyeLength(Ti[i]);

        data[i][0] = Ti[i];
        data[i][1] = P;
        data[i][2] = lambda;

        solver->computePartitionFunction(Ti[i], P, lambda);
        data[i][3] = solver->getPf();
    }

    // Restore initial mixture state
    gasmix->setT(T0);
    gasmix->restartComposition();
}

void PartitionFunctionCsv::PrintMessage(const std::string& filename) {
    std::cout << "Partition Function " << filename << " printed." << std::endl;
}

void PartitionFunctionCsv::Print(const std::string& filename, const std::vector<double>& x, GasMixture* gasmix) {
    
    DataPrinter::Print(filename, x, gasmix);
}

PartitionFunctionCsv::PartitionFunctionCsv(PFinterface* solver)
    : solver(solver) {}

PartitionFunctionCsv::PartitionFunctionCsv(PFinterface* solver, const std::string& folder)
    : solver(solver) {
    customFolder = folder;
}
