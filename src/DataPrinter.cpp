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

    std::filesystem::path base = std::filesystem::current_path() / "out";

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

    data.resize(temperatureRange.size(), std::vector<double>(mix->getN() + 2));

    for (int i = 0; i < temperatureRange.size(); i++) {
        
        mix->setT(temperatureRange[i]);
        solver->CompositionSolve(mix, mix);

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
    header = "Th [K], λₑ [W/(m·K)], λₕ [W/(m·K)], μ [Pa·s], σ [S/m]";
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
    data.resize(temperatureRange.size(), std::vector<double>(5));

    for (size_t i = 0; i < temperatureRange.size(); ++i) {
        gasmix->setT(temperatureRange[i]);
        solver->computeTransport(gasmix);

        data[i][0] = temperatureRange[i];
        for (int j = 0; j < 4; ++j)
            data[i][j + 1] = solver->Tp[j];
    }

    gasmix->setT(T0);
    gasmix->restartComposition();
    solver->computeTransport(gasmix);
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

// ___________________________ TransportCrossSectionCsv __________________________ //

#include "TcsCalculator.h"

std::string TransportCrossSectionCsv::BuildFileName(const std::string& filename) const {
    return "TCS_" + filename + ".csv";
}

void TransportCrossSectionCsv::PrepareHeader() {
    header = "E [eV], Qe(1) [Å²], Qe(2) [Å²], Qe(3) [Å²], Qe(4) [Å²]";
}

void TransportCrossSectionCsv::PrepareInelasticHeader() {
    header = "E [eV], Qin(1) [Å²], Qin(2) [Å²], Qin(3) [Å²], Qin(4) [Å²]";
}

void TransportCrossSectionCsv::PrepareMultiHeader(MultiCs* cscalc, int i) {
    header = "state " + std::to_string(i + 1) + ", E [eV], Q(1) [Å²], Q(2) [Å²], Q(3) [Å²], Q(4) [Å²]";
}

void TransportCrossSectionCsv::PrepareData(const std::vector<double>&, GasMixture*) {
    auto& e = solver->TCScalculator->E;
    auto& q = solver->TCScalculator->Q;

    data.resize(e.size(), std::vector<double>(5));
    for (size_t i = 0; i < e.size(); ++i) {
        data[i][0] = e[i];
        for (size_t j = 0; j < 4; ++j)
            data[i][j + 1] = q[i][j];
    }
}

void TransportCrossSectionCsv::PrepareData(MultiCs* cscalc, int i) {
    auto cs = (*cscalc)[i];
    if (!cs) return;

    const auto& E = cs->E;
    const auto& Q = cs->Q;

    data.resize(E.size(), std::vector<double>(6, 0.0));
    for (size_t row = 0; row < E.size(); ++row) {
        data[row][0] = 0.0; // state placeholder
        data[row][1] = E[row];
        for (int j = 0; j < 4; ++j)
            data[row][2 + j] = Q[row][j];
    }
}

void TransportCrossSectionCsv::PrepareInelasticData(CsHolder* tcsElIn) {
    auto& e = tcsElIn->Qin->E;
    auto& q = tcsElIn->Qin->Q;

    data.resize(e.size(), std::vector<double>(5));
    for (size_t i = 0; i < e.size(); ++i) {
        data[i][0] = e[i];
        for (size_t j = 0; j < 4; ++j)
            data[i][j + 1] = q[i][j];
    }
}

void TransportCrossSectionCsv::Print(const std::string& filename, const std::vector<double>& x, GasMixture* gasmix) {
    if (!solver || !solver->TCScalculator)
        return;

    if (solver->TCScalculator->Q.empty())
        solver->TCScalculator->Compute();

    std::string tempFolder = customFolder;
    customFolder += "/TransportCrossSections_" + tempFolder;

    if (auto multi = dynamic_cast<MultiCs*>(solver->TCScalculator)) {
        for (int i = 0; i < multi->Size(); ++i) {
            PrepareMultiHeader(multi, i);
            PrepareData(multi, i);
            (i == 0) ? WriteData(filename) : AppendData(filename);
        }
    } 
    else if (auto holder = dynamic_cast<CsHolder*>(solver->TCScalculator)) {
        PrepareHeader();
        PrepareData(x, gasmix);
        WriteData(filename);

        PrepareInelasticHeader();
        PrepareInelasticData(holder);
        AppendData(filename);
    } 
    else {
        DataPrinter::Print(filename, x, gasmix);
    }

    PrintMessage(filename);
    customFolder = tempFolder;
}

void TransportCrossSectionCsv::PrintMessage(const std::string& filename) {
    std::cout << "Transport Cross Sections "
              << solver->GetIntInterface()->InteractionName()
              << " printed (" << filename << ")." << std::endl;
}

TransportCrossSectionCsv::TransportCrossSectionCsv(TcsInterface* solver)
    : solver(solver) {}

TransportCrossSectionCsv::TransportCrossSectionCsv(TcsInterface* solver, const std::string& folder)
    : solver(solver) {
    customFolder = folder;
}

// ___________________________ CollisionIntegralCsv __________________________ //

std::string CollisionIntegralCsv::BuildFileName(const std::string& filename) const {
    return "CI_" + filename + ".csv";
}

void CollisionIntegralCsv::PrepareHeader() {
    header =
        "Tij* [#],"
        " Ω(1;1)[#], Ω(1;2)[#], Ω(1;3)[#], Ω(1;4)[#], Ω(1;5)[#], Ω(1;6)[#], Ω(1;7)[#],"
        " Ω(2;2)[#], Ω(2;3)[#], Ω(2;4)[#], Ω(2;5)[#], Ω(2;6)[#],"
        " Ω(3;3)[#], Ω(3;4)[#], Ω(3;5)[#],"
        " Ω(4;4)[#]";
}

void CollisionIntegralCsv::PrepareData(const std::vector<double>& x, GasMixture* gasmix) {
    data.resize(x.size(), std::vector<double>(17));

    for (size_t i = 0; i < x.size(); ++i) {
        gasmix->setT(x[i]);

        double lambda = gasmix->getCompositionObj()->getDebyeLength(x[i]);
        double Te = x[i] * gasmix->theta->get();

        solver->ComputeCollisionIntegral(Te, x[i], lambda);

        data[i][0] = solver->TijStar;

        const auto& omega = solver->omega4th;
        for (size_t j = 0; j < omega.size(); ++j)
            data[i][j + 1] = omega[j];
    }

    // Reset state
    gasmix->setT(x.front());
    gasmix->restartComposition();
}

void CollisionIntegralCsv::PrintMessage(const std::string& filename) {
    std::cout << "Collision Integral "
              << solver->InteractionName()
              << " printed (" << filename << ")." << std::endl;
}

void CollisionIntegralCsv::Print(const std::string& filename, const std::vector<double>& x, GasMixture* gasmix) {
    std::string tempFolder = customFolder;
    customFolder += "/CollisionIntegrals_" + tempFolder;

    DataPrinter::Print(filename, x, gasmix);

    customFolder = tempFolder;
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
    std::string tempFolder = customFolder;
    customFolder += "/PartitionFunctions_" + tempFolder;

    DataPrinter::Print(filename, x, gasmix);

    customFolder = tempFolder;
}

PartitionFunctionCsv::PartitionFunctionCsv(PFinterface* solver)
    : solver(solver) {}

PartitionFunctionCsv::PartitionFunctionCsv(PFinterface* solver, const std::string& folder)
    : solver(solver) {
    customFolder = folder;
}
