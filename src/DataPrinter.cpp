 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "DataPrinter.h"

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
