#include"ParfCalculator.h"
#include <filesystem>
#include <regex>
#include <stdexcept>
#include <numbers>

#ifdef PPFM_USE_CURL
    #include <curl/curl.h>
    void ElectronicAtomicPF::CsvDownloadFromNISTAtomicSpectraDatabase(const std::string& hyperref, 
        const std::string& filename) {
        
        CURL* curl;
        CURLcode res;
        std::stringstream raw_data;

        curl = curl_easy_init();
        if (curl) {

            curl_easy_setopt(curl, CURLOPT_URL, hyperref.c_str());
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, &raw_data);

            res = curl_easy_perform(curl);
            if (res != CURLE_OK) {
            
                curl_easy_cleanup(curl);
                throw std::runtime_error("CURL error: " + std::string(curl_easy_strerror(res)));
            
            }
            
            curl_easy_cleanup(curl);

        } else {
            
            throw std::runtime_error("Failed to initialize CURL.");

        }

        if (raw_data.str().empty()) {
            
            throw std::runtime_error("Downloaded content is empty. Check the provided link.");

        }

        // Apertura del file di output con precisione aumentata
        std::ofstream clean_file(filename, std::ios::out | std::ios::trunc);
        if (!clean_file.is_open()) {
            
            throw std::runtime_error("Unable to open file for writing: " + filename);

        }

        // Configurazione precisione
        clean_file << std::fixed << std::setprecision(8);

        // Scrive l'header corretto
        clean_file << "g,Level (eV)\n";

        std::string line;
        std::regex valid_line_regex(R"(^\s*(\d+)\t.*\t\"([\d\.]+)\"\t.*$)");
        std::smatch match;

        int valid_lines_count = 0;

        // Salta l'header
        std::getline(raw_data, line);

        // Parsing delle righe
        while (std::getline(raw_data, line)) {

            if (std::regex_match(line, match, valid_line_regex)) {
            
                try {
            
                    int g = std::stoi(match[1].str());
                    double level = std::stod(match[2].str());
                    clean_file << g << "," << level << "\n";
                    valid_lines_count++;
            
                } catch (const std::exception&) {
            
                    // Ignora righe con valori non convertibili
                    continue;
                }
            }
        }

        clean_file.close();

        // Controlla se ci sono righe valide
        if (valid_lines_count == 0) {
            
            throw std::runtime_error("No valid data found in the downloaded content. "
            "Ensure the link points to a valid NIST data file.");
        }

        const int minimum_valid_lines = 5;
        if (valid_lines_count < minimum_valid_lines) {
            
            throw std::runtime_error(
                "Insufficient valid data in the downloaded content. "
                "The link may not point to a complete or correct NIST data file."
            );
        }
    }

#else 

    void ElectronicAtomicPF::CsvDownloadFromNISTAtomicSpectraDatabase(const std::string& hyperref, 
        const std::string& filename) {

        std::cerr << "[PPFM] libcurl not available: DownloadFromNist() skipped.\n";
    }

#endif

double PfPoly::compute(double T, double P, double debye) {
    
    if (coefficients.size() == 0) {
        
        throw std::invalid_argument(
            "Invalid Partition Function Generation, \n"
            "\tcheck for hard-coded specializations in the HCParFun.tpp file \n"
            "\tfor polynomial approximation or try different method.\n"
            "\tsee ParfCalculator.h for methods"
        ); 
        abort();
    }

    double Qint = 0.;
    for (int i = coefficients.size() - 1, j = 0; i > -1 && j < coefficients.size(); i--, j++) 
        Qint += coefficients[j] * pow(T, i);
    
    return Qint;
}

void PfTtable::ParseFile(std::ifstream& file) {

    std::string line;
    // Skip header
    std::getline(file, line); 
    
    double temp, q;

    while (std::getline(file, line)) {
    
        std::stringstream ss(line);
        std::string temp_str, q_str;
        
        if (std::getline(ss, temp_str, ',') && std::getline(ss, q_str, ',')) {
    
            try {
    
                temp = std::stod(temp_str);
                q = std::stod(q_str);
                temperatures.push_back(temp);
                QQ.push_back(q);
    
            } catch (...) {
    
                throw std::runtime_error("Error in T table data");
    
            }
        }
    }
}

// Constructor without prefix
PfTtable::PfTtable(Species* sp) {

    LoadData("Partition_Functions", sp->getFormula());
    this->metodo = "Table(T)";

}

// Constructor with prefix
PfTtable::PfTtable(Species* sp, const std::string& prefix) {

    this->customPrefix = prefix + "_";
    LoadData("Partition_Functions", sp->getFormula());
    this->metodo = "Table(T)";

}

double PfTtable::compute(double T, double P, double debye) {

    std::vector<double> logQQ(temperatures.size());

    std::transform(QQ.begin(), QQ.end(), logQQ.begin(), [](double val) {

        if (val > 0) return std::log(val);

        throw std::runtime_error("Negative value in partition function");

    });

    return std::exp(interpolateSpline(temperatures, logQQ, T));

}

void PfTPtable::ParseFile(std::ifstream& file) {

    std::string line;

    std::getline(file, line); // Skip header
    std::getline(file, line); // Read pressure values
    
    std::stringstream ss(line);
    std::string cell;
    std::getline(ss, cell, ','); // Skip first label cell

    while (std::getline(ss, cell, ',')) 
        pressures.push_back(std::stod(cell) * 1e5); // Convert bar to Pascal
    
    while (std::getline(file, line)) {
    
        std::stringstream ss(line);
        std::vector<double> qs_row;
        std::getline(ss, cell, ',');
        temperatures.push_back(std::stod(cell));
        while (std::getline(ss, cell, ',')) 
            qs_row.push_back(std::stod(cell));
        QQ2D.push_back(qs_row);
    }
}

// Constructor without prefix
PfTPtable::PfTPtable(Species* sp) {

    LoadData("Partition_Functions", sp->getFormula());
    this->metodo = "Table(T,P)";

}

// Constructor with prefix
PfTPtable::PfTPtable(Species* sp, const std::string& prefix) {

    this->customPrefix = prefix + "_";
    LoadData("Partition_Functions", sp->getFormula());
    this->metodo = "Table(T,P)";

}

double PfTPtable::compute(double T, double P, double debye) {

    std::vector<double> QonTemp(pressures.size());

    for (int j = 0; j < pressures.size(); j++) {

        std::vector<double> logQQ(temperatures.size());

        for (int i = 0; i < temperatures.size(); i++) {

            if (QQ2D[i][j] > 0)
                logQQ[i] = std::log(QQ2D[i][j]);
            else
                throw std::runtime_error("Negative value in partition function");
        
        }
        
        QonTemp[j] = std::exp(interpolateSpline(temperatures, logQQ, T));
    
    }
    
    return interp(pressures, QonTemp, P);

}

std::string PfTtable::BuildFileName(const std::string& name) {

    return "PF_" + customPrefix + name + "_T.csv";

}

std::string PfTPtable::BuildFileName(const std::string& name) {

    return "PF_" + customPrefix + name + "_TP.csv";

}


// _________________________ ElectronicAtomicPf Implementation _________________________

void ElectronicAtomicPF::ParseFile(std::ifstream& file) {

    std::string line;
    std::getline(file, line); // Skip header
    
    std::vector<double> gs;
    std::vector<double> energies;
    
    while (std::getline(file, line)) {

        std::stringstream ss(line);
        std::string gStr, energyStr;
        
        if (std::getline(ss, gStr, ',') && std::getline(ss, energyStr, ',')) {

            try {

                double g = std::stod(gStr);
                double energy = std::stod(energyStr);
                gs.push_back(g);
                energies.push_back(energy * eVtoJ); // Convert eV to J

            } catch (...) {

                throw std::runtime_error("Error in electronic configuration data");

            }
        }
    }
    
    Levels.resize(gs.size());
    for (int i = 0; i < gs.size(); ++i) 
        Levels[i] = {gs[i], energies[i]};
}

// Constructor without prefix
ElectronicAtomicPF::ElectronicAtomicPF(Species* sp) : sp(sp) {

    LoadData("Electronic_Configurations", sp->getFormula());
    this->metodo = "From ElConfig";

}

// Constructor with prefix
ElectronicAtomicPF::ElectronicAtomicPF ( const std::string& prefix, Species* sp ) : sp(sp) {

    this->customPrefix = prefix + "_";
    LoadData("Electronic_Configurations", sp->getFormula());
    this->metodo = "From ElConfig";

}

// Constructor for downloading data from NIST
ElectronicAtomicPF::ElectronicAtomicPF(Species* sp, const std::string& hyperref ) : sp(sp) {

    std::string speciesName = sp->getFormula();
    std::filesystem::path datadir = std::filesystem::current_path() / "data/Electronic_Configurations";
    std::string filename = (datadir / (speciesName + "_ElConfig.csv")).string();
    
    try {

        CsvDownloadFromNISTAtomicSpectraDatabase(hyperref, filename);
        new (this) ElectronicAtomicPF(sp); // Call standard constructor

    } catch (const std::exception& e) {

        throw std::runtime_error("Error in downloading or processing electronic configuration: " + std::string(e.what()));

    }
}

double ElectronicAtomicPF::compute(double T, double P, double debye) {

    // From Capitelli Fundamental Aspects of Chemical Plasma Physics (FACPP)

    // FACPP formula 8.3: Lowering Ionization Potential
    double DeltaIs = (std::pow(qe,2)*(sp->getCharge()+1)) / (4*std::numbers::pi*eps0*debye);

    // Debye-Hückel cutoff for partition function
    double epsMax = sp->IonLim() - DeltaIs;
    
    // FACPP formula 4.17
    double Qint = 0; 
    int i = 0;
    
    while (i < Levels.size() && Levels[i][1] < epsMax) {

        Qint += Levels[i][0] * exp(-Levels[i][1] / (KB*T));
        i++;

    }
    
    // FACPP formula 4.18
    return Qint;
}

std::string ElectronicAtomicPF::BuildFileName(const std::string& name) {

    return customPrefix + name + "_ElConfig.csv";

}

// _________________________ Hard-coded polinomial partition functions _________________________

// Ar
PfPoly::PfPoly(Argon* argon ) {
    Init( { 
        9.486e-31 , -1.561e-25 ,  9.383e-21 ,
        -2.399e-16 ,  2.929e-12 , -1.694e-08 ,
        4.02e-05 , 0.9746 
    } );
}
// Ar+
PfPoly::PfPoly(ArgonI* argonI ) {
    Init ( { 
        5.469e-31, -1.026e-25, 7.884e-21,
        -3.129e-16, 6.956e-12, -8.793e-08,
        0.0006138, 3.764 
    } );
}
// Ar+2
PfPoly::PfPoly(ArgonII* argonII ) {
    Init ( { 
        4.783e-39,  -1.071e-33 , 1.023e-28 ,
        -5.438e-24 , 1.761e-19 , -3.574e-15,
        4.504e-11 , -3.414e-07 , 0.001633  ,
        4.425
    } );
}
// e-
PfPoly::PfPoly(Electron* argonII ) {
    Init ( { 
        2
    } );
}
// H
PfPoly::PfPoly(Hydrogen* elettrone ) {
    Init ( { 

        -1.053e-21 , 1.093e-16 , -2.737e-12 , 
        2.489e-08 , -7.648e-05 , 2.048
        
    } );
}
// N2
PfPoly::PfPoly(MolecularNitrogen* ptr ) {
    Init ( { 
        1.9398e-61  , -6.5948e-56 , 1.0136e-50 , 
        -9.3112e-46 , 5.6905e-41  , -2.4344e-36 , 
        7.4594e-32  , -1.6446e-27 , 2.5844e-23 , 
        -2.8421e-19 , 2.1549e-15  , -1.1051e-11 , 
        3.7645e-08 , -2.5860e-05 , 1.7916e-01 , 
        1.7105e-01 
    } );
}
// N
PfPoly::PfPoly(Nitrogen* ptr ) {
    Init ( { 
            1.0640e-54 , -2.7824e-49 ,  3.2152e-44 , 
        -2.1701e-39 ,  9.5130e-35 , -2.8401e-30 , 
            5.8671e-26 , -8.3264e-22 ,  7.9088e-18 , 
        -4.8090e-14 ,  1.7522e-10 , -3.3505e-07 , 
            2.6244e-04 ,  3.9470e+00
    } );
}
// N+
PfPoly::PfPoly(NitrogenI* ptr ) {
    Init ( { 
        -3.8088e-50 ,  1.0557e-44 , -1.2818e-39 ,
            8.9679e-35 , -4.0007e-30 ,  1.1896e-25 ,
        -2.3937e-21 ,  3.2435e-17 , -2.8896e-13 ,
            1.6147e-09 , -5.2124e-06 ,  8.4801e-03 ,
            3.5065e+00
    } );
}
// N+2
PfPoly::PfPoly(NitrogenII* ptr ) {
    Init ( { 
                    
            2.2272e-54 , -6.5250e-49 ,  8.4580e-44 ,
        -6.3954e-39 ,  3.1311e-34 , -1.0417e-29 ,
            2.4043e-25 , -3.8616e-21 ,  4.2646e-17 ,
        -3.1500e-13 ,  1.4843e-09 , -4.1373e-06 ,
            6.0587e-03 ,  2.0937e+00          
                
    } );
}
// N+3
PfPoly::PfPoly(NitrogenIII* ptr ) {
    Init ( { 
        -4.9536e-19 ,  5.0112e-14 , -7.5443e-10 ,  
            2.8509e-06 ,  9.9856e-01    
    } );
}
// H2
PfPoly::PfPoly(MolecularHydrogen* elettrone ) {
    Init ( { 
            -2.169e-24 , 2.831e-19 , -1.348e-14 , 
            2.473e-10 , 2.717e-08 ,  0.004224  ,  
            2.619
    } );
}
// H+
PfPoly::PfPoly(HydrogenI* elettrone ) {
    Init ( { 
        1.
    } );
}
// H-
PfPoly::PfPoly(HydrogenAnion* ptr ) {
    Init ( { 
        1.
    } );
}
