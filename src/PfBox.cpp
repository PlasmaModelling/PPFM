 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "PfBox.h"
#include "PartitionFunction.h"
#include "Mixture.h"

PfBox::PfBox(Mixture* mix) {
    
    for (int i = 0; i < mix->cispecies.size() ; i++) {
        std::visit([this](auto&& sp1){
            PFinterface* partitionFunction = new PartitionFunction(sp1) ; 
            partitionfunctions.push_back(partitionFunction) ;
        }, mix->cispecies[i] );
    }   
}

PFinterface* PfBox::operator[](int i){
    if (i<0 || i>= partitionfunctions.size())
        throw std::out_of_range("Index out of PfBox::partitionfunctions range") ;
    else
        return partitionfunctions[i] ;    
}

double PfBox::operator()(int i){
    if (i<0 || i>= partitionfunctions.size())
        throw std::out_of_range("Index out of PfBox::partitionfunctions range") ;
    else
        return partitionfunctions[i]->getPf() ;     
}

void PfBox::PrintPartitionFunctions(const std::vector<double>& Ti, GasMixture* gasmix, const std::string& folder) {
 
    for (auto& pf : partitionfunctions) {

        PartitionFunctionCsv writer(pf);

        writer.customFolder = folder;

        writer.Print(pf->getSp()->getFormula(), Ti, gasmix);

    }

}

void PfBox::computePartitionFunctions(double temperature, double pressure, double lambda) {

    bool exception_occurred = false;
    std::vector<std::string> error_messages(partitionfunctions.size());

    #pragma omp parallel for
    for (int i = 0; i < partitionfunctions.size(); i++) {
        try {

            partitionfunctions[i]->computePartitionFunction(temperature, pressure, lambda);

        } catch (const std::exception& e) {

            #pragma omp critical
            {
                exception_occurred = true;
                error_messages[i] = "Error in the " + std::to_string(i) + "-th partition function: "
                                    + partitionfunctions[i]->getSp()->getFormula() + "\n\n\t" + e.what() + "\n\n"
                                    + "- see PartitionFunction.h interface to choose for methods.\n";
            }
        }
    }

    if (exception_occurred) {
        for (const auto& msg : error_messages) {
            if (!msg.empty()) {
                std::cerr << msg << std::endl;
            }
        }
        std::exit(EXIT_FAILURE);
    }
}


void PfBox::info() {
    
    std::cout << std::left;
    std::cout<<"[i] | Partition | Method\n    | Function  |"<< std::endl;
    std::cout<<"‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"<< std::endl;
    for (int i = 0; i < partitionfunctions.size(); i++) {
        std::cout << std::setw(6) << i; partitionfunctions[i]->info() ;
    }
    std::cout << std::endl ;
}

void PfBox::AllAbInitio() { 
    
    for ( int i = 0; i < partitionfunctions.size()-1; i++ )         
        partitionfunctions[i]->setAbInitio(); 

}