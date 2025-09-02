 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "DataLoader.h"

void DataLoader::LoadData(const std::string& folderName, const std::string& name) {

    try {

        std::filesystem::path filePath = BuildFilePath(folderName);
        std::string fileName = BuildFileName(name);
        std::filesystem::path fullPath = filePath / fileName;  

        std::ifstream file(fullPath);
        if (!file.is_open()) 
            throw std::runtime_error("File " + fileName + " not found in:\n" + filePath.string());

        ParseFile ( file );
        file.close();

    } catch (const std::exception& e) {

        throw std::runtime_error("Error loading file: " + std::string(e.what()));
    
    }

    loaded = true ; 

}

std::filesystem::path DataLoader::BuildFilePath(const std::string& folderName) {
    
    std::filesystem::path datadir = std::filesystem::current_path();
    
    datadir = datadir / "data" / folderName;
    
    return datadir;

}
