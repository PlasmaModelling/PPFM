#ifndef DATALOADER_H
#define DATALOADER_H

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
/**
 * @brief Abstract base class for loading data from files
 * @details Provides a generic interface to build file paths, parse data, 
 * and handle loading logic. Derived classes must implement file naming and parsing. */
class DataLoader {
    
    protected:

    /// @brief Flag indicating whether the data has been successfully loaded
    bool loaded = false;  

    /// @brief Prefix used for custom file naming in derived classes
    std::string customPrefix; 

    /**
     * @brief Initializes internal state (to be defined in derived classes) */
    virtual void Init() = 0;

    /**
     * @brief Constructs the full path to the data folder
     * @param folderName Name of the subfolder inside the "data" directory
     * @return Full filesystem path */
    virtual std::filesystem::path BuildFilePath(const std::string& folderName);

    /**
     * @brief Constructs the file name from a given base name
     * @param name Base name to complete into a file name
     * @return File name string */
    virtual std::string BuildFileName(const std::string& name) = 0;

    /**
     * @brief Parses the file content
     * @param file Input file stream */
    virtual void ParseFile(std::ifstream& file) = 0;

    public:
    /**
     * @brief Loads data from a file given folder and base name
     * @param folderName Name of the folder inside the "data" directory
     * @param name Base name used to construct the file name
     * @details Builds the full file path, opens the file, and calls ParseFile. 
     * Throws if the file is not found or parsing fails. */
    void LoadData(const std::string& folderName, const std::string& name);

};

#endif
