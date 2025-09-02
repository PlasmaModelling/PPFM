 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef DATAPRINTER_H
#define DATAPRINTER_H

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

class GasMixture ; 

class DataPrinter {

    protected:

    /// @brief Data matrix to print in Csv file    
    std::vector<std::vector<double>> data ;

    /** @brief Header file to construct in the printable object, 
     remember to use "," as separator for the printables.  */
    std::string header ;

    /** @brief Implement to return a Default filename with coherent prefixes and extensions.  
     * @param filename Reference string for the file name ex: species->getFormula().  
    */
    virtual std::string BuildFileName( const std::string& filename ) const = 0 ;

    /** @brief Return the complete absolute path ( default: PathTo/PPFMfolder/out/" directory ). 
     * Override to specify other path when deriving this class. 
     * * @param filename Reference string for the file name ex: species->getFormula(). */
    virtual std::string BuildFilePath( const std::string& filename ) const ;

    /** @brief Override to implement data calculations and storing in this class member 
    vector<vector<double>> data. 
    * @param x Independent variable 
    * @param gasmix Pointer to class GasMixture  */
    virtual void PrepareData(const std::vector<double>& x, GasMixture* gasmix) = 0 ;

    // Override to implement the header to the data file
    virtual void PrepareHeader() = 0 ;

    /// @brief FileWriting centralized method.
    void WriteData(const std::string& filename ) ; 
    
    /// @brief FileAppendingCentralizedMethod
    void AppendData(const std::string& filename) ;

    public:

    /// @brief customizable folder, just assign it implementing constructor in derived classes
    std::string customFolder ;
    
    /// @brief Centralized printing method, override to customize Printing
    virtual void Print ( const std::string& filename, const std::vector<double>& x, GasMixture* gasmix ) ;

    virtual void PrintMessage(const std::string& filename) = 0 ; 

};


#endif