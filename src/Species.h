 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#ifndef SPECIES_H
#define SPECIES_H

#include "Tools.h"
#include <variant>
#include <stdexcept>

/**
 * @file Species.h
 * @brief This file contains the Species class hierarchy, 
 * which represents different chemical species and their properties.
 * The Species class is an abstract base class that defines the interface for all species,
 * then, species separates in Elements, which are neutral and basis for plasma composition, 
 * Charged Species, which are charged and have a reference to their neutral constituent,
 * and molecules, divided into polyatomic, diatomic etc... 
*/

#define amuKg 1.66054e-27 
#define eVtoJ 1.60217663e-19 
#define eamu  5.4858e-4 

/**
 * @brief Abstract base class for representing chemical species.
 * @details Defines the interface for all chemical species, 
 * including their properties and behavior. 
 * Each chemical specie is implemented as a class and realize interface 
 * methods when defining mass, charge, energy, ionization limit etc...
*/
class Species {

    public:
    
    /// virtual method to get the chemical formula of the species
    virtual std::string getFormula()        = 0 ;
    
    /// @brief virtual method to get the mass of the species in kg 
    virtual double      getMass()           = 0 ;           ///< unit: kg 

    /// @brief virtual method to get the charge of the species in elementary charge units   
    virtual int         getCharge()         { return 0 ; }  ///< unit: # 
    
    /// @brief virtual method to get the formation energy of the species in joules
    virtual double      formationEnergy()   = 0 ;           ///< unit: J 

    /// @brief virtual method to get the ionization limit of the species in joules
    virtual double      IonLim()            = 0 ;           ///< unit: J

} ;

/** @brief Class representing a chemical element.
 * @details Elements are neutral species and serve as the basis for plasma composition.
*/
class Element : public Species {};

/** @brief Class representing a charged chemical species.
 * @details Charged species are particles that have a net electric charge
 *  and are derived from neutral elements.
*/
class ChargedSpecies : public virtual Species {

    public : 

    /** @brief virtual method to get the neutral constituent of the charged species
     * @return pointer to the neutral constituent element, fundamental for chemical 
     * composition algorithms.
    */
    virtual Element* Constituent() = 0 ;
    
} ;

/// @brief Class representing a generic polyatomic chemical molecule.
class PolyAtomicMolecule : public virtual Species {

    protected : 

    /// @brief vector of pointers to the constituent elements of the molecule.
    std::vector<Element*>   costituents ; 

    /// @brief vector of integers representing the abundancy of each constituent element in the molecule.
    std::vector<int>        abundancy ;

    public : 

    /// @brief virtual method to get the number of constituent elements in the molecule
    virtual int             numberOfCostituents()   = 0 ;

    /// @brief virtual method to get the abundancy of the i-th constituent element in the molecule
    virtual int             operator()(int i )      {return abundancy[i] ; } ;

    /// @brief virtual method to get a pointer to the i-th constituent element in the molecule
    virtual Element*        operator[](int i )      {return costituents[i] ;}

} ;

/** @brief Class representing a generic biatomic chemical molecule, derived from PolyAtomicMolecule.
 * @details class meant to lock the number of costituents to 2, 
 * and to implement the operator() and operator[] methods so that derived classes do not have to. 
*/
class BiatomicMolecule : public PolyAtomicMolecule { 

    public : 
    
    int     operator()(int i)   override { 
        if ( i < 0 || i > 1 )
            throw std::invalid_argument("Error in Species costituents") ; 
        else
            return abundancy[i] ; 
    } 

    Element* operator[](int i)  override {
        if ( i < 0 || i > 1 )
            throw std::invalid_argument("Error in Species costituents") ; 
        else
            return costituents[i] ; 
    }

} ;

/**
 * @brief Class representing a homonuclear biatomic molecule, derived from BiatomicMolecule.
 * @details This class represents a molecule composed of two identical atoms.
 * It is intended to lock the number of constituents to 1 and to set the abundancy to 2.
*/
class HomoNuclearBiatomicMolecule : public BiatomicMolecule {
    
    // behavior of energy levels ond so on 
    public : 
    
    HomoNuclearBiatomicMolecule() {
        costituents.resize(1) ;
        abundancy.resize(1) ; 
        abundancy[0] = 2 ;
    }

    int numberOfCostituents()   override { 
        return  1 ; 
    } 

} ;

/**
 * @brief Class representing a heteronuclear biatomic molecule, derived from BiatomicMolecule.
 * @details This class represents a molecule composed of two different atoms.
 * It is intended to lock the number of constituents to 2 and to set the abundancy to 1 for each constituent.
*/
class HeteroNuclearBiatomicMolecule : public BiatomicMolecule {

    public : 
    
    HeteroNuclearBiatomicMolecule() {
        costituents.resize(2) ;
        abundancy.resize(2) ; 
        abundancy[0] = 1 ;
        abundancy[1] = 1 ;
    }

    int numberOfCostituents()   override { 
        return  2 ; 
    } 

} ;


/// @brief e-
class Electron : public ChargedSpecies {
    
    public:
    
    double      getMass()           override { return 9.1093837015e-31 ; }
    int         getCharge()         override { return -1 ; }
    Element*    Constituent()       override { return nullptr ; }
    std::string getFormula()        override { return "e-"; } 
    double      formationEnergy()   override { return 0.0;}
    double      IonLim()            override { std::cerr<<"Asked for electron ionization limit. aborting...\n";
        abort() ; return 999999999999 ; }

} ;

class Aluminum : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    
    double IonLim() override ;
} ;

class AluminumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class AluminumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class AluminumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Argon : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class ArgonI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ArgonII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ArgonIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ArgonIV : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Arsenic : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class ArsenicI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ArsenicII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ArsenicIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Barium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class BariumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class BariumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class BariumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Beryllium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class BerylliumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class BerylliumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class BerylliumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Boron : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class BoronI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class BoronII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class BoronIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Bromine : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class BromineI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class BromineII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class BromineIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Calcium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class CalciumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class CalciumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class CalciumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Carbon : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class CarbonI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class CarbonII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class CarbonIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;

class CarbonIV : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;

class Cesium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class CesiumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class CesiumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class CesiumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Chlorine : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class ChlorineI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ChlorineII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ChlorineIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Chromium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class ChromiumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ChromiumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ChromiumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Cobalt : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class CobaltI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class CobaltII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class CobaltIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Copper : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class CopperI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class CopperII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class CopperIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Fluorine : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class FluorineI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class FluorineII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class FluorineIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Gallium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class GalliumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class GalliumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class GalliumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Germanium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class GermaniumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class GermaniumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class GermaniumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Gold : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class GoldI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class GoldII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class GoldIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Helium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class HeliumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class HeliumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Hydrogen : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class HydrogenI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Iodine : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class IodineI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class IodineII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class IodineIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Iron : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class IronI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class IronII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class IronIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Krypton : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class KryptonI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class KryptonII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class KryptonIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;

class KryptonIV : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Lead : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class LeadI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class LeadII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class LeadIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Lithium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class LithiumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class LithiumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class LithiumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Magnesium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class MagnesiumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class MagnesiumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class MagnesiumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Manganese : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class ManganeseI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ManganeseII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ManganeseIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Mercury : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class MercuryI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class MercuryII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class MercuryIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Molybdenum : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class MolybdenumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class MolybdenumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class MolybdenumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Neon : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class NeonI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class NeonII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class NeonIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Nickel : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class NickelI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class NickelII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class NickelIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Niobium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class NiobiumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class NiobiumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class NiobiumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Nitrogen : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class NitrogenI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class NitrogenII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class NitrogenIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Oxygen : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class OxygenI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class OxygenII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class OxygenIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;

class OxygenIV : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;

class Phosphorus : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class PhosphorusI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class PhosphorusII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class PhosphorusIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Platinum : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class PlatinumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class PlatinumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class PlatinumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Potassium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class PotassiumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class PotassiumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class PotassiumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Radon : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class RadonI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class RadonII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class RadonIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Rhenium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class RheniumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class RheniumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class RheniumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Rubidium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class RubidiumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class RubidiumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class RubidiumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Silicon : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class SiliconI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class SiliconII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class SiliconIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Silver : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class SilverI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class SilverII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class SilverIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Sodium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class SodiumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class SodiumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class SodiumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Strontium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class StrontiumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class StrontiumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class StrontiumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Sulfur : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class SulfurI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class SulfurII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class SulfurIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Tantalum : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class TantalumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class TantalumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class TantalumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Titanium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class TitaniumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class TitaniumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class TitaniumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Tungsten : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class TungstenI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class TungstenII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class TungstenIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Vanadium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class VanadiumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class VanadiumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class VanadiumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Xenon : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class XenonI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class XenonII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class XenonIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;

class XenonIV : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;

class Zinc : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class ZincI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ZincII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ZincIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class Zirconium : public Element {
public:
    double getMass() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
} ;


class ZirconiumI : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ZirconiumII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class ZirconiumIII : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class HydrogenAnion : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class NitrogenAnion : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    double IonLim() override ;
    Element* Constituent() override ;
} ;


class OxygenAnion : public ChargedSpecies {
public:
    double getMass() override ;
    int getCharge() override ;
    std::string getFormula() override ;
    double formationEnergy() override ;
    
    double IonLim() override ;
    Element* Constituent() override ;
} ;

class CarbonMonoxide : public HeteroNuclearBiatomicMolecule {
    
    public:
    
    CarbonMonoxide();

    int numberOfCostituents() override;
    int operator()(int i) override;
    Element* operator[](int i) override;

    double      getMass() override;
    int         getCharge() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    double      IonLim() override;

};

class CarbonMonoxideI : public HeteroNuclearBiatomicMolecule, public ChargedSpecies {
    
    public:
    
    CarbonMonoxideI();

    int numberOfCostituents() override;
    int operator()(int i) override;
    Element* operator[](int i) override;

    double      getMass() override;
    int         getCharge() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    double      IonLim() override;
    Element*    Constituent() override;

};

class CarbonDioxide : public PolyAtomicMolecule {
    
    public:
    
    CarbonDioxide();

    int numberOfCostituents() override;
    int operator()(int i) override;
    Element* operator[](int i) override;

    double      getMass() override;
    int         getCharge() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    double      IonLim() override;
};

class CarbonDioxideI : public PolyAtomicMolecule, public ChargedSpecies {
    
    public:
    
    CarbonDioxideI();

    int numberOfCostituents() override;
    int operator()(int i) override;
    Element* operator[](int i) override;

    double      getMass() override;
    int         getCharge() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    double      IonLim() override;
    Element*    Constituent() override;

};

class CarbonDioxideAnion : public PolyAtomicMolecule, public ChargedSpecies {
    
    public:
    
    CarbonDioxideAnion();

    int numberOfCostituents() override;
    int operator()(int i) override;
    Element* operator[](int i) override;

    double      getMass() override;
    int         getCharge() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    double      IonLim() override;
    Element*    Constituent() override;

};

class Ozone : public PolyAtomicMolecule {

    public:

    Ozone();

    int numberOfCostituents() override;
    int operator()(int i) override;
    Element* operator[](int i) override;

    double      getMass() override;
    int         getCharge() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    double      IonLim() override;
};

class MolecularHydrogen : public HomoNuclearBiatomicMolecule {
public:
    MolecularHydrogen();

    double      getMass() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    double      IonLim() override;
};

class MolecularNitrogen : public HomoNuclearBiatomicMolecule {
public:
    MolecularNitrogen();

    double      getMass() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    double      IonLim() override;
};

class MolecularOxygen : public HomoNuclearBiatomicMolecule {
public:
    MolecularOxygen();

    double      getMass() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    double      IonLim() override;
};

class MolecularHydrogenI : public HomoNuclearBiatomicMolecule, public ChargedSpecies {
public:
    MolecularHydrogenI();

    double      getMass() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    int         getCharge() override;
    Element*    Constituent() override;
    double      IonLim() override;
};

class MolecularOxygenI : public HomoNuclearBiatomicMolecule, public ChargedSpecies {
public:
    MolecularOxygenI();

    double      getMass() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    int         getCharge() override;
    Element*    Constituent() override;
    double      IonLim() override;
};

class MolecularOxygenAnion : public HomoNuclearBiatomicMolecule, public ChargedSpecies {
public:
    MolecularOxygenAnion();

    double      getMass() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    int         getCharge() override;
    Element*    Constituent() override;
    double      IonLim() override;
};

class MolecularNitrogenI : public HomoNuclearBiatomicMolecule, public ChargedSpecies {
public:
    MolecularNitrogenI();
    
    double      getMass() override;
    std::string getFormula() override;
    double      formationEnergy() override;
    int         getCharge() override;
    Element*    Constituent() override;
    double      IonLim() override;
};

class CarbonAnion : public ChargedSpecies {
public:
    double getMass() override;
    int getCharge() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
    Element* Constituent() override;
};

class DiCarbon : public HomoNuclearBiatomicMolecule {
public:
    DiCarbon();

    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class DiCarbonI : public HomoNuclearBiatomicMolecule, public ChargedSpecies {
public:
    DiCarbonI();

    double getMass() override;
    int getCharge() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
    Element* Constituent() override;
};

class DiCarbonAnion : public HomoNuclearBiatomicMolecule, public ChargedSpecies {
public:
    DiCarbonAnion();

    double getMass() override;
    int getCharge() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
    Element* Constituent() override;
};

class EthynilRadical : public PolyAtomicMolecule {
public:
    EthynilRadical();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class Acetylene : public PolyAtomicMolecule {
public:
    Acetylene();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class Ethylene : public PolyAtomicMolecule {
public:
    Ethylene();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class DiCarbonMonoxide : public PolyAtomicMolecule {
public:
    DiCarbonMonoxide();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class TriCarbon : public PolyAtomicMolecule {
public:
    TriCarbon();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class MethylidineRadical : public PolyAtomicMolecule {
public:
    MethylidineRadical();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class Methylene : public PolyAtomicMolecule {
public:
    Methylene();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class MethylRadical : public PolyAtomicMolecule {
public:
    MethylRadical();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class Methane : public PolyAtomicMolecule {
public:
    Methane();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class FormylRadical : public PolyAtomicMolecule {
public:
    FormylRadical();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class FormylRadicalI : public PolyAtomicMolecule, public ChargedSpecies {
public:
    FormylRadicalI();

    int numberOfCostituents() override;
    double getMass() override;
    int getCharge() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
    Element* Constituent() override;
};

class Formaldehyde : public PolyAtomicMolecule {
public:
    Formaldehyde();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class Water : public PolyAtomicMolecule {
public:
    Water();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class HydroPeroxyRadical : public PolyAtomicMolecule {
public:
    HydroPeroxyRadical();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class HydroxilRadical : public HeteroNuclearBiatomicMolecule {
public:
    HydroxilRadical();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class HydroxilRadicalI : public HeteroNuclearBiatomicMolecule, public ChargedSpecies {
public:
    HydroxilRadicalI();

    int numberOfCostituents() override;
    double getMass() override;
    int getCharge() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
    Element* Constituent() override;
};

class PropadieneDione : public PolyAtomicMolecule {
public:
    PropadieneDione();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class PropadieneDioneI : public PropadieneDione , public ChargedSpecies {
public:
    PropadieneDioneI();

    std::string getFormula() override;
    double formationEnergy() override;
    Element* Constituent() override;

};

class TetraCarbon : public PolyAtomicMolecule {
public: 
    TetraCarbon();
    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class PentaCarbon : public PolyAtomicMolecule {
public: 
    PentaCarbon();
    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class EsaCarbon : public PolyAtomicMolecule {
public: 
    EsaCarbon();
    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class CarbonTrioxide : public PolyAtomicMolecule {
public: 

    CarbonTrioxide() ; 
    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class CarbonTrioxideAnion : public CarbonTrioxide, public ChargedSpecies {
public: 
    
    CarbonTrioxideAnion() ; 

    int getCharge() override;
    std::string getFormula() override;
    double formationEnergy() override;
    Element* Constituent() override;
};
class CNCRadical : public PolyAtomicMolecule {
public:
    
    CNCRadical();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class Cyanogen : public PolyAtomicMolecule {
public:
    
    Cyanogen();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class CyanoEthynilRadical : public PolyAtomicMolecule {
public:
    
    CyanoEthynilRadical();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class TwoButyneDiNitrile : public PolyAtomicMolecule {
public:
    
    TwoButyneDiNitrile();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class HydrogenCyanide : public PolyAtomicMolecule {
public:
    
    HydrogenCyanide();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class CyanoRadical : public PolyAtomicMolecule {
public:
    
    CyanoRadical();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class IsoCyanatoRadical : public PolyAtomicMolecule {
public:
    
    IsoCyanatoRadical();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class Imidogen : public PolyAtomicMolecule {
public:
    
    Imidogen();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class NitrusOxide : public PolyAtomicMolecule {
public:
    
    NitrusOxide();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class Ammonia : public PolyAtomicMolecule {
public:
    
    Ammonia();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};

class NitricOxide : public PolyAtomicMolecule {
public:

    NitricOxide();

    int numberOfCostituents() override;
    double getMass() override;
    std::string getFormula() override;
    double formationEnergy() override;
    double IonLim() override;
};
class NitricOxideI : public NitricOxide, public ChargedSpecies {
public:
    NitricOxideI();

    int getCharge() override;
    std::string getFormula() override;
    double formationEnergy() override;
    Element* Constituent() override;
};

#endif
