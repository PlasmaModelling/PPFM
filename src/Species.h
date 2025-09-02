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


#define amuKg 1.66054e-27 
#define eVtoJ 1.60217663e-19 
#define eamu  5.4858e-4 

class Species {

    public:
    
    ///@todo implement electronic configuration placeholder

    virtual std::string getFormula()        = 0 ;
    virtual double      getMass()           = 0 ;           ///< unit: kg 
    virtual int         getCharge()         { return 0 ; }  ///< unit: # 
    virtual double      formationEnergy()   = 0 ;           ///< unit: J 
    virtual double      IonLim()            = 0 ;           ///< unit: J
} ;

class Element : public Species {};

class ChargedSpecies : public virtual Species {

    public : 

    virtual Element* Constituent() = 0 ;
    
} ;

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


class PolyAtomicMolecule : public virtual Species {

    protected : 

    std::vector<Element*>   costituents ; 
    std::vector<int>        abundancy ;

    public : 

    virtual int             numberOfCostituents()   = 0 ;
    virtual int             operator()(int i )      = 0 ;
    virtual Element*        operator[](int i )      = 0 ; 

} ;

class Ozone : public PolyAtomicMolecule{

    public:

    Ozone(){
        costituents.push_back(new Oxygen) ; 
        abundancy.push_back(3) ;
    }
    int numberOfCostituents() override {return 3;}
    int operator()(int i) override {return abundancy[0];}
    Element* operator[](int i) override {return costituents[0];}
    
    double      getMass()           override { return (47.9982) * amuKg ; } 
    int         getCharge()         override { return 0 ; }
    std::string getFormula()        override { return "O3" ; }
    double      formationEnergy()   override { return (-1.478671) * eVtoJ ; }
    double      IonLim()            override { return 999999999999 ; }

} ; 

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

class MolecularHydrogen : public HomoNuclearBiatomicMolecule {
    
    public : 

    MolecularHydrogen() {
        costituents[0] = new Hydrogen ;
    }

    double      getMass()           override { return 2.01588 * amuKg ; }
    std::string getFormula()        override { return "H2" ; }
    double      formationEnergy()   override { return -7.174875e-19 /* -7.240152773e-19 */ ; }
    double      IonLim()            override { return 999999999999 ; }

} ; 

class MolecularNitrogen : public HomoNuclearBiatomicMolecule {
    
    public : 

    MolecularNitrogen() {
        costituents[0] = new Nitrogen ;
    }

    double      getMass()           override { return 28.0134 * amuKg ; }
    std::string getFormula()        override { return "N2" ; }
    double      formationEnergy()   override { return -1.563614e-18 /* -7.240152773e-19 */ ; }
    double      IonLim()            override { return 999999999999 ; }

} ; 

class MolecularOxygen : public HomoNuclearBiatomicMolecule {
    
    public : 

    MolecularOxygen() {
        costituents[0] = new Oxygen ;
    }

    double      getMass()           override { return 31.9988 * amuKg ; }
    std::string getFormula()        override { return "O2" ; }
    double      formationEnergy()   override { return -5.11672987 * eVtoJ ; }
    double      IonLim()            override { return 999999999999 ; }

} ; 

// charged molecules
class MolecularHydrogenI : public HomoNuclearBiatomicMolecule, public ChargedSpecies {
    
    public : 

    MolecularHydrogenI() {
        costituents[0] = new Hydrogen ;
    }

    double      getMass()           override { return ( 2.01588 - eamu ) * amuKg ; }
    std::string getFormula()        override { return "H2+" ; }
    double      formationEnergy()   override { return -7.174875e-19 + (15.42593 * eVtoJ) ; }
    
    int         getCharge()         override { return 1 ; }
    Element*    Constituent()       override { return new Hydrogen ; }
    double      IonLim()            override { return 999999999999 ; }
    
} ; 


// charged molecules
class MolecularOxygenI : public HomoNuclearBiatomicMolecule, public ChargedSpecies {
    
    public : 

    MolecularOxygenI() {
        costituents[0] = new Oxygen ;
    }

    double      getMass()           override { return ( 31.9988 - eamu ) * amuKg ; }
    std::string getFormula()        override { return "O2+" ; }
    double      formationEnergy()   override { return (new MolecularOxygen)->formationEnergy() + (12.0697 * eVtoJ); }
    
    int         getCharge()         override { return 1 ; }
    Element*    Constituent()       override { return new Oxygen ; }
    double      IonLim()            override { return 999999999999 ; }
    
} ; 

class MolecularOxygenAnion : public HomoNuclearBiatomicMolecule, public ChargedSpecies {
    
    public : 

    MolecularOxygenAnion() {
        costituents[0] = new Oxygen ;
    }

    double      getMass()           override { return ( 31.9988 + eamu ) * amuKg ; }
    std::string getFormula()        override { return "O2-" ; }
    double      formationEnergy()   override { return (new MolecularOxygen)->formationEnergy() + (-0.448 * eVtoJ); }
    
    int         getCharge()         override { return -1 ; }
    Element*    Constituent()       override { return new Oxygen ; }
    double      IonLim()            override { return 999999999999 ; }
    
} ; 


// charged molecules
class MolecularNitrogenI : public HomoNuclearBiatomicMolecule, public ChargedSpecies {
    
    public : 

    MolecularNitrogenI() {
        costituents[0] = new Nitrogen ;
    }

    double      getMass()           override { return (new MolecularNitrogen)->getMass()-(eamu*amuKg) ; }
    std::string getFormula()        override { return "N2+" ; }
    double      formationEnergy()   override { return (new MolecularNitrogen)->formationEnergy() + (15.581*eVtoJ); }
    
    int         getCharge()         override { return 1 ; }
    Element*    Constituent()       override { return new Nitrogen ; }
    double      IonLim()            override { return 999999999999 ; }
    
} ; 

// using AcceptedSpecies = 

//     std::variant

//         <   
//             /* Electron */

//             /* Aluminum*, AluminumI*, AluminumII*, AluminumIII*, */
            
//             /* Argon*, ArgonI*, ArgonII*, ArgonIII*, ArgonIV*, */
            
//             /* Arsenic*, ArsenicI*, ArsenicII*, ArsenicIII*,
            
//             Barium*, BariumI*, BariumII*, BariumIII*,
            
//             Beryllium*, BerylliumI*, BerylliumII*, BerylliumIII*,
            
//             Boron*, BoronI*, BoronII*, BoronIII*,
            
//             Bromine*, BromineI*, BromineII*, BromineIII*,
            
//             Calcium*, CalciumI*, CalciumII*, CalciumIII*,
            
//             Carbon*, CarbonI*, CarbonII*, CarbonIII*,
            
//             Cesium*, CesiumI*, CesiumII*, CesiumIII*,
            
//             Chlorine*, ChlorineI*, ChlorineII*, ChlorineIII*,
            
//             Chromium*, ChromiumI*, ChromiumII*, ChromiumIII*,
            
//             Cobalt*, CobaltI*, CobaltII*, CobaltIII*,
            
//             Copper*, CopperI*, CopperII*, CopperIII*,
            
//             Fluorine*, FluorineI*, FluorineII*, FluorineIII*,
            
//             Gallium*, GalliumI*, GalliumII*, GalliumIII*,
            
//             Germanium*, GermaniumI*, GermaniumII*, GermaniumIII*,
            
//             Gold*, GoldI*, GoldII*, GoldIII*,
            
//             Helium*, HeliumI*, HeliumII*,
            
//             Hydrogen*, HydrogenI*,
            
//             Iodine*, IodineI*, IodineII*, IodineIII*,
            
//             Iron*, IronI*, IronII*, IronIII*,
            
//             Krypton*, KryptonI*, KryptonII*, KryptonIII*,
            
//             Lead*, LeadI*, LeadII*, LeadIII*,
            
//             Lithium*, LithiumI*, LithiumII*, LithiumIII*,
            
//             Magnesium*, MagnesiumI*, MagnesiumII*, MagnesiumIII*,
            
//             Manganese*, ManganeseI*, ManganeseII*, ManganeseIII*,
            
//             Mercury*, MercuryI*, MercuryII*, MercuryIII*,
            
//             Molybdenum*, MolybdenumI*, MolybdenumII*, MolybdenumIII*,
            
//             Neon*, NeonI*, NeonII*, NeonIII*,
            
//             Nickel*, NickelI*, NickelII*, NickelIII*,
            
//             Niobium*, NiobiumI*, NiobiumII*, NiobiumIII*,
            
//             Nitrogen*, NitrogenI*, NitrogenII*, NitrogenIII*,
            
//             Oxygen*, OxygenI*, OxygenII*, OxygenIII*,
            
//             Phosphorus*, PhosphorusI*, PhosphorusII*, PhosphorusIII*,
            
//             Platinum*, PlatinumI*, PlatinumII*, PlatinumIII*,
            
//             Potassium*, PotassiumI*, PotassiumII*, PotassiumIII*,
            
//             Radon*, RadonI*, RadonII*, RadonIII*,
            
//             Rhenium*, RheniumI*, RheniumII*, RheniumIII*,
            
//             Rubidium*, RubidiumI*, RubidiumII*, RubidiumIII*,
            
//             Silicon*, SiliconI*, SiliconII*, SiliconIII*,
            
//             Silver*, SilverI*, SilverII*, SilverIII*,
            
//             Sodium*, SodiumI*, SodiumII*, SodiumIII*,
            
//             Strontium*, StrontiumI*, StrontiumII*, StrontiumIII*,
            
//             Sulfur*, SulfurI*, SulfurII*, SulfurIII*,
            
//             Tantalum*, TantalumI*, TantalumII*, TantalumIII*,
            
//             Titanium*, TitaniumI*, TitaniumII*, TitaniumIII*,
            
//             Tungsten*, TungstenI*, TungstenII*, TungstenIII*,
            
//             Vanadium*, VanadiumI*, VanadiumII*, VanadiumIII*,
            
//             Xenon*, XenonI*, XenonII*, XenonIII*,
            
//             Zinc*, ZincI*, ZincII*, ZincIII*,
            
//             Zirconium*, ZirconiumI*, ZirconiumII*, ZirconiumIII*,
            
//             HydrogenAnion*, NitrogenAnion*, OxygenAnion*, 
            
//             MolecularHydrogen*, MolecularHydrogenI*, MolecularNitrogen*, 
            
//             MolecularNitrogenI*, MolecularOxygen*, MolecularOxygenAnion*, 
            
//             MolecularOxygenI*, */ /* Ozone*   */
        
//         >
//     ;
// ;

#endif
