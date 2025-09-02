 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include"Species.h"

double Carbon::getMass() { return 12.0107 * amuKg ; }

std::string Carbon::getFormula() { return "C" ; }

double Carbon::formationEnergy() { return 0.0 * eVtoJ ; }

double Carbon::IonLim() { return 11.2603 * eVtoJ ; }

double CarbonI::getMass() { return (12.01015 - 1*eamu) * amuKg ; }

int CarbonI::getCharge() { return 1 ; }

std::string CarbonI::getFormula() { return "C+" ; }

double CarbonI::formationEnergy() { return 11.2603 * eVtoJ ; }

double CarbonI::IonLim() { return 24.3833 * eVtoJ ; }

Element* CarbonI::Constituent() { return new Carbon; }

double CarbonII::getMass() { return (12.0096 - 2*eamu) * amuKg ; }

int CarbonII::getCharge() { return 2 ; }

std::string CarbonII::getFormula() { return "C+2" ; }

double CarbonII::formationEnergy() { return 35.6436 * eVtoJ ; }

double CarbonII::IonLim() { return 47.887 * eVtoJ ; }

Element* CarbonII::Constituent() { return new Carbon; }

double CarbonIII::getMass() { return (12.00905 - 3*eamu) * amuKg ; }

int CarbonIII::getCharge() { return 3 ; }

std::string CarbonIII::getFormula() { return "C+3" ; }

double CarbonIII::formationEnergy() { return 83.5306 * eVtoJ ; }

double CarbonIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* CarbonIII::Constituent() { return new Carbon; }

double Bromine::getMass() { return 79.904 * amuKg ; }

std::string Bromine::getFormula() { return "Br" ; }

double Bromine::formationEnergy() { return 0.0 * eVtoJ ; }

double Bromine::IonLim() { return 11.81381 * eVtoJ ; }

double BromineI::getMass() { return (79.90344 - 1*eamu) * amuKg ; }

int BromineI::getCharge() { return 1 ; }

std::string BromineI::getFormula() { return "Br+" ; }

double BromineI::formationEnergy() { return 11.81381 * eVtoJ ; }

double BromineI::IonLim() { return 21.8 * eVtoJ ; }

Element* BromineI::Constituent() { return new Bromine; }

double BromineII::getMass() { return (79.90288 - 2*eamu) * amuKg ; }

int BromineII::getCharge() { return 2 ; }

std::string BromineII::getFormula() { return "Br+2" ; }

double BromineII::formationEnergy() { return 33.61381 * eVtoJ ; }

double BromineII::IonLim() { return 36.0 * eVtoJ ; }

Element* BromineII::Constituent() { return new Bromine; }

double BromineIII::getMass() { return (79.90232 - 3*eamu) * amuKg ; }

int BromineIII::getCharge() { return 3 ; }

std::string BromineIII::getFormula() { return "Br+3" ; }

double BromineIII::formationEnergy() { return 70.61381 * eVtoJ ; }

double BromineIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* BromineIII::Constituent() { return new Bromine; }

double Iodine::getMass() { return 126.90447 * amuKg ; }

std::string Iodine::getFormula() { return "I" ; }

double Iodine::formationEnergy() { return 0.0 * eVtoJ ; }

double Iodine::IonLim() { return 10.45126 * eVtoJ ; }

double IodineI::getMass() { return (126.90392 - 1*eamu) * amuKg ; }

int IodineI::getCharge() { return 1 ; }

std::string IodineI::getFormula() { return "I+" ; }

double IodineI::formationEnergy() { return 10.45126 * eVtoJ ; }

double IodineI::IonLim() { return 19.10995 * eVtoJ ; }

Element* IodineI::Constituent() { return new Iodine; }

double IodineII::getMass() { return (126.90337 - 2*eamu) * amuKg ; }

int IodineII::getCharge() { return 2 ; }

std::string IodineII::getFormula() { return "I+2" ; }

double IodineII::formationEnergy() { return 29.56121 * eVtoJ ; }

double IodineII::IonLim() { return 35.0 * eVtoJ ; }

Element* IodineII::Constituent() { return new Iodine; }

double IodineIII::getMass() { return (126.90282 - 3*eamu) * amuKg ; }

int IodineIII::getCharge() { return 3 ; }

std::string IodineIII::getFormula() { return "I+3" ; }

double IodineIII::formationEnergy() { return 64.56121 * eVtoJ ; }

double IodineIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* IodineIII::Constituent() { return new Iodine; }

double Sulfur::getMass() { return 32.065 * amuKg ; }

std::string Sulfur::getFormula() { return "S" ; }

double Sulfur::formationEnergy() { return 0.0 * eVtoJ ; }

double Sulfur::IonLim() { return 10.36001 * eVtoJ ; }

double SulfurI::getMass() { return (32.06445 - 1*eamu) * amuKg ; }

int SulfurI::getCharge() { return 1 ; }

std::string SulfurI::getFormula() { return "S+" ; }

double SulfurI::formationEnergy() { return 10.36001 * eVtoJ ; }

double SulfurI::IonLim() { return 23.33064 * eVtoJ ; }

Element* SulfurI::Constituent() { return new Sulfur; }

double SulfurII::getMass() { return (32.0639 - 2*eamu) * amuKg ; }

int SulfurII::getCharge() { return 2 ; }

std::string SulfurII::getFormula() { return "S+2" ; }

double SulfurII::formationEnergy() { return 33.69065 * eVtoJ ; }

double SulfurII::IonLim() { return 34.83 * eVtoJ ; }

Element* SulfurII::Constituent() { return new Sulfur; }

double SulfurIII::getMass() { return (32.06335 - 3*eamu) * amuKg ; }

int SulfurIII::getCharge() { return 3 ; }

std::string SulfurIII::getFormula() { return "S+3" ; }

double SulfurIII::formationEnergy() { return 68.52065 * eVtoJ ; }

double SulfurIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* SulfurIII::Constituent() { return new Sulfur; }