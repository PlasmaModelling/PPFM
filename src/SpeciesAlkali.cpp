 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "Species.h"

double Lithium::getMass() { return 6.94108 * amuKg ; }

std::string Lithium::getFormula() { return "Li" ; }

double Lithium::formationEnergy() { return 0.0 * eVtoJ ; }

double Lithium::IonLim() { return 5.39172 * eVtoJ ; }

double LithiumI::getMass() { return (6.94053 - 1*eamu) * amuKg ; }

int LithiumI::getCharge() { return 1 ; }

std::string LithiumI::getFormula() { return "Li+" ; }

double LithiumI::formationEnergy() { return 5.39172 * eVtoJ ; }

double LithiumI::IonLim() { return 75.64 * eVtoJ ; }

Element* LithiumI::Constituent() { return new Lithium; }

double LithiumII::getMass() { return (6.94 - 2*eamu) * amuKg ; }

int LithiumII::getCharge() { return 2 ; }

std::string LithiumII::getFormula() { return "Li+2" ; }

double LithiumII::formationEnergy() { return 81.03172 * eVtoJ ; }

double LithiumII::IonLim() { return 122.45 * eVtoJ ; }

Element* LithiumII::Constituent() { return new Lithium; }

double LithiumIII::getMass() { return (6.93945 - 3*eamu) * amuKg ; }

int LithiumIII::getCharge() { return 3 ; }

std::string LithiumIII::getFormula() { return "Li+3" ; }

double LithiumIII::formationEnergy() { return 203.48172 * eVtoJ ; }

double LithiumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* LithiumIII::Constituent() { return new Lithium; }

double Sodium::getMass() { return 22.98976928 * amuKg ; }

std::string Sodium::getFormula() { return "Na" ; }

double Sodium::formationEnergy() { return 0.0 * eVtoJ ; }

double Sodium::IonLim() { return 5.139077 * eVtoJ ; }

double SodiumI::getMass() { return (22.989219 - 1*eamu) * amuKg ; }

int SodiumI::getCharge() { return 1 ; }

std::string SodiumI::getFormula() { return "Na+" ; }

double SodiumI::formationEnergy() { return 5.139077 * eVtoJ ; }

double SodiumI::IonLim() { return 47.28636 * eVtoJ ; }

Element* SodiumI::Constituent() { return new Sodium; }

double SodiumII::getMass() { return (22.98867 - 2*eamu) * amuKg ; }

int SodiumII::getCharge() { return 2 ; }

std::string SodiumII::getFormula() { return "Na+2" ; }

double SodiumII::formationEnergy() { return 52.425437 * eVtoJ ; }

double SodiumII::IonLim() { return 71.62 * eVtoJ ; }

Element* SodiumII::Constituent() { return new Sodium; }

double SodiumIII::getMass() { return (22.98812 - 3*eamu) * amuKg ; }

int SodiumIII::getCharge() { return 3 ; }

std::string SodiumIII::getFormula() { return "Na+3" ; }

double SodiumIII::formationEnergy() { return 124.045437 * eVtoJ ; }

double SodiumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* SodiumIII::Constituent() { return new Sodium; }

double Potassium::getMass() { return 39.0983 * amuKg ; }

std::string Potassium::getFormula() { return "K" ; }

double Potassium::formationEnergy() { return 0.0 * eVtoJ ; }

double Potassium::IonLim() { return 4.3407 * eVtoJ ; }

double PotassiumI::getMass() { return (39.09775 - 1*eamu) * amuKg ; }

int PotassiumI::getCharge() { return 1 ; }

std::string PotassiumI::getFormula() { return "K+" ; }

double PotassiumI::formationEnergy() { return 4.3407 * eVtoJ ; }

double PotassiumI::IonLim() { return 31.625 * eVtoJ ; }

Element* PotassiumI::Constituent() { return new Potassium; }

double PotassiumII::getMass() { return (39.0972 - 2*eamu) * amuKg ; }

int PotassiumII::getCharge() { return 2 ; }

std::string PotassiumII::getFormula() { return "K+2" ; }

double PotassiumII::formationEnergy() { return 35.9657 * eVtoJ ; }

double PotassiumII::IonLim() { return 45.8 * eVtoJ ; }

Element* PotassiumII::Constituent() { return new Potassium; }

double PotassiumIII::getMass() { return (39.09665 - 3*eamu) * amuKg ; }

int PotassiumIII::getCharge() { return 3 ; }

std::string PotassiumIII::getFormula() { return "K+3" ; }

double PotassiumIII::formationEnergy() { return 81.7657 * eVtoJ ; }

double PotassiumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* PotassiumIII::Constituent() { return new Potassium; }

double Rubidium::getMass() { return 85.4678 * amuKg ; }

std::string Rubidium::getFormula() { return "Rb" ; }

double Rubidium::formationEnergy() { return 0.0 * eVtoJ ; }

double Rubidium::IonLim() { return 4.1771 * eVtoJ ; }

double RubidiumI::getMass() { return (85.46725 - 1*eamu) * amuKg ; }

int RubidiumI::getCharge() { return 1 ; }

std::string RubidiumI::getFormula() { return "Rb+" ; }

double RubidiumI::formationEnergy() { return 4.1771 * eVtoJ ; }

double RubidiumI::IonLim() { return 27.856 * eVtoJ ; }

Element* RubidiumI::Constituent() { return new Rubidium; }

double RubidiumII::getMass() { return (85.4667 - 2*eamu) * amuKg ; }

int RubidiumII::getCharge() { return 2 ; }

std::string RubidiumII::getFormula() { return "Rb+2" ; }

double RubidiumII::formationEnergy() { return 32.0331 * eVtoJ ; }

double RubidiumII::IonLim() { return 40.9 * eVtoJ ; }

Element* RubidiumII::Constituent() { return new Rubidium; }

double RubidiumIII::getMass() { return (85.46615 - 3*eamu) * amuKg ; }

int RubidiumIII::getCharge() { return 3 ; }

std::string RubidiumIII::getFormula() { return "Rb+3" ; }

double RubidiumIII::formationEnergy() { return 72.9331 * eVtoJ ; }

double RubidiumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* RubidiumIII::Constituent() { return new Rubidium; }

double Cesium::getMass() { return 132.90545 * amuKg ; }

std::string Cesium::getFormula() { return "Cs" ; }

double Cesium::formationEnergy() { return 0.0 * eVtoJ ; }

double Cesium::IonLim() { return 3.8939 * eVtoJ ; }

double CesiumI::getMass() { return (132.9049 - 1*eamu) * amuKg ; }

int CesiumI::getCharge() { return 1 ; }

std::string CesiumI::getFormula() { return "Cs+" ; }

double CesiumI::formationEnergy() { return 3.8939 * eVtoJ ; }

double CesiumI::IonLim() { return 25.02 * eVtoJ ; }

Element* CesiumI::Constituent() { return new Cesium; }

double CesiumII::getMass() { return (132.90435 - 2*eamu) * amuKg ; }

int CesiumII::getCharge() { return 2 ; }

std::string CesiumII::getFormula() { return "Cs+2" ; }

double CesiumII::formationEnergy() { return 28.9139 * eVtoJ ; }

double CesiumII::IonLim() { return 34.21 * eVtoJ ; }

Element* CesiumII::Constituent() { return new Cesium; }

double CesiumIII::getMass() { return (132.9038 - 3*eamu) * amuKg ; }

int CesiumIII::getCharge() { return 3 ; }

std::string CesiumIII::getFormula() { return "Cs+3" ; }

double CesiumIII::formationEnergy() { return 63.1239 * eVtoJ ; }

double CesiumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* CesiumIII::Constituent() { return new Cesium; }