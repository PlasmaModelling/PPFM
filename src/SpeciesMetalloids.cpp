#include"Species.h"

double Boron::getMass() { return 10.811 * amuKg ; }

std::string Boron::getFormula() { return "B" ; }

double Boron::formationEnergy() { return 0.0 * eVtoJ ; }

double Boron::IonLim() { return 8.298025 * eVtoJ ; }

double BoronI::getMass() { return (10.81045 - 1*eamu) * amuKg ; }

int BoronI::getCharge() { return 1 ; }

std::string BoronI::getFormula() { return "B+" ; }

double BoronI::formationEnergy() { return 8.298025 * eVtoJ ; }

double BoronI::IonLim() { return 25.1549 * eVtoJ ; }

Element* BoronI::Constituent() { return new Boron; }

double BoronII::getMass() { return (10.8099 - 2*eamu) * amuKg ; }

int BoronII::getCharge() { return 2 ; }

std::string BoronII::getFormula() { return "B+2" ; }

double BoronII::formationEnergy() { return 33.452925 * eVtoJ ; }

double BoronII::IonLim() { return 37.93 * eVtoJ ; }

Element* BoronII::Constituent() { return new Boron; }

double BoronIII::getMass() { return (10.80935 - 3*eamu) * amuKg ; }

int BoronIII::getCharge() { return 3 ; }

std::string BoronIII::getFormula() { return "B+3" ; }

double BoronIII::formationEnergy() { return 71.382925 * eVtoJ ; }

double BoronIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* BoronIII::Constituent() { return new Boron; }

double Silicon::getMass() { return 28.0855 * amuKg ; }

std::string Silicon::getFormula() { return "Si" ; }

double Silicon::formationEnergy() { return 0.0 * eVtoJ ; }

double Silicon::IonLim() { return 8.15168 * eVtoJ ; }

double SiliconI::getMass() { return (28.08495 - 1*eamu) * amuKg ; }

int SiliconI::getCharge() { return 1 ; }

std::string SiliconI::getFormula() { return "Si+" ; }

double SiliconI::formationEnergy() { return 8.15168 * eVtoJ ; }

double SiliconI::IonLim() { return 16.34585 * eVtoJ ; }

Element* SiliconI::Constituent() { return new Silicon; }

double SiliconII::getMass() { return (28.0844 - 2*eamu) * amuKg ; }

int SiliconII::getCharge() { return 2 ; }

std::string SiliconII::getFormula() { return "Si+2" ; }

double SiliconII::formationEnergy() { return 24.49753 * eVtoJ ; }

double SiliconII::IonLim() { return 33.49302 * eVtoJ ; }

Element* SiliconII::Constituent() { return new Silicon; }

double SiliconIII::getMass() { return (28.08385 - 3*eamu) * amuKg ; }

int SiliconIII::getCharge() { return 3 ; }

std::string SiliconIII::getFormula() { return "Si+3" ; }

double SiliconIII::formationEnergy() { return 57.99055 * eVtoJ ; }

double SiliconIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* SiliconIII::Constituent() { return new Silicon; }

double Arsenic::getMass() { return 74.9216 * amuKg ; }

std::string Arsenic::getFormula() { return "As" ; }

double Arsenic::formationEnergy() { return 0.0 * eVtoJ ; }

double Arsenic::IonLim() { return 9.7886 * eVtoJ ; }

double ArsenicI::getMass() { return (74.92105 - 1*eamu) * amuKg ; }

int ArsenicI::getCharge() { return 1 ; }

std::string ArsenicI::getFormula() { return "As+" ; }

double ArsenicI::formationEnergy() { return 9.7886 * eVtoJ ; }

double ArsenicI::IonLim() { return 18.509 * eVtoJ ; }

Element* ArsenicI::Constituent() { return new Arsenic; }

double ArsenicII::getMass() { return (74.9205 - 2*eamu) * amuKg ; }

int ArsenicII::getCharge() { return 2 ; }

std::string ArsenicII::getFormula() { return "As+2" ; }

double ArsenicII::formationEnergy() { return 28.2976 * eVtoJ ; }

double ArsenicII::IonLim() { return 29.6 * eVtoJ ; }

Element* ArsenicII::Constituent() { return new Arsenic; }

double ArsenicIII::getMass() { return (74.91995 - 3*eamu) * amuKg ; }

int ArsenicIII::getCharge() { return 3 ; }

std::string ArsenicIII::getFormula() { return "As+3" ; }

double ArsenicIII::formationEnergy() { return 57.8976 * eVtoJ ; }

double ArsenicIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* ArsenicIII::Constituent() { return new Arsenic; }

double Germanium::getMass() { return 72.63 * amuKg ; }

std::string Germanium::getFormula() { return "Ge" ; }

double Germanium::formationEnergy() { return 0.0 * eVtoJ ; }

double Germanium::IonLim() { return 7.89941 * eVtoJ ; }

double GermaniumI::getMass() { return (72.62945 - 1*eamu) * amuKg ; }

int GermaniumI::getCharge() { return 1 ; }

std::string GermaniumI::getFormula() { return "Ge+" ; }

double GermaniumI::formationEnergy() { return 7.89941 * eVtoJ ; }

double GermaniumI::IonLim() { return 15.932 * eVtoJ ; }

Element* GermaniumI::Constituent() { return new Germanium; }

double GermaniumII::getMass() { return (72.6289 - 2*eamu) * amuKg ; }

int GermaniumII::getCharge() { return 2 ; }

std::string GermaniumII::getFormula() { return "Ge+2" ; }

double GermaniumII::formationEnergy() { return 23.83141 * eVtoJ ; }

double GermaniumII::IonLim() { return 30.7 * eVtoJ ; }

Element* GermaniumII::Constituent() { return new Germanium; }

double GermaniumIII::getMass() { return (72.62835 - 3*eamu) * amuKg ; }

int GermaniumIII::getCharge() { return 3 ; }

std::string GermaniumIII::getFormula() { return "Ge+3" ; }

double GermaniumIII::formationEnergy() { return 54.53141 * eVtoJ ; }

double GermaniumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* GermaniumIII::Constituent() { return new Germanium; }

double Phosphorus::getMass() { return 30.973762 * amuKg ; }

std::string Phosphorus::getFormula() { return "P" ; }

double Phosphorus::formationEnergy() { return 0.0 * eVtoJ ; }

double Phosphorus::IonLim() { return 10.48669 * eVtoJ ; }

double PhosphorusI::getMass() { return (30.97321 - 1*eamu) * amuKg ; }

int PhosphorusI::getCharge() { return 1 ; }

std::string PhosphorusI::getFormula() { return "P+" ; }

double PhosphorusI::formationEnergy() { return 10.48669 * eVtoJ ; }

double PhosphorusI::IonLim() { return 19.72511 * eVtoJ ; }

Element* PhosphorusI::Constituent() { return new Phosphorus; }

double PhosphorusII::getMass() { return (30.97266 - 2*eamu) * amuKg ; }

int PhosphorusII::getCharge() { return 2 ; }

std::string PhosphorusII::getFormula() { return "P+2" ; }

double PhosphorusII::formationEnergy() { return 30.2118 * eVtoJ ; }

double PhosphorusII::IonLim() { return 30.61 * eVtoJ ; }

Element* PhosphorusII::Constituent() { return new Phosphorus; }

double PhosphorusIII::getMass() { return (30.97211 - 3*eamu) * amuKg ; }

int PhosphorusIII::getCharge() { return 3 ; }

std::string PhosphorusIII::getFormula() { return "P+3" ; }

double PhosphorusIII::formationEnergy() { return 60.8218 * eVtoJ ; }

double PhosphorusIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* PhosphorusIII::Constituent() { return new Phosphorus; }