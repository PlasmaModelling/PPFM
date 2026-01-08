 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "Species.h"

// ---------------- CarbonMonoxide ----------------

CarbonMonoxide::CarbonMonoxide() {
    costituents.push_back(new Carbon);
    costituents.push_back(new Oxygen);
    abundancy.push_back(1);
    abundancy.push_back(1);
}

int CarbonMonoxide::numberOfCostituents() { return 2; }
int CarbonMonoxide::operator()(int i) { return abundancy[i]; }
Element* CarbonMonoxide::operator[](int i) { return costituents[i]; }

double CarbonMonoxide::getMass() { return (28.0101) * amuKg; }
int CarbonMonoxide::getCharge() { return 0; }
std::string CarbonMonoxide::getFormula() { return "CO"; }
double CarbonMonoxide::formationEnergy() { return -1.7874e-18; } // J
double CarbonMonoxide::IonLim() { return 14.014 * eVtoJ; }

// ---------------- CarbonMonoxideI ----------------

CarbonMonoxideI::CarbonMonoxideI() {
    costituents.push_back(new Carbon);
    costituents.push_back(new Oxygen);
    abundancy.push_back(1);
    abundancy.push_back(1);
}

int CarbonMonoxideI::numberOfCostituents() { return 2; }
int CarbonMonoxideI::operator()(int i) { return abundancy[i]; }
Element* CarbonMonoxideI::operator[](int i) { return costituents[i]; }

double CarbonMonoxideI::getMass() { return (28.0101 - 1. * eamu) * amuKg; }
int CarbonMonoxideI::getCharge() { return 1; }
std::string CarbonMonoxideI::getFormula() { return "CO+"; }
double CarbonMonoxideI::formationEnergy() { return 4.5792e-19; } // J (efCO + IonLim(CO))
double CarbonMonoxideI::IonLim() { return 14.014 * eVtoJ; }

// ---------------- CarbonDioxide ----------------

CarbonDioxide::CarbonDioxide() {
    costituents.push_back(new Carbon);
    costituents.push_back(new Oxygen);
    abundancy.push_back(1);
    abundancy.push_back(2);
}

int CarbonDioxide::numberOfCostituents() { return 2; }
int CarbonDioxide::operator()(int i) { return abundancy[i]; }
Element* CarbonDioxide::operator[](int i) { return costituents[i]; }

double CarbonDioxide::getMass() { return (44.0095) * amuKg; }
int CarbonDioxide::getCharge() { return 0; }
std::string CarbonDioxide::getFormula() { return "CO2"; }
double CarbonDioxide::formationEnergy() { return -8.8369e-19; } // J
double CarbonDioxide::IonLim() { return 13.777 * eVtoJ; }

// ---------------- CarbonDioxideI ----------------

CarbonDioxideI::CarbonDioxideI() {
    costituents.push_back(new Carbon);
    costituents.push_back(new Oxygen);
    abundancy.push_back(1);
    abundancy.push_back(2);
}

int CarbonDioxideI::numberOfCostituents() { return 2; }
int CarbonDioxideI::operator()(int i) { return abundancy[i]; }
Element* CarbonDioxideI::operator[](int i) { return costituents[i]; }

double CarbonDioxideI::getMass() { return (44.0095 - eamu) * amuKg; }
int CarbonDioxideI::getCharge() { return 1; }
std::string CarbonDioxideI::getFormula() { return "CO2+"; }
double CarbonDioxideI::formationEnergy() { return 1.3236e-18; } // J
double CarbonDioxideI::IonLim() { return 9999999999999999; }

// ---------------- Ozone ----------------

Ozone::Ozone() {
    costituents.push_back(new Oxygen);
    abundancy.push_back(3);
}

int Ozone::numberOfCostituents() { return 3; }
int Ozone::operator()(int i) { return abundancy[0]; }
Element* Ozone::operator[](int i) { return costituents[0]; }

double Ozone::getMass() { return (47.9982) * amuKg; }
int Ozone::getCharge() { return 0; }
std::string Ozone::getFormula() { return "O3"; }
double Ozone::formationEnergy() { return (-1.478671) * eVtoJ; }
double Ozone::IonLim() { return 999999999999; }

MolecularHydrogen::MolecularHydrogen() {
    costituents[0] = new Hydrogen;
}

double MolecularHydrogen::getMass() { return 2.01588 * amuKg; }
std::string MolecularHydrogen::getFormula() { return "H2"; }
double MolecularHydrogen::formationEnergy() { return -7.174875e-19; }
double MolecularHydrogen::IonLim() { return 999999999999; }


MolecularNitrogen::MolecularNitrogen() {
    costituents[0] = new Nitrogen;
}

double MolecularNitrogen::getMass() { return 28.0134 * amuKg; }
std::string MolecularNitrogen::getFormula() { return "N2"; }
double MolecularNitrogen::formationEnergy() { return -1.563614e-18; }
double MolecularNitrogen::IonLim() { return 999999999999; }


MolecularOxygen::MolecularOxygen() {
    costituents[0] = new Oxygen;
}

double MolecularOxygen::getMass() { return 31.9988 * amuKg; }
std::string MolecularOxygen::getFormula() { return "O2"; }
double MolecularOxygen::formationEnergy() { return -5.11672987 * eVtoJ; }
double MolecularOxygen::IonLim() { return 12.0697 * eVtoJ; }


MolecularHydrogenI::MolecularHydrogenI() {
    costituents[0] = new Hydrogen;
}

double MolecularHydrogenI::getMass() { return (2.01588 - eamu) * amuKg; }
std::string MolecularHydrogenI::getFormula() { return "H2+"; }
double MolecularHydrogenI::formationEnergy() { return -7.174875e-19 + (15.42593 * eVtoJ); }
int MolecularHydrogenI::getCharge() { return 1; }
Element* MolecularHydrogenI::Constituent() { return new Hydrogen; }
double MolecularHydrogenI::IonLim() { return 999999999999; }


MolecularOxygenI::MolecularOxygenI() {
    costituents[0] = new Oxygen;
}

double MolecularOxygenI::getMass() { return (31.9988 - eamu) * amuKg; }
std::string MolecularOxygenI::getFormula() { return "O2+"; }
double MolecularOxygenI::formationEnergy() {
    return (new MolecularOxygen)->formationEnergy() +
           (new MolecularOxygen)->IonLim();
}
int MolecularOxygenI::getCharge() { return 1; }
Element* MolecularOxygenI::Constituent() { return new Oxygen; }
double MolecularOxygenI::IonLim() { return 999999999999; }


MolecularOxygenAnion::MolecularOxygenAnion() {
    costituents[0] = new Oxygen;
}

double MolecularOxygenAnion::getMass() { return (31.9988 + eamu) * amuKg; }
std::string MolecularOxygenAnion::getFormula() { return "O2-"; }
double MolecularOxygenAnion::formationEnergy() {
    return (new MolecularOxygen)->formationEnergy() +
           (-0.448 * eVtoJ);
}
int MolecularOxygenAnion::getCharge() { return -1; }
Element* MolecularOxygenAnion::Constituent() { return new Oxygen; }
double MolecularOxygenAnion::IonLim() { return 999999999999; }


MolecularNitrogenI::MolecularNitrogenI() {
    costituents[0] = new Nitrogen;
}

double MolecularNitrogenI::getMass() {
    return (new MolecularNitrogen)->getMass() - (eamu * amuKg);
}
std::string MolecularNitrogenI::getFormula() { return "N2+"; }
double MolecularNitrogenI::formationEnergy() {
    return (new MolecularNitrogen)->formationEnergy() +
           (15.581 * eVtoJ);
}
int MolecularNitrogenI::getCharge() { return 1; }
Element* MolecularNitrogenI::Constituent() { return new Nitrogen; }
double MolecularNitrogenI::IonLim() { return 999999999999; }
