 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "Species.h"

// ---------------- CarbonMonoxide ----------------

CarbonMonoxide::CarbonMonoxide() {
    costituents[0] = (new Carbon);
    costituents[1] = (new Oxygen);
    abundancy[0] = (1);
    abundancy[1] = (1);
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
    costituents[0] = (new Carbon);
    costituents[1] = (new Oxygen);
    abundancy[0] = (1);
    abundancy[1] = (1);
}

int CarbonMonoxideI::numberOfCostituents() { return 2; }
int CarbonMonoxideI::operator()(int i) { return abundancy[i]; }
Element* CarbonMonoxideI::operator[](int i) { return costituents[i]; }

double CarbonMonoxideI::getMass() { return (28.0101 - 1. * eamu) * amuKg; }
int CarbonMonoxideI::getCharge() { return 1; }
std::string CarbonMonoxideI::getFormula() { return "CO+"; }
double CarbonMonoxideI::formationEnergy() { return 4.5792e-19; } // J (efCO + IonLim(CO))
double CarbonMonoxideI::IonLim() { return 14.014 * eVtoJ; }
Element* CarbonMonoxideI::Constituent() { std::cerr<<"Asked for constituent of an HeteroNuclearMolecule. aborting...\n"; abort(); return nullptr; }

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
Element* CarbonDioxideI::Constituent() { std::cerr<<"Asked for constituent of an PolyAtomicMolecule. aborting...\n"; abort(); return nullptr; }

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

/* 
 Dummy implementation of molecules missing from POM mixture 
 C C+ C+2 C- C2 C2+ C2H C2H2 C2H4 C2O C3 CH CH2 CH2O CH3 CH4 
 CHO CHO+ CO CO+ CO2 CO2+ H H+ H- H2 H2+ H2O HO2 
 O O+ O+2 O- O2 O2+ OH OH+ e-.
 implementation is called dummy for the missing of formation energies 
 which are not required for collision integrals calculations neither 
 if the chemical composition is loaded from external files.
 */

 // -------------------- CarbonAnion --------------------
double CarbonAnion::getMass() { return ((new Carbon)->getMass() + 1. * eamu); }
int CarbonAnion::getCharge() { return -1; }
std::string CarbonAnion::getFormula() { return "C-"; }
double CarbonAnion::formationEnergy() { return 0.0; }
double CarbonAnion::IonLim() { return 0.0; }
Element* CarbonAnion::Constituent() { return new Carbon(); }

// -------------------- DiCarbon --------------------
DiCarbon::DiCarbon() {
    costituents[0] = new Carbon();
}

double DiCarbon::getMass() { return 2. * (new Carbon)->getMass(); }
std::string DiCarbon::getFormula() { return "C2"; }
double DiCarbon::formationEnergy() { return 0.0; }
double DiCarbon::IonLim() { return 0.0; }

// -------------------- DiCarbonI --------------------
DiCarbonI::DiCarbonI() {
    costituents[0] = new Carbon();
}

double DiCarbonI::getMass() { return 2. * (new Carbon)->getMass(); }
int DiCarbonI::getCharge() { return 1; }
std::string DiCarbonI::getFormula() { return "C2+"; }
double DiCarbonI::formationEnergy() { return 0.0; }
double DiCarbonI::IonLim() { return 0.0; }
Element* DiCarbonI::Constituent() { return new Carbon(); }

// -------------------- EthynilRadical --------------------
EthynilRadical::EthynilRadical() {
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Hydrogen();
    abundancy[0] = 2;
    abundancy[1] = 1;
}

int EthynilRadical::numberOfCostituents() { return 2; }
double EthynilRadical::getMass() { return (2. * (new Carbon)->getMass()) + (new Hydrogen)->getMass(); }
std::string EthynilRadical::getFormula() { return "C2H"; }
double EthynilRadical::formationEnergy() { return 0.0; }
double EthynilRadical::IonLim() { return 0.0; }

// -------------------- Acetylene --------------------
Acetylene::Acetylene() {
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Hydrogen();
    abundancy[0] = 2;
    abundancy[1] = 2;
}

int Acetylene::numberOfCostituents() { return 2; }
double Acetylene::getMass() { return (2. * (new Carbon)->getMass()) + (2. * (new Hydrogen)->getMass()); }
std::string Acetylene::getFormula() { return "C2H2"; }
double Acetylene::formationEnergy() { return 0.0; }
double Acetylene::IonLim() { return 0.0; }

// -------------------- Ethylene --------------------
Ethylene::Ethylene() {
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Hydrogen();
    abundancy[0] = 2;
    abundancy[1] = 4;
}

int Ethylene::numberOfCostituents() { return 2; }
double Ethylene::getMass() { return (2. * (new Carbon)->getMass()) + (4. * (new Hydrogen)->getMass()); }
std::string Ethylene::getFormula() { return "C2H4"; }
double Ethylene::formationEnergy() { return 0.0; }
double Ethylene::IonLim() { return 0.0; }

// -------------------- DiCarbonMonoxide --------------------
DiCarbonMonoxide::DiCarbonMonoxide() {
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Oxygen();
    abundancy[0] = 2;
    abundancy[1] = 1;
}

int DiCarbonMonoxide::numberOfCostituents() { return 2; }
double DiCarbonMonoxide::getMass() { return (2. * (new Carbon)->getMass()) + (new Oxygen)->getMass(); }
std::string DiCarbonMonoxide::getFormula() { return "C2O"; }
double DiCarbonMonoxide::formationEnergy() { return 0.0; }
double DiCarbonMonoxide::IonLim() { return 0.0; }

// -------------------- TriCarbon --------------------
TriCarbon::TriCarbon() {
    costituents.resize(1);
    abundancy.resize(1);
    costituents[0] = new Carbon();
    abundancy[0] = 3;
}

int TriCarbon::numberOfCostituents() { return 1; }
double TriCarbon::getMass() { return 3. * (new Carbon)->getMass(); }
std::string TriCarbon::getFormula() { return "C3"; }
double TriCarbon::formationEnergy() { return 0.0; }
double TriCarbon::IonLim() { return 0.0; }

// -------------------- MethylidineRadical --------------------
MethylidineRadical::MethylidineRadical() {
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Hydrogen();
    abundancy[0] = 1;
    abundancy[1] = 1;
}

int MethylidineRadical::numberOfCostituents() { return 2; }
double MethylidineRadical::getMass() { return (new Carbon)->getMass() + (new Hydrogen)->getMass(); }
std::string MethylidineRadical::getFormula() { return "CH"; }
double MethylidineRadical::formationEnergy() { return 0.0; }
double MethylidineRadical::IonLim() { return 0.0; }

// -------------------- Methylene --------------------
Methylene::Methylene() {
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Hydrogen();
    abundancy[0] = 1;
    abundancy[1] = 2;
}

int Methylene::numberOfCostituents() { return 2; }
double Methylene::getMass() { return (new Carbon)->getMass() + 2. * (new Hydrogen)->getMass(); }
std::string Methylene::getFormula() { return "CH2"; }
double Methylene::formationEnergy() { return 0.0; }
double Methylene::IonLim() { return 0.0; }

// -------------------- MethylRadical --------------------
MethylRadical::MethylRadical() {
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Hydrogen();
    abundancy[0] = 1;
    abundancy[1] = 3;
}

int MethylRadical::numberOfCostituents() { return 2; }
double MethylRadical::getMass() { return (new Carbon)->getMass() + 3. * (new Hydrogen)->getMass(); }
std::string MethylRadical::getFormula() { return "CH3"; }
double MethylRadical::formationEnergy() { return 0.0; }
double MethylRadical::IonLim() { return 0.0; }

// -------------------- Methane --------------------
Methane::Methane() {
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Hydrogen();
    abundancy[0] = 1;
    abundancy[1] = 4;
}

int Methane::numberOfCostituents() { return 2; }
double Methane::getMass() { return (new Carbon)->getMass() + 4. * (new Hydrogen)->getMass(); }
std::string Methane::getFormula() { return "CH4"; }
double Methane::formationEnergy() { return 0.0; }
double Methane::IonLim() { return 0.0; }

// -------------------- FormylRadical --------------------
FormylRadical::FormylRadical() {
    costituents.resize(3);
    abundancy.resize(3);
    costituents[0] = new Carbon();
    costituents[1] = new Hydrogen();
    costituents[2] = new Oxygen();
    abundancy[0] = 1;
    abundancy[1] = 1;
    abundancy[2] = 1;
}

int FormylRadical::numberOfCostituents() { return 3; }
double FormylRadical::getMass() {
    return (new Carbon)->getMass() + (new Hydrogen)->getMass() + (new Oxygen)->getMass();
}
std::string FormylRadical::getFormula() { return "CHO"; }
double FormylRadical::formationEnergy() { return 0.0; }
double FormylRadical::IonLim() { return 0.0; }

// -------------------- FormylRadicalI --------------------
FormylRadicalI::FormylRadicalI() {
    costituents.resize(3);
    abundancy.resize(3);
    costituents[0] = new Carbon();
    costituents[1] = new Hydrogen();
    costituents[2] = new Oxygen();
    abundancy[0] = 1;
    abundancy[1] = 1;
    abundancy[2] = 1;
}

int FormylRadicalI::numberOfCostituents() { return 3; }
double FormylRadicalI::getMass() {
    return (new Carbon)->getMass() + (new Hydrogen)->getMass() + (new Oxygen)->getMass();
}
int FormylRadicalI::getCharge() { return 1; }
std::string FormylRadicalI::getFormula() { return "CHO+"; }
double FormylRadicalI::formationEnergy() { return 0.0; }
double FormylRadicalI::IonLim() { return 0.0; }
Element* FormylRadicalI::Constituent() { return new Carbon(); }

// -------------------- Formaldehyde --------------------
Formaldehyde::Formaldehyde() {
    costituents.resize(3);
    abundancy.resize(3);
    costituents[0] = new Carbon();
    costituents[1] = new Hydrogen();
    costituents[2] = new Oxygen();
    abundancy[0] = 1;
    abundancy[1] = 2;
    abundancy[2] = 1;
}

int Formaldehyde::numberOfCostituents() { return 3; }
double Formaldehyde::getMass() {
    return (new Carbon)->getMass() + (2. * (new Hydrogen)->getMass()) + (new Oxygen)->getMass();
}
std::string Formaldehyde::getFormula() { return "CH2O"; }
double Formaldehyde::formationEnergy() { return 0.0; }
double Formaldehyde::IonLim() { return 0.0; }

// -------------------- Water --------------------
Water::Water() {
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Hydrogen();
    costituents[1] = new Oxygen();
    abundancy[0] = 2;
    abundancy[1] = 1;
}

int Water::numberOfCostituents() { return 2; }
double Water::getMass() { return (2. * (new Hydrogen)->getMass()) + (new Oxygen)->getMass(); }
std::string Water::getFormula() { return "H2O"; }
double Water::formationEnergy() { return 0.0; }
double Water::IonLim() { return 0.0; }

// -------------------- HydroPeroxyRadical --------------------
HydroPeroxyRadical::HydroPeroxyRadical() {
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Hydrogen();
    costituents[1] = new Oxygen();
    abundancy[0] = 1;
    abundancy[1] = 2;
}

int HydroPeroxyRadical::numberOfCostituents() { return 2; }
double HydroPeroxyRadical::getMass() { return (new Hydrogen)->getMass() + (2. * (new Oxygen)->getMass()); }
std::string HydroPeroxyRadical::getFormula() { return "HO2"; }
double HydroPeroxyRadical::formationEnergy() { return 0.0; }
double HydroPeroxyRadical::IonLim() { return 0.0; }

// -------------------- HydroxilRadical --------------------
HydroxilRadical::HydroxilRadical() {
    costituents[0] = new Hydrogen();
    costituents[1] = new Oxygen();
}

int HydroxilRadical::numberOfCostituents() { return 2; }
double HydroxilRadical::getMass() { return (new Hydrogen)->getMass() + (new Oxygen)->getMass(); }
std::string HydroxilRadical::getFormula() { return "OH"; }
double HydroxilRadical::formationEnergy() { return 0.0; }
double HydroxilRadical::IonLim() { return 0.0; }

// -------------------- HydroxilRadicalI --------------------
HydroxilRadicalI::HydroxilRadicalI() {
    costituents[0] = new Hydrogen();
    costituents[1] = new Oxygen();
}

int HydroxilRadicalI::numberOfCostituents() { return 2; }
double HydroxilRadicalI::getMass() { return (new Hydrogen)->getMass() + (new Oxygen)->getMass(); }
int HydroxilRadicalI::getCharge() { return 1; }
std::string HydroxilRadicalI::getFormula() { return "OH+"; }
double HydroxilRadicalI::formationEnergy() { return 0.0; }
double HydroxilRadicalI::IonLim() { return 0.0; }

Element* HydroxilRadicalI::Constituent() {
    std::cerr << "Error: Constituent() not implemented for HydroxilRadicalI." << std::endl;
    return nullptr;
}