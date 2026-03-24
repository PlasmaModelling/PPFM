 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include "Species.h"

/* VERY IMPORTANT */
/* 
** The molecular formation energies used in the Saha equations, i.e. 
** "dissociation energies with negative sign", refer to atomization 
** reactions (ex:CO2 -> C + 2O) and at T = 0K.
** Please, when calculating dissociation energies refer to atomization
** reactions and apply correction for H-H0 that can be found for example 
** in JANAF tables. 
*/

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
double CarbonMonoxide::formationEnergy() { return -1.7797e-18; } // J
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
double CarbonMonoxideI::formationEnergy() { return 4.6554e-19; } // J (efCO + IonLim(CO))
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
double CarbonDioxide::formationEnergy() { return -2.6534e-18; } // J
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
double CarbonDioxideI::formationEnergy() { return -4.461e-19; } // J
double CarbonDioxideI::IonLim() { return 9999999999999999; }
Element* CarbonDioxideI::Constituent() { std::cerr<<"Asked for constituent of an PolyAtomicMolecule. aborting...\n"; abort(); return nullptr; }

// ---------------- CarbonDioxideAnion ----------------

CarbonDioxideAnion::CarbonDioxideAnion() {
    costituents.push_back(new Carbon);
    costituents.push_back(new Oxygen);
    abundancy.push_back(1);
    abundancy.push_back(2);
}

int CarbonDioxideAnion::numberOfCostituents() { return 2; }
int CarbonDioxideAnion::operator()(int i) { return abundancy[i]; }
Element* CarbonDioxideAnion::operator[](int i) { return costituents[i]; }

double CarbonDioxideAnion::getMass() { return (44.0095 + eamu) * amuKg; }
int CarbonDioxideAnion::getCharge() { return -1; }
std::string CarbonDioxideAnion::getFormula() { return "CO2-"; }
double CarbonDioxideAnion::formationEnergy() { return (new CarbonDioxide)->formationEnergy()+(0.6*eVtoJ); } // J
double CarbonDioxideAnion::IonLim() { return 9999999999999999; }
Element* CarbonDioxideAnion::Constituent() { std::cerr<<"Asked for constituent of an PolyAtomicMolecule. aborting...\n"; abort(); return nullptr; }

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
double Ozone::formationEnergy() { return -9.5368e-19; }
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
double MolecularOxygen::formationEnergy() { return -8.1963e-19; }
double MolecularOxygen::IonLim() { return 12.0697 * eVtoJ; }


MolecularHydrogenI::MolecularHydrogenI() {
    costituents[0] = new Hydrogen;
}

double MolecularHydrogenI::getMass() { return (2.01588 - eamu) * amuKg; }
std::string MolecularHydrogenI::getFormula() { return "H2+"; }
double MolecularHydrogenI::formationEnergy() { return ((new MolecularHydrogen)->formationEnergy()) + (15.42593 * eVtoJ); }
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

// -------------------- DiCarbon --------------------
DiCarbon::DiCarbon() {
    costituents[0] = new Carbon();
}

double DiCarbon::getMass() { return 2. * (new Carbon)->getMass(); }
std::string DiCarbon::getFormula() { return "C2"; }
double DiCarbon::formationEnergy() { return -9.8488e-19; }
double DiCarbon::IonLim() { return 999999999; }

// -------------------- DiCarbonI --------------------
DiCarbonI::DiCarbonI() {
    costituents[0] = new Carbon();
}

double DiCarbonI::getMass() { return 2. * (new Carbon)->getMass(); }
int DiCarbonI::getCharge() { return 1; }
std::string DiCarbonI::getFormula() { return "C2+"; }
double DiCarbonI::formationEnergy() { return 8.416e-19; }
double DiCarbonI::IonLim() { return 999999999; }
Element* DiCarbonI::Constituent() { return new Carbon(); }

// -------------------- DiCarbonAnion --------------------
DiCarbonAnion::DiCarbonAnion() {
    costituents[0] = new Carbon();
}

double DiCarbonAnion::getMass() { return 2. * (new Carbon)->getMass(); }
int DiCarbonAnion::getCharge() { return -1; }
std::string DiCarbonAnion::getFormula() { return "C2-"; }
double DiCarbonAnion::formationEnergy() { return 8.416e-19 - (3.273*eVtoJ); }
double DiCarbonAnion::IonLim() { return 999999999; }
Element* DiCarbonAnion::Constituent() { return new Carbon(); }

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
double EthynilRadical::formationEnergy() { return -1.9337e-18; }
double EthynilRadical::IonLim() { return 999999999; }

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
double Acetylene::formationEnergy() { return -2.7019e-18; }
double Acetylene::IonLim() { return 999999999; }

// -------------------- PropadieneDione --------------------
PropadieneDione::PropadieneDione(){
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Oxygen();
    abundancy[0] = 3 ; 
    abundancy[1] = 2 ;
}

int PropadieneDione::numberOfCostituents() { return 2; }
double PropadieneDione::getMass() { return 3.*(new Carbon)->getMass() + 2.*(new Oxygen)->getMass(); };
std::string PropadieneDione::getFormula() { return "C3O2"; };
double PropadieneDione::formationEnergy() { return -3.7926e-18; };
double PropadieneDione::IonLim() { return 99999999999.; };

// -------------------- PropadieneDioneI --------------------
PropadieneDioneI::PropadieneDioneI() : PropadieneDione() {}

std::string PropadieneDioneI::getFormula() { return "C3O2+"; };
double PropadieneDioneI::formationEnergy() { return -2.0935e-18; };

Element* PropadieneDioneI::Constituent() {
    std::cerr << "Error: Constituent() not implemented for HydroxilRadicalI." << std::endl;
    return nullptr;
}

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
double Ethylene::formationEnergy() { return -3.6956e-18; }
double Ethylene::IonLim() { return 999999999; }

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
double DiCarbonMonoxide::formationEnergy() { return -2.3025e-18; }
double DiCarbonMonoxide::IonLim() { return 999999999; }

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
double TriCarbon::formationEnergy() { return -2.1954e-18; }
double TriCarbon::IonLim() { return 999999999; }

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
double MethylidineRadical::formationEnergy() { return -5.5883e-19; }
double MethylidineRadical::IonLim() { return 999999999; }

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
double Methylene::formationEnergy() { return -1.2576e-18; }
double Methylene::IonLim() { return 999999999; }

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
double MethylRadical::formationEnergy() { return -2.0097e-18; }
double MethylRadical::IonLim() { return 999999999; }

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
double Methane::formationEnergy() { return -2.727e-18; }
double Methane::IonLim() { return 999999999; }

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
double FormylRadical::formationEnergy() { return -1.8778e-18; }
double FormylRadical::IonLim() { return 999999999; }

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
double FormylRadicalI::formationEnergy() { return -5.7687e-19; }
double FormylRadicalI::IonLim() { return 999999999; }
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
double Formaldehyde::formationEnergy() { return -2.4943e-18; }
double Formaldehyde::IonLim() { return 999999999; }

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
double Water::formationEnergy() { return -1.5784e-18; }
double Water::IonLim() { return 999999999; }

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
double HydroPeroxyRadical::formationEnergy() { return -1.17e-18; }
double HydroPeroxyRadical::IonLim() { return 999999999; }

// -------------------- HydroxilRadical --------------------
HydroxilRadical::HydroxilRadical() {
    costituents[0] = new Hydrogen();
    costituents[1] = new Oxygen();
}

int HydroxilRadical::numberOfCostituents() { return 2; }
double HydroxilRadical::getMass() { return (new Hydrogen)->getMass() + (new Oxygen)->getMass(); }
std::string HydroxilRadical::getFormula() { return "OH"; }
double HydroxilRadical::formationEnergy() { return -7.0479e-19; }
double HydroxilRadical::IonLim() { return 999999999; }

// -------------------- HydroxilRadicalI --------------------
HydroxilRadicalI::HydroxilRadicalI() {
    costituents[0] = new Hydrogen();
    costituents[1] = new Oxygen();
}

int HydroxilRadicalI::numberOfCostituents() { return 2; }
double HydroxilRadicalI::getMass() { return (new Hydrogen)->getMass() + (new Oxygen)->getMass(); }
int HydroxilRadicalI::getCharge() { return 1; }
std::string HydroxilRadicalI::getFormula() { return "OH+"; }
double HydroxilRadicalI::formationEnergy() { return 1.3808e-18; }
double HydroxilRadicalI::IonLim() { return 999999999; }

Element* HydroxilRadicalI::Constituent() {
    std::cerr << "Error: Constituent() not implemented for HydroxilRadicalI." << std::endl;
    return nullptr;
}

// -------------------- TetraCarbon --------------------

TetraCarbon::TetraCarbon() {
    costituents.resize(1);
    abundancy.resize(1);
    costituents[0] = new Carbon();
    abundancy[0] = 4;
}

int TetraCarbon::numberOfCostituents() { return 1; }
double TetraCarbon::getMass() { return 4. * (new Carbon)->getMass(); }
std::string TetraCarbon::getFormula() { return "C4"; }
double TetraCarbon::formationEnergy() { return -3.1229e-18; }
double TetraCarbon::IonLim() { return 99999999999; }


// -------------------- PentaCarbon --------------------

PentaCarbon::PentaCarbon() {
    costituents.resize(1);
    abundancy.resize(1);
    costituents[0] = new Carbon();
    abundancy[0] = 5;
}

int PentaCarbon::numberOfCostituents() { return 1; }
double PentaCarbon::getMass() { return 5. * (new Carbon)->getMass(); }
std::string PentaCarbon::getFormula() { return "C5"; }
double PentaCarbon::formationEnergy() { return -4.29e-18; }
double PentaCarbon::IonLim() { return 99999999999; }


// -------------------- EsaCarbon --------------------

EsaCarbon::EsaCarbon() {
    costituents.resize(1);
    abundancy.resize(1);
    costituents[0] = new Carbon();
    abundancy[0] = 6;
}

int EsaCarbon::numberOfCostituents() { return 1; }
double EsaCarbon::getMass() { return 6. * (new Carbon)->getMass(); }
std::string EsaCarbon::getFormula() { return "C6"; }
double EsaCarbon::formationEnergy() { return -5.1229e-18; }
double EsaCarbon::IonLim() { return 99999999999; }

// -------------------- CarbonTrioxide --------------------

CarbonTrioxide::CarbonTrioxide() {
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Oxygen();
    abundancy[0] = 1;
    abundancy[1] = 3;
}

int CarbonTrioxide::numberOfCostituents() { return 2; }
double CarbonTrioxide::getMass() { return (new Carbon)->getMass() + 3.*(new Oxygen)->getMass(); }
std::string CarbonTrioxide::getFormula() { return "CO3"; }
double CarbonTrioxide::formationEnergy() { return -2.2042e-18; }
double CarbonTrioxide::IonLim() { return 99999999999; }

// -------------------- CarbonTrioxideAnion --------------------

CarbonTrioxideAnion::CarbonTrioxideAnion() : CarbonTrioxide() {}

int CarbonTrioxideAnion::getCharge(){return -1;}
std::string CarbonTrioxideAnion::getFormula() { return "CO3-"; }
double CarbonTrioxideAnion::formationEnergy() { return -2.797e-18; }

Element* CarbonTrioxideAnion::Constituent() {
    std::cerr << "Error: Constituent() not implemented for HydroxilRadicalI." << std::endl;
    return nullptr;
}

CNCRadical::CNCRadical(){
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Nitrogen();
    abundancy[0] = 2;
    abundancy[1] = 1;
}

int CNCRadical::numberOfCostituents(){
    return 2;
}

double CNCRadical::getMass(){
    return costituents[0]->getMass() * abundancy[0] + costituents[1]->getMass() * abundancy[1];
}

std::string CNCRadical::getFormula(){
    return "C2N";
}

double CNCRadical::formationEnergy(){
    return -1.8228e-18;
}

double CNCRadical::IonLim(){
    return 9999999999;
}



Cyanogen::Cyanogen(){
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Nitrogen();
    abundancy[0] = 2;
    abundancy[1] = 2;
}

int Cyanogen::numberOfCostituents(){
    return 2;
}

double Cyanogen::getMass(){
    return costituents[0]->getMass() * abundancy[0] + costituents[1]->getMass() * abundancy[1];
}

std::string Cyanogen::getFormula(){
    return "C2N2";
}

double Cyanogen::formationEnergy(){
    return -3.4154e-18;
}

double Cyanogen::IonLim(){
    return 9999999999;
}



CyanoEthynilRadical::CyanoEthynilRadical(){
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Nitrogen();
    abundancy[0] = 3;
    abundancy[1] = 1;
}

int CyanoEthynilRadical::numberOfCostituents(){
    return 2;
}

double CyanoEthynilRadical::getMass(){
    return costituents[0]->getMass() * abundancy[0] + costituents[1]->getMass() * abundancy[1];
}

std::string CyanoEthynilRadical::getFormula(){
    return "C3N";
}

double CyanoEthynilRadical::formationEnergy(){
    return -3.113e-18;
}

double CyanoEthynilRadical::IonLim(){
    return 9999999999;
}



TwoButyneDiNitrile::TwoButyneDiNitrile(){
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Nitrogen();
    abundancy[0] = 4;
    abundancy[1] = 2;
}

int TwoButyneDiNitrile::numberOfCostituents(){
    return 2;
}

double TwoButyneDiNitrile::getMass(){
    return costituents[0]->getMass() * abundancy[0] + costituents[1]->getMass() * abundancy[1];
}

std::string TwoButyneDiNitrile::getFormula(){
    return "C4N2";
}

double TwoButyneDiNitrile::formationEnergy(){
    return -5.4097e-18;
}

double TwoButyneDiNitrile::IonLim(){
    return 9999999999;
}



HydrogenCyanide::HydrogenCyanide(){
    costituents.resize(3);
    abundancy.resize(3);
    costituents[0] = new Hydrogen();
    costituents[1] = new Carbon();
    costituents[2] = new Nitrogen();
    abundancy[0] = 1;
    abundancy[1] = 1;
    abundancy[2] = 1;
}

int HydrogenCyanide::numberOfCostituents(){
    return 3;
}

double HydrogenCyanide::getMass(){
    return costituents[0]->getMass() * abundancy[0] + costituents[1]->getMass() * abundancy[1] + costituents[2]->getMass() * abundancy[2];
}

std::string HydrogenCyanide::getFormula(){
    return "CHN";
}

double HydrogenCyanide::formationEnergy(){
    return -2.0965e-18;
}

double HydrogenCyanide::IonLim(){
    return 9999999999;
}



CyanoRadical::CyanoRadical(){
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Carbon();
    costituents[1] = new Nitrogen();
    abundancy[0] = 1;
    abundancy[1] = 1;
}

int CyanoRadical::numberOfCostituents(){
    return 2;
}

double CyanoRadical::getMass(){
    return costituents[0]->getMass() * abundancy[0] + costituents[1]->getMass() * abundancy[1];
}

std::string CyanoRadical::getFormula(){
    return "CN";
}

double CyanoRadical::formationEnergy(){
    return -1.2457e-18;
}

double CyanoRadical::IonLim(){
    return 9999999999;
}



IsoCyanatoRadical::IsoCyanatoRadical(){
    costituents.resize(3);
    abundancy.resize(3);
    costituents[0] = new Carbon();
    costituents[1] = new Nitrogen();
    costituents[2] = new Oxygen();
    abundancy[0] = 1;
    abundancy[1] = 1;
    abundancy[2] = 1;
}

int IsoCyanatoRadical::numberOfCostituents(){
    return 3;
}

double IsoCyanatoRadical::getMass(){
    return costituents[0]->getMass() * abundancy[0] + costituents[1]->getMass() * abundancy[1] + costituents[2]->getMass() * abundancy[2];
}

std::string IsoCyanatoRadical::getFormula(){
    return "CNO";
}

double IsoCyanatoRadical::formationEnergy(){
    return -2.1086e-18;
}

double IsoCyanatoRadical::IonLim(){
    return 9999999999;
}



Imidogen::Imidogen(){
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Nitrogen();
    costituents[1] = new Hydrogen();
    abundancy[0] = 1;
    abundancy[1] = 1;
}

int Imidogen::numberOfCostituents(){
    return 2;
}

double Imidogen::getMass(){
    return costituents[0]->getMass() * abundancy[0] + costituents[1]->getMass() * abundancy[1];
}

std::string Imidogen::getFormula(){
    return "HN";
}

double Imidogen::formationEnergy(){
    return -5.1534e-19;
}

double Imidogen::IonLim(){
    return 9999999999;
}



NitrusOxide::NitrusOxide(){
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Nitrogen();
    costituents[1] = new Oxygen();
    abundancy[0] = 2;
    abundancy[1] = 1;
}

int NitrusOxide::numberOfCostituents(){
    return 2;
}

double NitrusOxide::getMass(){
    return costituents[0]->getMass() * abundancy[0] + costituents[1]->getMass() * abundancy[1];
}

std::string NitrusOxide::getFormula(){
    return "N2O";
}

double NitrusOxide::formationEnergy(){
    return -1.8315e-18;
}

double NitrusOxide::IonLim(){
    return 9999999999;
}



Ammonia::Ammonia(){
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Nitrogen();
    costituents[1] = new Hydrogen();
    abundancy[0] = 1;
    abundancy[1] = 3;
}

int Ammonia::numberOfCostituents(){
    return 2;
}

double Ammonia::getMass(){
    return costituents[0]->getMass() * abundancy[0] + costituents[1]->getMass() * abundancy[1];
}

std::string Ammonia::getFormula(){
    return "NH3";
}

double Ammonia::formationEnergy(){
    return -1.9226e-18;
}

double Ammonia::IonLim(){
    return 9999999999;
}



NitricOxide::NitricOxide(){
    costituents.resize(2);
    abundancy.resize(2);
    costituents[0] = new Nitrogen();
    costituents[1] = new Oxygen();
    abundancy[0] = 1;
    abundancy[1] = 1;
}

int NitricOxide::numberOfCostituents(){
    return 2;
}

double NitricOxide::getMass(){
    return costituents[0]->getMass() * abundancy[0] + costituents[1]->getMass() * abundancy[1];
}

std::string NitricOxide::getFormula(){
    return "NO";
}

double NitricOxide::formationEnergy(){
    return -1.0425e-18;
}

double NitricOxide::IonLim(){
    return 9999999999;
}

NitricOxideI::NitricOxideI() : NitricOxide() {}

int NitricOxideI::getCharge(){
    return 1;
}

std::string NitricOxideI::getFormula(){
    return "NO+";
}

double NitricOxideI::formationEnergy(){
    return NitricOxide::formationEnergy() + (9.26*eVtoJ);
}

Element* NitricOxideI::Constituent(){
    return costituents[0];
}