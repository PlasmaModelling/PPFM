 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include"Species.h"

double Aluminum::getMass() { return 26.9815385 * amuKg ; }

std::string Aluminum::getFormula() { return "Al" ; }

double Aluminum::formationEnergy() { return 0.0 * eVtoJ ; }

double Aluminum::IonLim() { return 5.9858 * eVtoJ ; }

double AluminumI::getMass() { return (26.98099 - 1*eamu) * amuKg ; }

int AluminumI::getCharge() { return 1 ; }

std::string AluminumI::getFormula() { return "Al+" ; }

double AluminumI::formationEnergy() { return 5.9858 * eVtoJ ; }

double AluminumI::IonLim() { return 18.828 * eVtoJ ; }

Element* AluminumI::Constituent() { return new Aluminum; }

double AluminumII::getMass() { return (26.98044 - 2*eamu) * amuKg ; }

int AluminumII::getCharge() { return 2 ; }

std::string AluminumII::getFormula() { return "Al+2" ; }

double AluminumII::formationEnergy() { return 24.8138 * eVtoJ ; }

double AluminumII::IonLim() { return 28.447 * eVtoJ ; }

Element* AluminumII::Constituent() { return new Aluminum; }

double AluminumIII::getMass() { return (26.97989 - 3*eamu) * amuKg ; }

int AluminumIII::getCharge() { return 3 ; }

std::string AluminumIII::getFormula() { return "Al+3" ; }

double AluminumIII::formationEnergy() { return 53.2618 * eVtoJ ; }

double AluminumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* AluminumIII::Constituent() { return new Aluminum; }

double Titanium::getMass() { return 47.867 * amuKg ; }

std::string Titanium::getFormula() { return "Ti" ; }

double Titanium::formationEnergy() { return 0.0 * eVtoJ ; }

double Titanium::IonLim() { return 6.8281 * eVtoJ ; }

double TitaniumI::getMass() { return (47.86645 - 1*eamu) * amuKg ; }

int TitaniumI::getCharge() { return 1 ; }

std::string TitaniumI::getFormula() { return "Ti+" ; }

double TitaniumI::formationEnergy() { return 6.8281 * eVtoJ ; }

double TitaniumI::IonLim() { return 13.5755 * eVtoJ ; }

Element* TitaniumI::Constituent() { return new Titanium; }

double TitaniumII::getMass() { return (47.8659 - 2*eamu) * amuKg ; }

int TitaniumII::getCharge() { return 2 ; }

std::string TitaniumII::getFormula() { return "Ti+2" ; }

double TitaniumII::formationEnergy() { return 20.4036 * eVtoJ ; }

double TitaniumII::IonLim() { return 27.4917 * eVtoJ ; }

Element* TitaniumII::Constituent() { return new Titanium; }

double TitaniumIII::getMass() { return (47.86535 - 3*eamu) * amuKg ; }

int TitaniumIII::getCharge() { return 3 ; }

std::string TitaniumIII::getFormula() { return "Ti+3" ; }

double TitaniumIII::formationEnergy() { return 47.8953 * eVtoJ ; }

double TitaniumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* TitaniumIII::Constituent() { return new Titanium; }

double Vanadium::getMass() { return 50.9415 * amuKg ; }

std::string Vanadium::getFormula() { return "V" ; }

double Vanadium::formationEnergy() { return 0.0 * eVtoJ ; }

double Vanadium::IonLim() { return 6.7462 * eVtoJ ; }

double VanadiumI::getMass() { return (50.94095 - 1*eamu) * amuKg ; }

int VanadiumI::getCharge() { return 1 ; }

std::string VanadiumI::getFormula() { return "V+" ; }

double VanadiumI::formationEnergy() { return 6.7462 * eVtoJ ; }

double VanadiumI::IonLim() { return 14.66 * eVtoJ ; }

Element* VanadiumI::Constituent() { return new Vanadium; }

double VanadiumII::getMass() { return (50.9404 - 2*eamu) * amuKg ; }

int VanadiumII::getCharge() { return 2 ; }

std::string VanadiumII::getFormula() { return "V+2" ; }

double VanadiumII::formationEnergy() { return 21.4062 * eVtoJ ; }

double VanadiumII::IonLim() { return 29.311 * eVtoJ ; }

Element* VanadiumII::Constituent() { return new Vanadium; }

double VanadiumIII::getMass() { return (50.93985 - 3*eamu) * amuKg ; }

int VanadiumIII::getCharge() { return 3 ; }

std::string VanadiumIII::getFormula() { return "V+3" ; }

double VanadiumIII::formationEnergy() { return 50.7172 * eVtoJ ; }

double VanadiumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* VanadiumIII::Constituent() { return new Vanadium; }

double Chromium::getMass() { return 51.9961 * amuKg ; }

std::string Chromium::getFormula() { return "Cr" ; }

double Chromium::formationEnergy() { return 0.0 * eVtoJ ; }

double Chromium::IonLim() { return 6.7665 * eVtoJ ; }

double ChromiumI::getMass() { return (51.99555 - 1*eamu) * amuKg ; }

int ChromiumI::getCharge() { return 1 ; }

std::string ChromiumI::getFormula() { return "Cr+" ; }

double ChromiumI::formationEnergy() { return 6.7665 * eVtoJ ; }

double ChromiumI::IonLim() { return 16.4857 * eVtoJ ; }

Element* ChromiumI::Constituent() { return new Chromium; }

double ChromiumII::getMass() { return (51.995 - 2*eamu) * amuKg ; }

int ChromiumII::getCharge() { return 2 ; }

std::string ChromiumII::getFormula() { return "Cr+2" ; }

double ChromiumII::formationEnergy() { return 23.2522 * eVtoJ ; }

double ChromiumII::IonLim() { return 30.96 * eVtoJ ; }

Element* ChromiumII::Constituent() { return new Chromium; }

double ChromiumIII::getMass() { return (51.99445 - 3*eamu) * amuKg ; }

int ChromiumIII::getCharge() { return 3 ; }

std::string ChromiumIII::getFormula() { return "Cr+3" ; }

double ChromiumIII::formationEnergy() { return 54.2122 * eVtoJ ; }

double ChromiumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* ChromiumIII::Constituent() { return new Chromium; }

double Manganese::getMass() { return 54.938044 * amuKg ; }

std::string Manganese::getFormula() { return "Mn" ; }

double Manganese::formationEnergy() { return 0.0 * eVtoJ ; }

double Manganese::IonLim() { return 7.43402 * eVtoJ ; }

double ManganeseI::getMass() { return (54.93749 - 1*eamu) * amuKg ; }

int ManganeseI::getCharge() { return 1 ; }

std::string ManganeseI::getFormula() { return "Mn+" ; }

double ManganeseI::formationEnergy() { return 7.43402 * eVtoJ ; }

double ManganeseI::IonLim() { return 15.64 * eVtoJ ; }

Element* ManganeseI::Constituent() { return new Manganese; }

double ManganeseII::getMass() { return (54.93694 - 2*eamu) * amuKg ; }

int ManganeseII::getCharge() { return 2 ; }

std::string ManganeseII::getFormula() { return "Mn+2" ; }

double ManganeseII::formationEnergy() { return 23.07402 * eVtoJ ; }

double ManganeseII::IonLim() { return 33.668 * eVtoJ ; }

Element* ManganeseII::Constituent() { return new Manganese; }

double ManganeseIII::getMass() { return (54.93639 - 3*eamu) * amuKg ; }

int ManganeseIII::getCharge() { return 3 ; }

std::string ManganeseIII::getFormula() { return "Mn+3" ; }

double ManganeseIII::formationEnergy() { return 56.74202 * eVtoJ ; }

double ManganeseIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* ManganeseIII::Constituent() { return new Manganese; }

double Iron::getMass() { return 55.845 * amuKg ; }

std::string Iron::getFormula() { return "Fe" ; }

double Iron::formationEnergy() { return 0.0 * eVtoJ ; }

double Iron::IonLim() { return 7.9024 * eVtoJ ; }

double IronI::getMass() { return (55.84445 - 1*eamu) * amuKg ; }

int IronI::getCharge() { return 1 ; }

std::string IronI::getFormula() { return "Fe+" ; }

double IronI::formationEnergy() { return 7.9024 * eVtoJ ; }

double IronI::IonLim() { return 16.1878 * eVtoJ ; }

Element* IronI::Constituent() { return new Iron; }

double IronII::getMass() { return (55.8439 - 2*eamu) * amuKg ; }

int IronII::getCharge() { return 2 ; }

std::string IronII::getFormula() { return "Fe+2" ; }

double IronII::formationEnergy() { return 24.0902 * eVtoJ ; }

double IronII::IonLim() { return 30.652 * eVtoJ ; }

Element* IronII::Constituent() { return new Iron; }

double IronIII::getMass() { return (55.84335 - 3*eamu) * amuKg ; }

int IronIII::getCharge() { return 3 ; }

std::string IronIII::getFormula() { return "Fe+3" ; }

double IronIII::formationEnergy() { return 54.7422 * eVtoJ ; }

double IronIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* IronIII::Constituent() { return new Iron; }

double Cobalt::getMass() { return 58.933194 * amuKg ; }

std::string Cobalt::getFormula() { return "Co" ; }

double Cobalt::formationEnergy() { return 0.0 * eVtoJ ; }

double Cobalt::IonLim() { return 7.88101 * eVtoJ ; }

double CobaltI::getMass() { return (58.93264 - 1*eamu) * amuKg ; }

int CobaltI::getCharge() { return 1 ; }

std::string CobaltI::getFormula() { return "Co+" ; }

double CobaltI::formationEnergy() { return 7.88101 * eVtoJ ; }

double CobaltI::IonLim() { return 17.0844 * eVtoJ ; }

Element* CobaltI::Constituent() { return new Cobalt; }

double CobaltII::getMass() { return (58.93208 - 2*eamu) * amuKg ; }

int CobaltII::getCharge() { return 2 ; }

std::string CobaltII::getFormula() { return "Co+2" ; }

double CobaltII::formationEnergy() { return 24.96541 * eVtoJ ; }

double CobaltII::IonLim() { return 33.5 * eVtoJ ; }

Element* CobaltII::Constituent() { return new Cobalt; }

double CobaltIII::getMass() { return (58.93152 - 3*eamu) * amuKg ; }

int CobaltIII::getCharge() { return 3 ; }

std::string CobaltIII::getFormula() { return "Co+3" ; }

double CobaltIII::formationEnergy() { return 58.46541 * eVtoJ ; }

double CobaltIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* CobaltIII::Constituent() { return new Cobalt; }

double Nickel::getMass() { return 58.6934 * amuKg ; }

std::string Nickel::getFormula() { return "Ni" ; }

double Nickel::formationEnergy() { return 0.0 * eVtoJ ; }

double Nickel::IonLim() { return 7.63978 * eVtoJ ; }

double NickelI::getMass() { return (58.69285 - 1*eamu) * amuKg ; }

int NickelI::getCharge() { return 1 ; }

std::string NickelI::getFormula() { return "Ni+" ; }

double NickelI::formationEnergy() { return 7.63978 * eVtoJ ; }

double NickelI::IonLim() { return 18.16806 * eVtoJ ; }

Element* NickelI::Constituent() { return new Nickel; }

double NickelII::getMass() { return (58.6923 - 2*eamu) * amuKg ; }

int NickelII::getCharge() { return 2 ; }

std::string NickelII::getFormula() { return "Ni+2" ; }

double NickelII::formationEnergy() { return 25.80784 * eVtoJ ; }

double NickelII::IonLim() { return 35.19745 * eVtoJ ; }

Element* NickelII::Constituent() { return new Nickel; }

double NickelIII::getMass() { return (58.69175 - 3*eamu) * amuKg ; }

int NickelIII::getCharge() { return 3 ; }

std::string NickelIII::getFormula() { return "Ni+3" ; }

double NickelIII::formationEnergy() { return 61.00529 * eVtoJ ; }

double NickelIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* NickelIII::Constituent() { return new Nickel; }

double Copper::getMass() { return 63.546 * amuKg ; }

std::string Copper::getFormula() { return "Cu" ; }

double Copper::formationEnergy() { return 0.0 * eVtoJ ; }

double Copper::IonLim() { return 7.72638 * eVtoJ ; }

double CopperI::getMass() { return (63.54544 - 1*eamu) * amuKg ; }

int CopperI::getCharge() { return 1 ; }

std::string CopperI::getFormula() { return "Cu+" ; }

double CopperI::formationEnergy() { return 7.72638 * eVtoJ ; }

double CopperI::IonLim() { return 20.29 * eVtoJ ; }

Element* CopperI::Constituent() { return new Copper; }

double CopperII::getMass() { return (63.54488 - 2*eamu) * amuKg ; }

int CopperII::getCharge() { return 2 ; }

std::string CopperII::getFormula() { return "Cu+2" ; }

double CopperII::formationEnergy() { return 28.01638 * eVtoJ ; }

double CopperII::IonLim() { return 36.84 * eVtoJ ; }

Element* CopperII::Constituent() { return new Copper; }

double CopperIII::getMass() { return (63.54432 - 3*eamu) * amuKg ; }

int CopperIII::getCharge() { return 3 ; }

std::string CopperIII::getFormula() { return "Cu+3" ; }

double CopperIII::formationEnergy() { return 64.85638 * eVtoJ ; }

double CopperIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* CopperIII::Constituent() { return new Copper; }

double Zinc::getMass() { return 65.38 * amuKg ; }

std::string Zinc::getFormula() { return "Zn" ; }

double Zinc::formationEnergy() { return 0.0 * eVtoJ ; }

double Zinc::IonLim() { return 9.39418 * eVtoJ ; }

double ZincI::getMass() { return (65.37963 - 1*eamu) * amuKg ; }

int ZincI::getCharge() { return 1 ; }

std::string ZincI::getFormula() { return "Zn+" ; }

double ZincI::formationEnergy() { return 9.39418 * eVtoJ ; }

double ZincI::IonLim() { return 17.9648 * eVtoJ ; }

Element* ZincI::Constituent() { return new Zinc; }

double ZincII::getMass() { return (65.37926 - 2*eamu) * amuKg ; }

int ZincII::getCharge() { return 2 ; }

std::string ZincII::getFormula() { return "Zn+2" ; }

double ZincII::formationEnergy() { return 27.35898 * eVtoJ ; }

double ZincII::IonLim() { return 39.722 * eVtoJ ; }

Element* ZincII::Constituent() { return new Zinc; }

double ZincIII::getMass() { return (65.37889 - 3*eamu) * amuKg ; }

int ZincIII::getCharge() { return 3 ; }

std::string ZincIII::getFormula() { return "Zn+3" ; }

double ZincIII::formationEnergy() { return 67.08098 * eVtoJ ; }

double ZincIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* ZincIII::Constituent() { return new Zinc; }

double Gallium::getMass() { return 69.723 * amuKg ; }

std::string Gallium::getFormula() { return "Ga" ; }

double Gallium::formationEnergy() { return 0.0 * eVtoJ ; }

double Gallium::IonLim() { return 5.999305 * eVtoJ ; }

double GalliumI::getMass() { return (69.72245 - 1*eamu) * amuKg ; }

int GalliumI::getCharge() { return 1 ; }

std::string GalliumI::getFormula() { return "Ga+" ; }

double GalliumI::formationEnergy() { return 5.999305 * eVtoJ ; }

double GalliumI::IonLim() { return 20.5152 * eVtoJ ; }

Element* GalliumI::Constituent() { return new Gallium; }

double GalliumII::getMass() { return (69.7219 - 2*eamu) * amuKg ; }

int GalliumII::getCharge() { return 2 ; }

std::string GalliumII::getFormula() { return "Ga+2" ; }

double GalliumII::formationEnergy() { return 26.514505 * eVtoJ ; }

double GalliumII::IonLim() { return 30.7 * eVtoJ ; }

Element* GalliumII::Constituent() { return new Gallium; }

double GalliumIII::getMass() { return (69.72135 - 3*eamu) * amuKg ; }

int GalliumIII::getCharge() { return 3 ; }

std::string GalliumIII::getFormula() { return "Ga+3" ; }

double GalliumIII::formationEnergy() { return 57.214505 * eVtoJ ; }

double GalliumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* GalliumIII::Constituent() { return new Gallium; }