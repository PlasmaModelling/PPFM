 // PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
 // (University of Bologna, Italy)                   // 
 // Licensed under CC BY 4.0.                        // 
 // To view a copy of this license, visit:           // 
 // https://creativecommons.org/licenses/by/4.0/     // 

#include"Species.h"

double Silver::getMass() { return 107.8682 * amuKg ; }

std::string Silver::getFormula() { return "Ag" ; }

double Silver::formationEnergy() { return 0.0 * eVtoJ ; }

double Silver::IonLim() { return 7.57624 * eVtoJ ; }

double SilverI::getMass() { return (107.86765 - 1*eamu) * amuKg ; }

int SilverI::getCharge() { return 1 ; }

std::string SilverI::getFormula() { return "Ag+" ; }

double SilverI::formationEnergy() { return 7.57624 * eVtoJ ; }

double SilverI::IonLim() { return 21.4798 * eVtoJ ; }

Element* SilverI::Constituent() { return new Silver; }

double SilverII::getMass() { return (107.8671 - 2*eamu) * amuKg ; }

int SilverII::getCharge() { return 2 ; }

std::string SilverII::getFormula() { return "Ag+2" ; }

double SilverII::formationEnergy() { return 29.05604 * eVtoJ ; }

double SilverII::IonLim() { return 34.83 * eVtoJ ; }

Element* SilverII::Constituent() { return new Silver; }

double SilverIII::getMass() { return (107.86655 - 3*eamu) * amuKg ; }

int SilverIII::getCharge() { return 3 ; }

std::string SilverIII::getFormula() { return "Ag+3" ; }

double SilverIII::formationEnergy() { return 63.88604 * eVtoJ ; }

double SilverIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* SilverIII::Constituent() { return new Silver; }

double Gold::getMass() { return 196.966569 * amuKg ; }

std::string Gold::getFormula() { return "Au" ; }

double Gold::formationEnergy() { return 0.0 * eVtoJ ; }

double Gold::IonLim() { return 9.22551 * eVtoJ ; }

double GoldI::getMass() { return (196.96602 - 1*eamu) * amuKg ; }

int GoldI::getCharge() { return 1 ; }

std::string GoldI::getFormula() { return "Au+" ; }

double GoldI::formationEnergy() { return 9.22551 * eVtoJ ; }

double GoldI::IonLim() { return 20.5147 * eVtoJ ; }

Element* GoldI::Constituent() { return new Gold; }

double GoldII::getMass() { return (196.96547 - 2*eamu) * amuKg ; }

int GoldII::getCharge() { return 2 ; }

std::string GoldII::getFormula() { return "Au+2" ; }

double GoldII::formationEnergy() { return 29.74021 * eVtoJ ; }

double GoldII::IonLim() { return 33.88 * eVtoJ ; }

Element* GoldII::Constituent() { return new Gold; }

double GoldIII::getMass() { return (196.96492 - 3*eamu) * amuKg ; }

int GoldIII::getCharge() { return 3 ; }

std::string GoldIII::getFormula() { return "Au+3" ; }

double GoldIII::formationEnergy() { return 63.62121 * eVtoJ ; }

double GoldIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* GoldIII::Constituent() { return new Gold; }

double Mercury::getMass() { return 200.592 * amuKg ; }

std::string Mercury::getFormula() { return "Hg" ; }

double Mercury::formationEnergy() { return 0.0 * eVtoJ ; }

double Mercury::IonLim() { return 10.4384 * eVtoJ ; }

double MercuryI::getMass() { return (200.59145 - 1*eamu) * amuKg ; }

int MercuryI::getCharge() { return 1 ; }

std::string MercuryI::getFormula() { return "Hg+" ; }

double MercuryI::formationEnergy() { return 10.4384 * eVtoJ ; }

double MercuryI::IonLim() { return 18.7556 * eVtoJ ; }

Element* MercuryI::Constituent() { return new Mercury; }

double MercuryII::getMass() { return (200.5909 - 2*eamu) * amuKg ; }

int MercuryII::getCharge() { return 2 ; }

std::string MercuryII::getFormula() { return "Hg+2" ; }

double MercuryII::formationEnergy() { return 29.194 * eVtoJ ; }

double MercuryII::IonLim() { return 34.0 * eVtoJ ; }

Element* MercuryII::Constituent() { return new Mercury; }

double MercuryIII::getMass() { return (200.59035 - 3*eamu) * amuKg ; }

int MercuryIII::getCharge() { return 3 ; }

std::string MercuryIII::getFormula() { return "Hg+3" ; }

double MercuryIII::formationEnergy() { return 63.194 * eVtoJ ; }

double MercuryIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* MercuryIII::Constituent() { return new Mercury; }

double Platinum::getMass() { return 195.084 * amuKg ; }

std::string Platinum::getFormula() { return "Pt" ; }

double Platinum::formationEnergy() { return 0.0 * eVtoJ ; }

double Platinum::IonLim() { return 9.00041 * eVtoJ ; }

double PlatinumI::getMass() { return (195.08345 - 1*eamu) * amuKg ; }

int PlatinumI::getCharge() { return 1 ; }

std::string PlatinumI::getFormula() { return "Pt+" ; }

double PlatinumI::formationEnergy() { return 9.00041 * eVtoJ ; }

double PlatinumI::IonLim() { return 19.176 * eVtoJ ; }

Element* PlatinumI::Constituent() { return new Platinum; }

double PlatinumII::getMass() { return (195.0829 - 2*eamu) * amuKg ; }

int PlatinumII::getCharge() { return 2 ; }

std::string PlatinumII::getFormula() { return "Pt+2" ; }

double PlatinumII::formationEnergy() { return 28.17641 * eVtoJ ; }

double PlatinumII::IonLim() { return 31.8 * eVtoJ ; }

Element* PlatinumII::Constituent() { return new Platinum; }

double PlatinumIII::getMass() { return (195.08235 - 3*eamu) * amuKg ; }

int PlatinumIII::getCharge() { return 3 ; }

std::string PlatinumIII::getFormula() { return "Pt+3" ; }

double PlatinumIII::formationEnergy() { return 59.97641 * eVtoJ ; }

double PlatinumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* PlatinumIII::Constituent() { return new Platinum; }

double Lead::getMass() { return 207.2 * amuKg ; }

std::string Lead::getFormula() { return "Pb" ; }

double Lead::formationEnergy() { return 0.0 * eVtoJ ; }

double Lead::IonLim() { return 7.41667 * eVtoJ ; }

double LeadI::getMass() { return (207.19945 - 1*eamu) * amuKg ; }

int LeadI::getCharge() { return 1 ; }

std::string LeadI::getFormula() { return "Pb+" ; }

double LeadI::formationEnergy() { return 7.41667 * eVtoJ ; }

double LeadI::IonLim() { return 15.032 * eVtoJ ; }

Element* LeadI::Constituent() { return new Lead; }

double LeadII::getMass() { return (207.1989 - 2*eamu) * amuKg ; }

int LeadII::getCharge() { return 2 ; }

std::string LeadII::getFormula() { return "Pb+2" ; }

double LeadII::formationEnergy() { return 22.44867 * eVtoJ ; }

double LeadII::IonLim() { return 31.4 * eVtoJ ; }

Element* LeadII::Constituent() { return new Lead; }

double LeadIII::getMass() { return (207.19835 - 3*eamu) * amuKg ; }

int LeadIII::getCharge() { return 3 ; }

std::string LeadIII::getFormula() { return "Pb+3" ; }

double LeadIII::formationEnergy() { return 53.84867 * eVtoJ ; }

double LeadIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* LeadIII::Constituent() { return new Lead; }

double Tantalum::getMass() { return 180.94788 * amuKg ; }

std::string Tantalum::getFormula() { return "Ta" ; }

double Tantalum::formationEnergy() { return 0.0 * eVtoJ ; }

double Tantalum::IonLim() { return 7.54957 * eVtoJ ; }

double TantalumI::getMass() { return (180.94733 - 1*eamu) * amuKg ; }

int TantalumI::getCharge() { return 1 ; }

std::string TantalumI::getFormula() { return "Ta+" ; }

double TantalumI::formationEnergy() { return 7.54957 * eVtoJ ; }

double TantalumI::IonLim() { return 14.001 * eVtoJ ; }

Element* TantalumI::Constituent() { return new Tantalum; }

double TantalumII::getMass() { return (180.94678 - 2*eamu) * amuKg ; }

int TantalumII::getCharge() { return 2 ; }

std::string TantalumII::getFormula() { return "Ta+2" ; }

double TantalumII::formationEnergy() { return 21.55057 * eVtoJ ; }

double TantalumII::IonLim() { return 31.0 * eVtoJ ; }

Element* TantalumII::Constituent() { return new Tantalum; }

double TantalumIII::getMass() { return (180.94623 - 3*eamu) * amuKg ; }

int TantalumIII::getCharge() { return 3 ; }

std::string TantalumIII::getFormula() { return "Ta+3" ; }

double TantalumIII::formationEnergy() { return 52.55057 * eVtoJ ; }

double TantalumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* TantalumIII::Constituent() { return new Tantalum; }

double Tungsten::getMass() { return 183.84 * amuKg ; }

std::string Tungsten::getFormula() { return "W" ; }

double Tungsten::formationEnergy() { return 0.0 * eVtoJ ; }

double Tungsten::IonLim() { return 7.864 * eVtoJ ; }

double TungstenI::getMass() { return (183.83945 - 1*eamu) * amuKg ; }

int TungstenI::getCharge() { return 1 ; }

std::string TungstenI::getFormula() { return "W+" ; }

double TungstenI::formationEnergy() { return 7.864 * eVtoJ ; }

double TungstenI::IonLim() { return 17.964 * eVtoJ ; }

Element* TungstenI::Constituent() { return new Tungsten; }

double TungstenII::getMass() { return (183.8389 - 2*eamu) * amuKg ; }

int TungstenII::getCharge() { return 2 ; }

std::string TungstenII::getFormula() { return "W+2" ; }

double TungstenII::formationEnergy() { return 25.828 * eVtoJ ; }

double TungstenII::IonLim() { return 36.0 * eVtoJ ; }

Element* TungstenII::Constituent() { return new Tungsten; }

double TungstenIII::getMass() { return (183.83835 - 3*eamu) * amuKg ; }

int TungstenIII::getCharge() { return 3 ; }

std::string TungstenIII::getFormula() { return "W+3" ; }

double TungstenIII::formationEnergy() { return 61.828 * eVtoJ ; }

double TungstenIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* TungstenIII::Constituent() { return new Tungsten; }

double Rhenium::getMass() { return 186.207 * amuKg ; }

std::string Rhenium::getFormula() { return "Re" ; }

double Rhenium::formationEnergy() { return 0.0 * eVtoJ ; }

double Rhenium::IonLim() { return 7.83352 * eVtoJ ; }

double RheniumI::getMass() { return (186.20645 - 1*eamu) * amuKg ; }

int RheniumI::getCharge() { return 1 ; }

std::string RheniumI::getFormula() { return "Re+" ; }

double RheniumI::formationEnergy() { return 7.83352 * eVtoJ ; }

double RheniumI::IonLim() { return 14.54 * eVtoJ ; }

Element* RheniumI::Constituent() { return new Rhenium; }

double RheniumII::getMass() { return (186.2059 - 2*eamu) * amuKg ; }

int RheniumII::getCharge() { return 2 ; }

std::string RheniumII::getFormula() { return "Re+2" ; }

double RheniumII::formationEnergy() { return 22.37352 * eVtoJ ; }

double RheniumII::IonLim() { return 33.5 * eVtoJ ; }

Element* RheniumII::Constituent() { return new Rhenium; }

double RheniumIII::getMass() { return (186.20535 - 3*eamu) * amuKg ; }

int RheniumIII::getCharge() { return 3 ; }

std::string RheniumIII::getFormula() { return "Re+3" ; }

double RheniumIII::formationEnergy() { return 55.87352 * eVtoJ ; }

double RheniumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* RheniumIII::Constituent() { return new Rhenium; }

double Molybdenum::getMass() { return 95.95 * amuKg ; }

std::string Molybdenum::getFormula() { return "Mo" ; }

double Molybdenum::formationEnergy() { return 0.0 * eVtoJ ; }

double Molybdenum::IonLim() { return 7.09238 * eVtoJ ; }

double MolybdenumI::getMass() { return (95.94945 - 1*eamu) * amuKg ; }

int MolybdenumI::getCharge() { return 1 ; }

std::string MolybdenumI::getFormula() { return "Mo+" ; }

double MolybdenumI::formationEnergy() { return 7.09238 * eVtoJ ; }

double MolybdenumI::IonLim() { return 16.151 * eVtoJ ; }

Element* MolybdenumI::Constituent() { return new Molybdenum; }

double MolybdenumII::getMass() { return (95.9489 - 2*eamu) * amuKg ; }

int MolybdenumII::getCharge() { return 2 ; }

std::string MolybdenumII::getFormula() { return "Mo+2" ; }

double MolybdenumII::formationEnergy() { return 23.24338 * eVtoJ ; }

double MolybdenumII::IonLim() { return 27.21 * eVtoJ ; }

Element* MolybdenumII::Constituent() { return new Molybdenum; }

double MolybdenumIII::getMass() { return (95.94835 - 3*eamu) * amuKg ; }

int MolybdenumIII::getCharge() { return 3 ; }

std::string MolybdenumIII::getFormula() { return "Mo+3" ; }

double MolybdenumIII::formationEnergy() { return 50.45338 * eVtoJ ; }

double MolybdenumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* MolybdenumIII::Constituent() { return new Molybdenum; }

double Niobium::getMass() { return 92.90637 * amuKg ; }

std::string Niobium::getFormula() { return "Nb" ; }

double Niobium::formationEnergy() { return 0.0 * eVtoJ ; }

double Niobium::IonLim() { return 6.75885 * eVtoJ ; }

double NiobiumI::getMass() { return (92.90582 - 1*eamu) * amuKg ; }

int NiobiumI::getCharge() { return 1 ; }

std::string NiobiumI::getFormula() { return "Nb+" ; }

double NiobiumI::formationEnergy() { return 6.75885 * eVtoJ ; }

double NiobiumI::IonLim() { return 14.315 * eVtoJ ; }

Element* NiobiumI::Constituent() { return new Niobium; }

double NiobiumII::getMass() { return (92.90527 - 2*eamu) * amuKg ; }

int NiobiumII::getCharge() { return 2 ; }

std::string NiobiumII::getFormula() { return "Nb+2" ; }

double NiobiumII::formationEnergy() { return 21.07385 * eVtoJ ; }

double NiobiumII::IonLim() { return 25.04 * eVtoJ ; }

Element* NiobiumII::Constituent() { return new Niobium; }

double NiobiumIII::getMass() { return (92.90472 - 3*eamu) * amuKg ; }

int NiobiumIII::getCharge() { return 3 ; }

std::string NiobiumIII::getFormula() { return "Nb+3" ; }

double NiobiumIII::formationEnergy() { return 46.11385 * eVtoJ ; }

double NiobiumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* NiobiumIII::Constituent() { return new Niobium; }

double Zirconium::getMass() { return 91.224 * amuKg ; }

std::string Zirconium::getFormula() { return "Zr" ; }

double Zirconium::formationEnergy() { return 0.0 * eVtoJ ; }

double Zirconium::IonLim() { return 6.6339 * eVtoJ ; }

double ZirconiumI::getMass() { return (91.22345 - 1*eamu) * amuKg ; }

int ZirconiumI::getCharge() { return 1 ; }

std::string ZirconiumI::getFormula() { return "Zr+" ; }

double ZirconiumI::formationEnergy() { return 6.6339 * eVtoJ ; }

double ZirconiumI::IonLim() { return 13.13 * eVtoJ ; }

Element* ZirconiumI::Constituent() { return new Zirconium; }

double ZirconiumII::getMass() { return (91.2229 - 2*eamu) * amuKg ; }

int ZirconiumII::getCharge() { return 2 ; }

std::string ZirconiumII::getFormula() { return "Zr+2" ; }

double ZirconiumII::formationEnergy() { return 19.7639 * eVtoJ ; }

double ZirconiumII::IonLim() { return 22.99 * eVtoJ ; }

Element* ZirconiumII::Constituent() { return new Zirconium; }

double ZirconiumIII::getMass() { return (91.22235 - 3*eamu) * amuKg ; }

int ZirconiumIII::getCharge() { return 3 ; }

std::string ZirconiumIII::getFormula() { return "Zr+3" ; }

double ZirconiumIII::formationEnergy() { return 42.7539 * eVtoJ ; }

double ZirconiumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* ZirconiumIII::Constituent() { return new Zirconium; }