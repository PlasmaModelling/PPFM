#include"Species.h"

double Beryllium::getMass() { return 9.012183 * amuKg ; }

std::string Beryllium::getFormula() { return "Be" ; }

double Beryllium::formationEnergy() { return 0.0 * eVtoJ ; }

double Beryllium::IonLim() { return 9.3227 * eVtoJ ; }

double BerylliumI::getMass() { return (9.01163 - 1*eamu) * amuKg ; }

int BerylliumI::getCharge() { return 1 ; }

std::string BerylliumI::getFormula() { return "Be+" ; }

double BerylliumI::formationEnergy() { return 9.3227 * eVtoJ ; }

double BerylliumI::IonLim() { return 18.211 * eVtoJ ; }

Element* BerylliumI::Constituent() { return new Beryllium; }

double BerylliumII::getMass() { return (9.01108 - 2*eamu) * amuKg ; }

int BerylliumII::getCharge() { return 2 ; }

std::string BerylliumII::getFormula() { return "Be+2" ; }

double BerylliumII::formationEnergy() { return 27.5337 * eVtoJ ; }

double BerylliumII::IonLim() { return 153.896 * eVtoJ ; }

Element* BerylliumII::Constituent() { return new Beryllium; }

double BerylliumIII::getMass() { return (9.01053 - 3*eamu) * amuKg ; }

int BerylliumIII::getCharge() { return 3 ; }

std::string BerylliumIII::getFormula() { return "Be+3" ; }

double BerylliumIII::formationEnergy() { return 181.4297 * eVtoJ ; }

double BerylliumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* BerylliumIII::Constituent() { return new Beryllium; }

double Magnesium::getMass() { return 24.305 * amuKg ; }

std::string Magnesium::getFormula() { return "Mg" ; }

double Magnesium::formationEnergy() { return 0.0 * eVtoJ ; }

double Magnesium::IonLim() { return 7.6462 * eVtoJ ; }

double MagnesiumI::getMass() { return (24.30445 - 1*eamu) * amuKg ; }

int MagnesiumI::getCharge() { return 1 ; }

std::string MagnesiumI::getFormula() { return "Mg+" ; }

double MagnesiumI::formationEnergy() { return 7.6462 * eVtoJ ; }

double MagnesiumI::IonLim() { return 15.0353 * eVtoJ ; }

Element* MagnesiumI::Constituent() { return new Magnesium; }

double MagnesiumII::getMass() { return (24.3039 - 2*eamu) * amuKg ; }

int MagnesiumII::getCharge() { return 2 ; }

std::string MagnesiumII::getFormula() { return "Mg+2" ; }

double MagnesiumII::formationEnergy() { return 22.6815 * eVtoJ ; }

double MagnesiumII::IonLim() { return 80.143 * eVtoJ ; }

Element* MagnesiumII::Constituent() { return new Magnesium; }

double MagnesiumIII::getMass() { return (24.30335 - 3*eamu) * amuKg ; }

int MagnesiumIII::getCharge() { return 3 ; }

std::string MagnesiumIII::getFormula() { return "Mg+3" ; }

double MagnesiumIII::formationEnergy() { return 102.8245 * eVtoJ ; }

double MagnesiumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* MagnesiumIII::Constituent() { return new Magnesium; }

double Calcium::getMass() { return 40.078 * amuKg ; }

std::string Calcium::getFormula() { return "Ca" ; }

double Calcium::formationEnergy() { return 0.0 * eVtoJ ; }

double Calcium::IonLim() { return 6.11316 * eVtoJ ; }

double CalciumI::getMass() { return (40.07754 - 1*eamu) * amuKg ; }

int CalciumI::getCharge() { return 1 ; }

std::string CalciumI::getFormula() { return "Ca+" ; }

double CalciumI::formationEnergy() { return 6.11316 * eVtoJ ; }

double CalciumI::IonLim() { return 11.87172 * eVtoJ ; }

Element* CalciumI::Constituent() { return new Calcium; }

double CalciumII::getMass() { return (40.07704 - 2*eamu) * amuKg ; }

int CalciumII::getCharge() { return 2 ; }

std::string CalciumII::getFormula() { return "Ca+2" ; }

double CalciumII::formationEnergy() { return 17.98488 * eVtoJ ; }

double CalciumII::IonLim() { return 50.9131 * eVtoJ ; }

Element* CalciumII::Constituent() { return new Calcium; }

double CalciumIII::getMass() { return (40.07654 - 3*eamu) * amuKg ; }

int CalciumIII::getCharge() { return 3 ; }

std::string CalciumIII::getFormula() { return "Ca+3" ; }

double CalciumIII::formationEnergy() { return 68.89798 * eVtoJ ; }

double CalciumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* CalciumIII::Constituent() { return new Calcium; }

double Strontium::getMass() { return 87.62 * amuKg ; }

std::string Strontium::getFormula() { return "Sr" ; }

double Strontium::formationEnergy() { return 0.0 * eVtoJ ; }

double Strontium::IonLim() { return 5.69488 * eVtoJ ; }

double StrontiumI::getMass() { return (87.61943 - 1*eamu) * amuKg ; }

int StrontiumI::getCharge() { return 1 ; }

std::string StrontiumI::getFormula() { return "Sr+" ; }

double StrontiumI::formationEnergy() { return 5.69488 * eVtoJ ; }

double StrontiumI::IonLim() { return 11.03036 * eVtoJ ; }

Element* StrontiumI::Constituent() { return new Strontium; }

double StrontiumII::getMass() { return (87.61887 - 2*eamu) * amuKg ; }

int StrontiumII::getCharge() { return 2 ; }

std::string StrontiumII::getFormula() { return "Sr+2" ; }

double StrontiumII::formationEnergy() { return 16.72524 * eVtoJ ; }

double StrontiumII::IonLim() { return 31.93 * eVtoJ ; }

Element* StrontiumII::Constituent() { return new Strontium; }

double StrontiumIII::getMass() { return (87.6183 - 3*eamu) * amuKg ; }

int StrontiumIII::getCharge() { return 3 ; }

std::string StrontiumIII::getFormula() { return "Sr+3" ; }

double StrontiumIII::formationEnergy() { return 48.65524 * eVtoJ ; }

double StrontiumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* StrontiumIII::Constituent() { return new Strontium; }

double Barium::getMass() { return 137.327 * amuKg ; }

std::string Barium::getFormula() { return "Ba" ; }

double Barium::formationEnergy() { return 0.0 * eVtoJ ; }

double Barium::IonLim() { return 5.2117 * eVtoJ ; }

double BariumI::getMass() { return (137.32615 - 1*eamu) * amuKg ; }

int BariumI::getCharge() { return 1 ; }

std::string BariumI::getFormula() { return "Ba+" ; }

double BariumI::formationEnergy() { return 5.2117 * eVtoJ ; }

double BariumI::IonLim() { return 10.0039 * eVtoJ ; }

Element* BariumI::Constituent() { return new Barium; }

double BariumII::getMass() { return (137.3253 - 2*eamu) * amuKg ; }

int BariumII::getCharge() { return 2 ; }

std::string BariumII::getFormula() { return "Ba+2" ; }

double BariumII::formationEnergy() { return 15.2156 * eVtoJ ; }

double BariumII::IonLim() { return 37.055 * eVtoJ ; }

Element* BariumII::Constituent() { return new Barium; }

double BariumIII::getMass() { return (137.32445 - 3*eamu) * amuKg ; }

int BariumIII::getCharge() { return 3 ; }

std::string BariumIII::getFormula() { return "Ba+3" ; }

double BariumIII::formationEnergy() { return 52.2706 * eVtoJ ; }

double BariumIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* BariumIII::Constituent() { return new Barium; }