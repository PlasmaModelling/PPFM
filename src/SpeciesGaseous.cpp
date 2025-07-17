#include"Species.h"

double Hydrogen::getMass() { return 1.00782503223 * amuKg ; }

std::string Hydrogen::getFormula() { return "H" ; }

double Hydrogen::formationEnergy() { return 0.0 * eVtoJ ; }

double Hydrogen::IonLim() { return 13.598434599702 * eVtoJ ; }

double HydrogenI::getMass() { return (1.007259 - 1*eamu) * amuKg ; }

int HydrogenI::getCharge() { return 1 ; }

std::string HydrogenI::getFormula() { return "H+" ; }

double HydrogenI::formationEnergy() { return 13.598434599702 * eVtoJ ; }

double HydrogenI::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* HydrogenI::Constituent() { return new Hydrogen; }

double HydrogenAnion::getMass() { return 1.6737e-27 /* 1.00849 * amuKg  */; }

int HydrogenAnion::getCharge() { return -1 ; }

std::string HydrogenAnion::getFormula() { return "H-" ; }

double HydrogenAnion::formationEnergy() { return -1.208400e-19 /* -0.75497 * eVtoJ */ ; }

double HydrogenAnion::IonLim() { return (-formationEnergy()) ; }

Element* HydrogenAnion::Constituent() { return new Hydrogen; }

double Helium::getMass() { return 4.002602 * amuKg ; }

std::string Helium::getFormula() { return "He" ; }

double Helium::formationEnergy() { return 0.0 * eVtoJ ; }

double Helium::IonLim() { return 24.587389011 * eVtoJ ; }

double HeliumI::getMass() { return (4.002096 - 1*eamu) * amuKg ; }

int HeliumI::getCharge() { return 1 ; }

std::string HeliumI::getFormula() { return "He+" ; }

double HeliumI::formationEnergy() { return 24.587389011 * eVtoJ ; }

double HeliumI::IonLim() { return 54.417760582 * eVtoJ ; }

Element* HeliumI::Constituent() { return new Helium; }

double HeliumII::getMass() { return (4.001589 - 2*eamu) * amuKg ; }

int HeliumII::getCharge() { return 2 ; }

std::string HeliumII::getFormula() { return "He+2" ; }

double HeliumII::formationEnergy() { return 79.005149593 * eVtoJ ; }

double HeliumII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* HeliumII::Constituent() { return new Helium; }

double Nitrogen::getMass() { return 14.0067 * amuKg ; }

std::string Nitrogen::getFormula() { return "N" ; }

double Nitrogen::formationEnergy() { return 0.0 * eVtoJ ; }

double Nitrogen::IonLim() { return 14.534139 * eVtoJ ; }

double NitrogenI::getMass() { return (14.00615 - 1*eamu) * amuKg ; }

int NitrogenI::getCharge() { return 1 ; }

std::string NitrogenI::getFormula() { return "N+" ; }

double NitrogenI::formationEnergy() { return 14.534139 * eVtoJ ; }

double NitrogenI::IonLim() { return 29.6013 * eVtoJ ; }

Element* NitrogenI::Constituent() { return new Nitrogen; }

double NitrogenII::getMass() { return (14.00561 - 2*eamu) * amuKg ; }

int NitrogenII::getCharge() { return 2 ; }

std::string NitrogenII::getFormula() { return "N+2" ; }

double NitrogenII::formationEnergy() { return 44.135439 * eVtoJ ; }

Element* NitrogenII::Constituent() { return new Nitrogen; }

double NitrogenII::IonLim() { return 47.44924 * eVtoJ ; }

double NitrogenIII::getMass() { return (14.00512 - 3*eamu) * amuKg ; }

int NitrogenIII::getCharge() { return 3 ; }

std::string NitrogenIII::getFormula() { return "N+3" ; }

double NitrogenIII::formationEnergy() { return 91.584679 * eVtoJ ; }

double NitrogenIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* NitrogenIII::Constituent() { return new Nitrogen; }

double NitrogenAnion::getMass() { return 14.0067 * amuKg ; }

int NitrogenAnion::getCharge() { return -1 ; }

std::string NitrogenAnion::getFormula() { return "N-" ; }

double NitrogenAnion::formationEnergy() { return 0.07 * eVtoJ ; }

double NitrogenAnion::IonLim() { return (-formationEnergy()) ; }

Element* NitrogenAnion::Constituent() { return new Nitrogen; }

double Oxygen::getMass() { return 15.999 * amuKg ; }

std::string Oxygen::getFormula() { return "O" ; }

double Oxygen::formationEnergy() { return 0.0 * eVtoJ ; }

double Oxygen::IonLim() { return 13.61806 * eVtoJ ; }

double OxygenI::getMass() { return (15.99845 - 1*eamu) * amuKg ; }

int OxygenI::getCharge() { return 1 ; }

std::string OxygenI::getFormula() { return "O+" ; }

double OxygenI::formationEnergy() { return 13.61806 * eVtoJ ; }

double OxygenI::IonLim() { return 35.1173 * eVtoJ ; }

Element* OxygenI::Constituent() { return new Oxygen; }

double OxygenII::getMass() { return (15.99791 - 2*eamu) * amuKg ; }

int OxygenII::getCharge() { return 2 ; }

std::string OxygenII::getFormula() { return "O+2" ; }

double OxygenII::formationEnergy() { return 48.73536 * eVtoJ ; }

double OxygenII::IonLim() { return 54.9355 * eVtoJ ; }

Element* OxygenII::Constituent() { return new Oxygen; }

double OxygenIII::getMass() { return (15.99743 - 3*eamu) * amuKg ; }

int OxygenIII::getCharge() { return 3 ; }

std::string OxygenIII::getFormula() { return "O+3" ; }

double OxygenIII::formationEnergy() { return 103.67086 * eVtoJ ; }

double OxygenIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* OxygenIII::Constituent() { return new Oxygen; }

double OxygenAnion::getMass() { return (15.9994 + eamu) * amuKg ; }

int OxygenAnion::getCharge() { return -1 ; }

std::string OxygenAnion::getFormula() { return "O-" ; }

double OxygenAnion::formationEnergy() { return -1.439157 * eVtoJ ; }

double OxygenAnion::IonLim() { return (-formationEnergy()) ; }

Element* OxygenAnion::Constituent() { return new Oxygen; }

double Fluorine::getMass() { return 18.99840316 * amuKg ; }

std::string Fluorine::getFormula() { return "F" ; }

double Fluorine::formationEnergy() { return 0.0 * eVtoJ ; }

double Fluorine::IonLim() { return 17.42282 * eVtoJ ; }

double FluorineI::getMass() { return (18.997915 - 1*eamu) * amuKg ; }

int FluorineI::getCharge() { return 1 ; }

std::string FluorineI::getFormula() { return "F+" ; }

double FluorineI::formationEnergy() { return 17.42282 * eVtoJ ; }

double FluorineI::IonLim() { return 34.97082 * eVtoJ ; }

double FluorineII::getMass() { return (18.99743 - 2*eamu) * amuKg ; }

Element* FluorineI::Constituent() { return new Fluorine; }

int FluorineII::getCharge() { return 2 ; }

std::string FluorineII::getFormula() { return "F+2" ; }

double FluorineII::formationEnergy() { return 52.39364 * eVtoJ ; }

double FluorineII::IonLim() { return 62.7084 * eVtoJ ; }

Element* FluorineII::Constituent() { return new Fluorine; }

double FluorineIII::getMass() { return (18.99694 - 3*eamu) * amuKg ; }

int FluorineIII::getCharge() { return 3 ; }

std::string FluorineIII::getFormula() { return "F+3" ; }

double FluorineIII::formationEnergy() { return 115.10204 * eVtoJ ; }

double FluorineIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* FluorineIII::Constituent() { return new Fluorine; }

double Neon::getMass() { return 20.1797 * amuKg ; }

std::string Neon::getFormula() { return "Ne" ; }

double Neon::formationEnergy() { return 0.0 * eVtoJ ; }

double Neon::IonLim() { return 21.5646 * eVtoJ ; }

double NeonI::getMass() { return (20.17916 - 1*eamu) * amuKg ; }

int NeonI::getCharge() { return 1 ; }

std::string NeonI::getFormula() { return "Ne+" ; }

double NeonI::formationEnergy() { return 21.5646 * eVtoJ ; }

double NeonI::IonLim() { return 40.96328 * eVtoJ ; }

Element* NeonI::Constituent() { return new Neon; }

double NeonII::getMass() { return (20.17862 - 2*eamu) * amuKg ; }

int NeonII::getCharge() { return 2 ; }

std::string NeonII::getFormula() { return "Ne+2" ; }

double NeonII::formationEnergy() { return 62.52788 * eVtoJ ; }

double NeonII::IonLim() { return 63.45 * eVtoJ ; }

Element* NeonII::Constituent() { return new Neon; }

double NeonIII::getMass() { return (20.17808 - 3*eamu) * amuKg ; }

int NeonIII::getCharge() { return 3 ; }

std::string NeonIII::getFormula() { return "Ne+3" ; }

double NeonIII::formationEnergy() { return 125.97788 * eVtoJ ; }

double NeonIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* NeonIII::Constituent() { return new Neon; }

double Chlorine::getMass() { return 35.4527 * amuKg ; }

std::string Chlorine::getFormula() { return "Cl" ; }

double Chlorine::formationEnergy() { return 0.0 * eVtoJ ; }

double Chlorine::IonLim() { return 12.96764 * eVtoJ ; }

double ChlorineI::getMass() { return (35.45215 - 1*eamu) * amuKg ; }

int ChlorineI::getCharge() { return 1 ; }

std::string ChlorineI::getFormula() { return "Cl+" ; }

double ChlorineI::formationEnergy() { return 12.96764 * eVtoJ ; }

double ChlorineI::IonLim() { return 23.814 * eVtoJ ; }

Element* ChlorineI::Constituent() { return new Chlorine; }

double ChlorineII::getMass() { return (35.4516 - 2*eamu) * amuKg ; }

int ChlorineII::getCharge() { return 2 ; }

std::string ChlorineII::getFormula() { return "Cl+2" ; }

double ChlorineII::formationEnergy() { return 36.78164 * eVtoJ ; }

double ChlorineII::IonLim() { return 39.61 * eVtoJ ; }

Element* ChlorineII::Constituent() { return new Chlorine; }

double ChlorineIII::getMass() { return (35.45105 - 3*eamu) * amuKg ; }

int ChlorineIII::getCharge() { return 3 ; }

std::string ChlorineIII::getFormula() { return "Cl+3" ; }

double ChlorineIII::formationEnergy() { return 76.39164 * eVtoJ ; }

double ChlorineIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* ChlorineIII::Constituent() { return new Chlorine; }

double Argon::getMass() { return 39.948 * amuKg ; }

std::string Argon::getFormula() { return "Ar" ; }

double Argon::formationEnergy() { return 0.0 * eVtoJ ; }

double Argon::IonLim() { return 15.75962 * eVtoJ ; }

double ArgonI::getMass() { return (39.94745 - 1*eamu) * amuKg ; }

int ArgonI::getCharge() { return 1 ; }

std::string ArgonI::getFormula() { return "Ar+" ; }

double ArgonI::formationEnergy() { return 15.75962 * eVtoJ ; }

double ArgonI::IonLim() { return 27.62967 * eVtoJ ; }

Element* ArgonI::Constituent() { return new Argon; }

double ArgonII::getMass() { return (39.9469 - 2*eamu) * amuKg ; }

int ArgonII::getCharge() { return 2 ; }

std::string ArgonII::getFormula() { return "Ar+2" ; }

double ArgonII::formationEnergy() { return 43.38929 * eVtoJ ; }

double ArgonII::IonLim() { return 40.74 * eVtoJ ; }

Element* ArgonII::Constituent() { return new Argon; }

double ArgonIII::getMass() { return (39.94635 - 3*eamu) * amuKg ; }

int ArgonIII::getCharge() { return 3 ; }

std::string ArgonIII::getFormula() { return "Ar+3" ; }

double ArgonIII::formationEnergy() { return 84.12929 * eVtoJ ; }

double ArgonIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* ArgonIII::Constituent() { return new Argon; }

double ArgonIV::getMass() { return (39.947 - 3.*eamu ) * amuKg ; }

int ArgonIV::getCharge() { return 4 ; }

std::string ArgonIV::getFormula() { return "Ar+4" ; }

double ArgonIV::formationEnergy() { return (27.62967 + 15.7596119 + 40.735 + 59.58) * eVtoJ ; }

double ArgonIV::IonLim() { return 74.84 * eVtoJ ; }

Element* ArgonIV::Constituent() { return new Argon; }

double Krypton::getMass() { return 83.798 * amuKg ; }

std::string Krypton::getFormula() { return "Kr" ; }

double Krypton::formationEnergy() { return 0.0 * eVtoJ ; }

double Krypton::IonLim() { return 13.99961 * eVtoJ ; }

double KryptonI::getMass() { return (83.79728 - 1*eamu) * amuKg ; }

int KryptonI::getCharge() { return 1 ; }

std::string KryptonI::getFormula() { return "Kr+" ; }

double KryptonI::formationEnergy() { return 13.99961 * eVtoJ ; }

double KryptonI::IonLim() { return 24.35984 * eVtoJ ; }

Element* KryptonI::Constituent() { return new Krypton; }

double KryptonII::getMass() { return (83.79656 - 2*eamu) * amuKg ; }

int KryptonII::getCharge() { return 2 ; }

std::string KryptonII::getFormula() { return "Kr+2" ; }

double KryptonII::formationEnergy() { return 38.35945 * eVtoJ ; }

double KryptonII::IonLim() { return 35.0 * eVtoJ ; }

Element* KryptonII::Constituent() { return new Krypton; }

double KryptonIII::getMass() { return (83.79584 - 3*eamu) * amuKg ; }

int KryptonIII::getCharge() { return 3 ; }

std::string KryptonIII::getFormula() { return "Kr+3" ; }

double KryptonIII::formationEnergy() { return 73.35945 * eVtoJ ; }

double KryptonIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* KryptonIII::Constituent() { return new Krypton; }

double Xenon::getMass() { return 131.293 * amuKg ; }

std::string Xenon::getFormula() { return "Xe" ; }

double Xenon::formationEnergy() { return 0.0 * eVtoJ ; }

double Xenon::IonLim() { return 12.1298 * eVtoJ ; }

double XenonI::getMass() { return (131.2923 - 1*eamu) * amuKg ; }

int XenonI::getCharge() { return 1 ; }

std::string XenonI::getFormula() { return "Xe+" ; }

double XenonI::formationEnergy() { return 12.1298 * eVtoJ ; }

double XenonI::IonLim() { return 21.20979 * eVtoJ ; }

Element* XenonI::Constituent() { return new Xenon; }

double XenonII::getMass() { return (131.2916 - 2*eamu) * amuKg ; }

int XenonII::getCharge() { return 2 ; }

std::string XenonII::getFormula() { return "Xe+2" ; }

double XenonII::formationEnergy() { return 33.33959 * eVtoJ ; }

double XenonII::IonLim() { return 32.123 * eVtoJ ; }

Element* XenonII::Constituent() { return new Xenon; }

double XenonIII::getMass() { return (131.2909 - 3*eamu) * amuKg ; }

int XenonIII::getCharge() { return 3 ; }

std::string XenonIII::getFormula() { return "Xe+3" ; }

double XenonIII::formationEnergy() { return 65.46259 * eVtoJ ; }

double XenonIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* XenonIII::Constituent() { return new Xenon; }

double Radon::getMass() { return 222.01758 * amuKg ; }

std::string Radon::getFormula() { return "Rn" ; }

double Radon::formationEnergy() { return 0.0 * eVtoJ ; }

double Radon::IonLim() { return 10.7485 * eVtoJ ; }

double RadonI::getMass() { return (222.0166 - 1*eamu) * amuKg ; }

int RadonI::getCharge() { return 1 ; }

std::string RadonI::getFormula() { return "Rn+" ; }

double RadonI::formationEnergy() { return 10.7485 * eVtoJ ; }

double RadonI::IonLim() { return 18.99 * eVtoJ ; }

Element* RadonI::Constituent() { return new Radon; }

double RadonII::getMass() { return (222.0156 - 2*eamu) * amuKg ; }

int RadonII::getCharge() { return 2 ; }

std::string RadonII::getFormula() { return "Rn+2" ; }

double RadonII::formationEnergy() { return 29.7385 * eVtoJ ; }

double RadonII::IonLim() { return 29.4 * eVtoJ ; }

Element* RadonII::Constituent() { return new Radon; }

double RadonIII::getMass() { return (222.0146 - 3*eamu) * amuKg ; }

int RadonIII::getCharge() { return 3 ; }

std::string RadonIII::getFormula() { return "Rn+3" ; }

double RadonIII::formationEnergy() { return 59.1385 * eVtoJ ; }

double RadonIII::IonLim() { return 9999999999.0 * eVtoJ ; }

Element* RadonIII::Constituent() { return new Radon; }
