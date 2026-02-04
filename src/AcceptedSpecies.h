#ifndef ACCEPTED_SPECIES_H
#define ACCEPTED_SPECIES_H

/* FILE DYNAMICALLY GENERATED AT COMPILE-TIME,
   CHECK AcceptedSpeciesGuard.py FOR DETAILS */

// PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni
// (University of Bologna, Italy)
// Licensed under CC BY 4.0.
// https://creativecommons.org/licenses/by/4.0/

#include "Species.h"

using AcceptedSpecies = std::variant<
    Electron*,
    Carbon*,
    CarbonI*,
    CarbonII*,
    Hydrogen*,
    HydrogenI*,
    Oxygen*,
    OxygenI*,
    OxygenII*,
    HydrogenAnion*,
    OxygenAnion*,
    MolecularHydrogen*,
    MolecularHydrogenI*,
    MolecularOxygen*,
    MolecularOxygenI*,
    CarbonMonoxide*,
    CarbonDioxide*,
    CarbonMonoxideI*,
    CarbonDioxideI*,
    CarbonAnion*,
    DiCarbon*,
    DiCarbonI*,
    EthynilRadical*,
    Acetylene*,
    Ethylene*,
    Formaldehyde*,
    DiCarbonMonoxide*,
    TriCarbon*,
    MethylidineRadical*,
    Methylene*,
    MethylRadical*,
    Methane*,
    FormylRadical*,
    FormylRadicalI*,
    Water*,
    HydroPeroxyRadical*,
    HydroxilRadical*,
    HydroxilRadicalI*
>;

#endif // ACCEPTED_SPECIES_H
