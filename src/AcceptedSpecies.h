#ifndef ACCEPTED_SPECIES_H
#define ACCEPTED_SPECIES_H

/* FILE DINAMICALLY GENERATED AT COMPILE-TIME,
 CHECK AcceptedSpeciesGuard.py FOR DETAILS */

// PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni // 
// (University of Bologna, Italy) // 
// Licensed under CC BY 4.0. // 
// To view a copy of this license, visit: // 
// https://creativecommons.org/licenses/by/4.0/ // 

#include "Species.h"

using AcceptedSpecies = std::variant<
    Electron*,
    Argon*,
    ArgonI*,
    ArgonII*
>;

#endif // ACCEPTED_SPECIES_H
