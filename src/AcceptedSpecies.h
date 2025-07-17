#ifndef ACCEPTED_SPECIES_H
#define ACCEPTED_SPECIES_H

/* FILE DINAMICALLY GENERATED AT COMPILE-TIME,
 CHECK AcceptedSpeciesGuard.py FOR DETAILS */

#include "Species.h"

using AcceptedSpecies = std::variant<
    Electron*,
    Argon*,
    ArgonI*,
    ArgonII*
>;

#endif // ACCEPTED_SPECIES_H
