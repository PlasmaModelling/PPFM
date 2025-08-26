#ifndef ACCEPTED_SPECIES_H
#define ACCEPTED_SPECIES_H

/* FILE DINAMICALLY GENERATED AT COMPILE-TIME,
 CHECK AcceptedSpeciesGuard.py FOR DETAILS */

#include "Species.h"

/* !! NO_MIXTURE_MODE activated: 
 AcceptedSpecies generated but not used. */ 
using AcceptedSpecies = std::variant<
    Argon*,
    ArgonI*
>;

#endif // ACCEPTED_SPECIES_H
