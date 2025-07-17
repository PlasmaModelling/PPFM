#!/usr/bin/env python3
# Python file to dinamically generate the AcceptedSpecies.h variant file 
# This ensures minimum compile-time while having very large chemical Species 
# hard-coded datasets. 
# variant too big would saturate RAM

# if your implementing the main file without using Mixture,
# or any of the classes that depends on Mixture, 
# then write :" NO_MIXTURE " in your main file

# EVERY MODIFICATION TO THIS SCRIPT REQUIRE DELETING build FOLDER AND
# RE-BUILD FOR THE FIRST TIME


import re
import os
from pathlib import Path

# === PATHS ===
src_dir = Path(os.getcwd())
main_dir = src_dir.parent 
main_file = main_dir / "main.cpp"
output_file = src_dir / "AcceptedSpecies.h"

# === COMPLETE SPACE SEPARATED LIST OF CHEMICAL SPECIES IMPLEMENTED IN SPECIES.H ===
species_block = """
Electron
Aluminum AluminumI AluminumII AluminumIII
Argon ArgonI ArgonII ArgonIII ArgonIV
Arsenic ArsenicI ArsenicII ArsenicIII
Barium BariumI BariumII BariumIII
Beryllium BerylliumI BerylliumII BerylliumIII
Boron BoronI BoronII BoronIII
Bromine BromineI BromineII BromineIII
Calcium CalciumI CalciumII CalciumIII
Carbon CarbonI CarbonII CarbonIII
Cesium CesiumI CesiumII CesiumIII
Chlorine ChlorineI ChlorineII ChlorineIII
Chromium ChromiumI ChromiumII ChromiumIII
Cobalt CobaltI CobaltII CobaltIII
Copper CopperI CopperII CopperIII
Fluorine FluorineI FluorineII FluorineIII
Gallium GalliumI GalliumII GalliumIII
Germanium GermaniumI GermaniumII GermaniumIII
Gold GoldI GoldII GoldIII
Helium HeliumI HeliumII
Hydrogen HydrogenI
Iodine IodineI IodineII IodineIII
Iron IronI IronII IronIII
Krypton KryptonI KryptonII KryptonIII
Lead LeadI LeadII LeadIII
Lithium LithiumI LithiumII LithiumIII
Magnesium MagnesiumI MagnesiumII MagnesiumIII
Manganese ManganeseI ManganeseII ManganeseIII
Mercury MercuryI MercuryII MercuryIII
Molybdenum MolybdenumI MolybdenumII MolybdenumIII
Neon NeonI NeonII NeonIII
Nickel NickelI NickelII NickelIII
Niobium NiobiumI NiobiumII NiobiumIII
Nitrogen NitrogenI NitrogenII NitrogenIII
Oxygen OxygenI OxygenII OxygenIII
Phosphorus PhosphorusI PhosphorusII PhosphorusIII
Platinum PlatinumI PlatinumII PlatinumIII
Potassium PotassiumI PotassiumII PotassiumIII
Radon RadonI RadonII RadonIII
Rhenium RheniumI RheniumII RheniumIII
Rubidium RubidiumI RubidiumII RubidiumIII
Silicon SiliconI SiliconII SiliconIII
Silver SilverI SilverII SilverIII
Sodium SodiumI SodiumII SodiumIII
Strontium StrontiumI StrontiumII StrontiumIII
Sulfur SulfurI SulfurII SulfurIII
Tantalum TantalumI TantalumII TantalumIII
Titanium TitaniumI TitaniumII TitaniumIII
Tungsten TungstenI TungstenII TungstenIII
Vanadium VanadiumI VanadiumII VanadiumIII
Xenon XenonI XenonII XenonIII
Zinc ZincI ZincII ZincIII
Zirconium ZirconiumI ZirconiumII ZirconiumIII
HydrogenAnion NitrogenAnion OxygenAnion
MolecularHydrogen MolecularHydrogenI MolecularNitrogen
MolecularNitrogenI MolecularOxygen MolecularOxygenAnion
MolecularOxygenI Ozone
"""

# === PARSING ===
all_species = species_block.split()

# === READ USED SPECIES IN main.cpp ===
with open(main_file, "r") as f:
    main_code = f.read()

# === CHECKS FOR NO MIXTURE FLAG ===
no_mixture_mode = "NO_MIXTURE" in main_code

# === MINIMAL VARIANT ===
if no_mixture_mode:
    variant_species = ["Argon", "ArgonI"]
    comment = "/* !! NO_MIXTURE_MODE activated: \n AcceptedSpecies generated but not used. */ \n"
    print(" !! NO_MIXTURE_MODE: Dummy variant AcceptedSpecies generated.")
else:
    used_species = set(re.findall(r'new\s+(\w+)', main_code))
    variant_species = [s for s in all_species if s in used_species]
    print(f" [OK] AcceptedSpecies variant generated with {len(variant_species)} species used.")
    comment = ""

# === ACCEPTEDSPECIES HEADER FILE WRITING ===
with open(output_file, "w") as f:
    f.write("#ifndef ACCEPTED_SPECIES_H\n"
    "#define ACCEPTED_SPECIES_H\n\n"
    "/* FILE DINAMICALLY GENERATED AT COMPILE-TIME,\n CHECK AcceptedSpeciesGuard.py FOR DETAILS */\n\n")
    f.write('#include "Species.h"\n\n')
    f.write(comment)
    f.write("using AcceptedSpecies = std::variant<\n")

    for i, specie in enumerate(variant_species):
        comma = "*," if i < len(variant_species) - 1 else "*"
        f.write(f"    {specie}{comma}\n")

    f.write(">;\n\n#endif // ACCEPTED_SPECIES_H\n")
