# !/usr/bin/env python3
# PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni (University of Bologna, Italy)
# Licensed under CC BY 4.0. https://creativecommons.org/licenses/by/4.0/
#
# Python file to dynamically generate the AcceptedSpecies.h variant file.
# It always rewrites the header from scratch, deleting any previous version.

import re
import os
from pathlib import Path

# === PATHS (use the script directory, not CWD) ===
src_dir = Path(__file__).resolve().parent
main_dir = src_dir.parent
main_file = main_dir / "main.cpp"
output_file = src_dir / "AcceptedSpecies.h"
tmp_file = src_dir / "AcceptedSpecies.h.tmp"

# === COMPLETE SPACE-SEPARATED LIST OF CHEMICAL SPECIES IMPLEMENTED IN SPECIES.H ===
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
Krypton KryptonI KryptonII KryptonIII KryptonIV
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
Xenon XenonI XenonII XenonIII XenonIV
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
    comment = "/* !! NO_MIXTURE_MODE activated:\n   AcceptedSpecies generated but not used. */\n"
    print(" !! NO_MIXTURE_MODE: Dummy variant AcceptedSpecies generated.")
else:
    used_species = set(re.findall(r'new\s+(\w+)', main_code))
    variant_species = [s for s in all_species if s in used_species]
    print(f" [OK] AcceptedSpecies variant generated with {len(variant_species)} species used.")
    comment = ""

# === WRITE TO TEMP AND ATOMICALLY REPLACE PREVIOUS HEADER ===
# (This effectively "deletes" the previous file and rewrites it new each time.)
with open(tmp_file, "w") as f:
    f.write(
        "#ifndef ACCEPTED_SPECIES_H\n"
        "#define ACCEPTED_SPECIES_H\n\n"
        "/* FILE DYNAMICALLY GENERATED AT COMPILE-TIME,\n"
        "   CHECK AcceptedSpeciesGuard.py FOR DETAILS */\n\n"
        "// PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni\n"
        "// (University of Bologna, Italy)\n"
        "// Licensed under CC BY 4.0.\n"
        "// https://creativecommons.org/licenses/by/4.0/\n\n"
    )
    f.write('#include "Species.h"\n\n')
    f.write(comment)
    f.write("using AcceptedSpecies = std::variant<\n")
    for i, specie in enumerate(variant_species):
        comma = "*," if i < len(variant_species) - 1 else "*"
        f.write(f"    {specie}{comma}\n")
    f.write(">;\n\n#endif // ACCEPTED_SPECIES_H\n")

# Replace old file (if exists) with the new one
os.replace(tmp_file, output_file)
print(f" [OK] Wrote {output_file.name} (replaced previous if existed).")
