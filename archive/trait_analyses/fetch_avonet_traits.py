"""
Join the bird species list to AVONET (Tobias et al. 2022, Ecology Letters) to add
migration status and other ecological traits.

AVONET Migration codes:  1 = sedentary, 2 = partially migratory, 3 = migratory

Source: https://figshare.com/articles/dataset/16586228  (ELEData.zip -> AVONET1_BirdLife.csv)

Usage:
    python fetch_avonet_traits.py [INPUT_DIR]

Writes {INPUT_DIR}/species_traits_avonet.csv, one row per species.
"""
import csv, io, os, sys, urllib.request, zipfile

INPUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "/local/storage/kav67/clean_birds"
AVONET_URL = "https://ndownloader.figshare.com/files/38429873"
CACHE = os.path.join(INPUT_DIR, "AVONET1_BirdLife.csv")

# LatinName in our table -> AVONET BirdLife taxonomy.
# Synonyms, spelling variants, and subspecies trinomials that need collapsing.
SYNONYMS = {
    "Ammodramus caudacutus":            "Ammospiza caudacuta",
    "Anser cygnoides":                  "Anser cygnoid",
    "Balearica regulorum gibbericeps":  "Balearica regulorum",
    "Chroicocephalus ridibundus":       "Larus ridibundus",
    "Coloeus monedula":                 "Corvus monedula",
    "Dixiphia pipra":                   "Pseudopipra pipra",
    "Guaruba guaruba":                  "Guaruba guarouba",
    "Lonchura striata domestica":       "Lonchura striata",
    "Oenanthe melanoleuca":             "Oenanthe hispanica",
    "Psittacula echo":                  "Alexandrinus eques",
    "Strigops habroptilus":             "Strigops habroptila",
    "Struthio camelus australis":       "Struthio camelus",
    # BirdLife lumps both of these into A. californica; they are sedentary either way.
    "Aphelocoma insularis":             "Aphelocoma californica",
    "Aphelocoma woodhouseii":           "Aphelocoma californica",
}
LUMPED = {"Aphelocoma insularis", "Aphelocoma woodhouseii"}

MIGRATION_LABEL = {"1": "sedentary", "2": "partial", "3": "migratory"}

KEEP = ["Family1", "Order1", "Migration", "Habitat", "Habitat.Density", "Trophic.Level",
        "Trophic.Niche", "Primary.Lifestyle", "Mass", "Hand-Wing.Index", "Range.Size",
        "Min.Latitude", "Max.Latitude", "Centroid.Latitude", "Centroid.Longitude"]


def load_avonet():
    if not os.path.exists(CACHE):
        print(f"downloading AVONET -> {CACHE}")
        with urllib.request.urlopen(AVONET_URL) as r:
            blob = r.read()
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            data = z.read("ELEData/TraitData/AVONET1_BirdLife.csv")
        with open(CACHE, "wb") as f:
            f.write(data)
    with open(CACHE, encoding="utf-8-sig") as f:
        return {r["Species1"].strip().lower(): r for r in csv.DictReader(f)}


def main():
    avonet = load_avonet()

    species = {}
    with open(os.path.join(INPUT_DIR, "IGH_VGP_table.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            latin = " ".join(r["LatinName"].split())
            if latin and latin.lower() != "na":
                species[latin] = r["Species"]

    out_path = os.path.join(INPUT_DIR, "species_traits_avonet.csv")
    unmatched = []
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Species", "LatinName", "AvonetName", "MatchType",
                    "Migration", "MigrationLabel"] + KEEP[:2] + KEEP[3:])
        for latin, common in sorted(species.items()):
            target = SYNONYMS.get(latin, latin)
            match = "exact" if target == latin else ("lumped" if latin in LUMPED else "synonym")
            row = avonet.get(target.lower())
            if row is None:
                unmatched.append(latin)
                continue
            mig = row["Migration"].strip()
            w.writerow([common, latin, target, match, mig, MIGRATION_LABEL.get(mig, "")]
                       + [row["Family1"], row["Order1"]] + [row[k] for k in KEEP[3:]])

    print(f"wrote {out_path}  ({len(species) - len(unmatched)}/{len(species)} matched)")
    if unmatched:
        print("UNMATCHED:", *unmatched, sep="\n  ")


if __name__ == "__main__":
    main()
