"""
Classify each haplotype's data source (VGP / CCGP / one of the collaborator
pangenome projects / unpublished) and its haplotype type (maternal / paternal /
primary-alternate pseudohaplotype / Hi-C-phased hap1-hap2 / single assembly).

Source is resolved by querying the NCBI Datasets API for the assembly's
BioProject lineage, cross-referenced against /local/storage/kav67/VGP_details.csv
(the VGP master tracking sheet) by accession. That cross-reference matters:
NCBI's own lineage text chains *every* genome in the sheet through "Vertebrate
Genomes Project" regardless of which initiative actually produced it, so text
matching alone can't tell a true VGP submission apart from one of the other
sequencing initiatives the sheet also tracks (Darwin Tree of Life, AmaZoomics,
Sanger 25G, Revive & Restore, Colossal, etc.) -- only the sheet's own "Main
project" column can.

A subset of species are California Conservation Genomics Project (CCGP)
genomes -- for a few (e.g. Savannah Sparrow, Song Sparrow) we have *both* a
VGP-tracked assembly and a separate CCGP assembly of the same species,
distinguished only by directory name (the "2" suffix marks whichever was
added second, not the project). House finches, jays (Aphelocoma scrub-jays +
Yucatan Jay), and Capuchino seedeaters are multi-individual pangenome
datasets from three different collaborators, not single VGP/CCGP submissions.
A handful of Darwin's finches have no NCBI accession at all and are
unpublished.

Haplotype type is resolved from the assembly filename itself (VGP/NCBI
filenames embed .mat/.pat for trio-phased assemblies, .hap1/.hap2 for
Hi-C-only phasing, .pri/.alt for solo pseudohaplotypes -- these mean different
things and are not interchangeable). A handful of VGP assemblies were saved
locally under a generic "{tolid}_pri.fna"/"{tolid}_alt.fna" name that doesn't
preserve the original suffix (this can be misleading: e.g. Eurasian Collared
Dove's real assembly is Hi-C-phased hap1/hap2, not a pri/alt pseudohaplotype
pair, and Rock Dove is trio-phased with the "main"/pri-equivalent haplotype
actually being paternal). For those, and for any other case the filename
doesn't resolve, we cross-reference /local/storage/kav67/VGP_details.csv
(the VGP master species sheet) by scientific name -- only for the haplotype
suffix, not any of its other columns.

Usage:
    python assign_haplotype_source.py [INPUT_DIR]

Reads {INPUT_DIR}/bird_genome_paths.csv (Haplotype = "Order/Species/Haplotype", Genome = path).
Writes {INPUT_DIR}/haplotype_sources.csv with columns:
    Order, Species, Haplotype, Accession, Source, SourceDetail, HaplotypeType, Published

Caches raw NCBI responses in {INPUT_DIR}/ncbi_bioproject_cache.json so re-runs
(e.g. after adding new haplotypes) don't re-query already-seen accessions.
"""
import csv, json, os, re, subprocess, sys, time
from collections import defaultdict

INPUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "/local/storage/kav67/clean_birds"
GENOME_PATHS = os.path.join(INPUT_DIR, "bird_genome_paths.csv")
CACHE_PATH = os.path.join(INPUT_DIR, "ncbi_bioproject_cache.json")
OUT_PATH = os.path.join(INPUT_DIR, "haplotype_sources.csv")
VGP_DETAILS_PATH = "/local/storage/kav67/VGP_details.csv"

ACC_RE = re.compile(r"(GC[AF]_\d+\.\d+)")
ACC_PREFIX_RE = re.compile(r"^GC[AF]_\d+\.\d+_?")

# One-off assemblies whose NCBI BioProject title doesn't cleanly say VGP/CCGP/
# pangenome, resolved manually by inspecting bioproject_lineage.
EXCEPTIONS = {
    "GCA_000002315.5": ("Other (GRC reference)", "Chicken Genome Reference Consortium (GRCg6a)"),
    "GCA_902806625.1": ("Published (non-VGP/CCGP)", "Rubin lab, Uppsala -- Nanopore Camarhynchus parvulus assembly"),
    "GCF_047663525.1": ("Other (non-VGP RefSeq)", "Pekin duck T2T reference (CAAS)"),
    "GCA_040182565.1": ("Other (non-VGP GenBank)", "Anser cygnoides breed:goose (Jiangsu Agri-animal Husbandry)"),
    "GCA_054553165.1": ("Other (non-VGP)", "Philippine Eagle trio-binned ONT assembly (Cal Academy of Sciences)"),
}

# Species with no LatinName resolvable elsewhere (VGP doves saved locally
# under a generic filename with no accession) -- needed only to look them up
# in VGP_details.csv by scientific name.
MANUAL_LATIN = {
    "Eurasian_CollaredDove": "Streptopelia decaocto",
    "Pink_Pigeon": "Nesoenas mayeri",
    "Nicobar_Pigeon": "Caloenas nicobarica",
    "BandTailed_Pigeon": "Patagioenas fasciata",
    "Rock_Dove": "Columba livia",
    "European_Turtle_Dove": "Streptopelia turtur",
}

SUFFIX_LABEL = {
    "mat": "Maternal",
    "pat": "Paternal",
    "hap1": "Haplotype 1",
    "hap2": "Haplotype 2",
    "pri": "Primary",
    "alt": "Alternate",
    "p": "Primary",
    "a": "Alternate",
    "merged": "Merged",
}

# VGP_details.csv "Main project" values, shortened where they carry a long
# descriptive tail after " :".
def _short_project(name):
    return name.split(":")[0].strip() if name else name
COMPLEMENT = {"pri": "alt", "alt": "pri", "hap1": "hap2", "hap2": "hap1",
              "mat": "pat", "pat": "mat", "p": "a", "a": "p"}


def load_cache():
    if os.path.exists(CACHE_PATH):
        return json.load(open(CACHE_PATH))
    return {}


def fetch_lineage(accessions, cache):
    missing = [a for a in accessions if a not in cache]
    chunk = 25
    for i in range(0, len(missing), chunk):
        batch = missing[i:i + chunk]
        url = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
               + ",".join(batch) + "/dataset_report?page_size=100")
        out = subprocess.run(["curl", "-s", "--max-time", "60", url],
                              capture_output=True, text=True).stdout
        try:
            reports = json.loads(out).get("reports", [])
        except Exception:
            reports = []
        for rep in reports:
            acc = rep["accession"]
            ai = rep.get("assembly_info", {})
            titles = [bp.get("title") for lg in ai.get("bioproject_lineage", [])
                      for bp in lg.get("bioprojects", [])]
            cache[acc] = {"organism": rep.get("organism", {}).get("organism_name"),
                          "submitter": ai.get("submitter"), "lineage": titles}
        time.sleep(0.5)
    json.dump(cache, open(CACHE_PATH, "w"))
    return cache


def classify_source(acc, path, cache, species, by_species, by_accession, latin_lookup):
    """Returns (Source, SourceDetail) -- detail is empty except for one-off exceptions.

    CCGP and the three collaborator pangenome projects are identified from
    NCBI's own BioProject title text (unambiguous, and not tracked in
    VGP_details.csv at all). Everything else that traces back to the VGP
    master tracking sheet (VGP proper, but also Darwin Tree of Life,
    AmaZoomics, and other affiliated sequencing initiatives it also tracks)
    is resolved from that sheet's "Main project" column via accession -- NCBI's
    lineage chains all of these through "Vertebrate Genomes Project" regardless
    of which initiative actually produced the assembly, so that text alone
    can't tell them apart.
    """
    if acc and acc in EXCEPTIONS:
        return EXCEPTIONS[acc]
    if acc and acc in cache:
        info = cache[acc]
        titles = " | ".join(t for t in info["lineage"] if t)
        tl = titles.lower()
        if "california conservation genomics project" in tl:
            return "CCGP", ""
        if "scrub-jay" in tl and "pangenome" in tl:
            return "Jay pangenome (Edwards lab)", ""
        if "house finch assembly" in tl:
            return "House finch pangenome (Edwards lab)", ""
        if "capuchino_seedeaters" in path.lower() or "sporophila" in (info["organism"] or "").lower():
            return "Seedeater pangenome (Campagna lab)", ""
    if acc and acc in by_accession:
        return by_accession[acc], ""
    sci_name = (cache.get(acc) or {}).get("organism") or latin_lookup.get(species)
    if sci_name and sci_name in by_species:
        return by_species[sci_name]["project"], ""
    if acc and acc in cache and "vertebrate genomes project" in " ".join(cache[acc]["lineage"]).lower():
        return "VGP", ""
    if "#Kats_Birds" in path and "/finches/" in path:
        return "Unpublished Darwin finch", ""
    if "RedwingedBlackbird_new" in path:
        return "Unpublished / Pennell Lab assembly", ""
    if not acc and "#Birds" in path:
        return "VGP", ""
    if acc and acc in cache:
        return "UNCLASSIFIED -- needs manual review", " | ".join(t for t in cache[acc]["lineage"] if t)
    return "UNRESOLVED -- accession not found on NCBI", ""


def local_suffix_token(genome_path):
    b = os.path.basename(genome_path)
    for ext in (".fna", ".fa", ".fasta"):
        if b.lower().endswith(ext):
            b = b[:-len(ext)]
            break
    b = ACC_PREFIX_RE.sub("", b)
    b = re.sub(r"_genomic$", "", b, flags=re.I).lower()
    if re.search(r"haplotype[_.]?1|hap1\b", b):
        return "hap1"
    if re.search(r"haplotype[_.]?2|hap2\b", b):
        return "hap2"
    if re.search(r"(?<![a-z])mat(?![a-z])", b):
        return "mat"
    if re.search(r"(?<![a-z])pat[w]?(?![a-z])", b):
        return "pat"
    if "merged" in b:
        return "merged"
    if "alternate" in b or re.search(r"(?<![a-z])alt(?![a-z])", b):
        return "alt"
    if "primary" in b or "principal" in b or re.search(r"(?<![a-z])pri(?![a-z])", b):
        return "pri"
    if re.search(r"(^|[._])p([._]|$)", b):
        return "p"
    if re.search(r"(^|[._])a([._]|$)", b):
        return "a"
    return None


def load_vgp_details():
    """Returns (by_species, by_accession):
    by_species: LatinName -> {"main": suffix, "second": suffix, "project": Main project}
    by_accession: accession -> Main project (for every accession the sheet lists for
    that species, main or alternate haplotype alike -- used only to pin down which
    sequencing initiative (VGP / DToL / etc.) actually produced the assembly, since
    NCBI's own BioProject lineage chains all of these through "Vertebrate Genomes
    Project" regardless of which one did the work).
    """
    if not os.path.exists(VGP_DETAILS_PATH):
        return {}, {}

    def suffix(assembly_id):
        assembly_id = assembly_id.strip().split(",")[0].strip()  # first if multiple listed
        return assembly_id.split(".")[-1].lower() if assembly_id else None

    by_species, by_accession = {}, {}
    with open(VGP_DETAILS_PATH, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            if row.get("Lineage") != "Birds":
                continue
            project = _short_project(row["Main project"].strip()) or "VGP"
            main = suffix(row["Assembly ID main haplotype"])
            second = suffix(row["Assembly IDs other high-quality haplotypes"]) or suffix(row["Assembly IDs alternate haplotypes"])
            by_species[row["Scientific Name"].strip()] = {"main": main, "second": second, "project": project}
            blob = " ".join([row["Accession # for main haplotype"],
                              row["Accession #s other high-quality haplotypes"],
                              row["Accession #s alternate haplotypes"]])
            for acc in ACC_RE.findall(blob):
                by_accession[acc] = project
    return by_species, by_accession


def pair_key(hap_label):
    """Strip a trailing haplotype-suffix word so sibling haplotypes of the same
    individual (e.g. bAytFul2_pri / bAytFul2_alt) group together, while a
    lone assembly with an unrelated tolid (e.g. bAytFul3) stays its own group."""
    return re.sub(r"[._]?(pri|alt|hap1|hap2|mat|pat|principal|alternate)$", "", hap_label, flags=re.I)


def classify_haplotype_type(order, species, hap_label, acc, genome_path, source, vgp_by_species, latin_lookup):
    # Prefer the real filename suffix; only fall back to VGP_details.csv when
    # the local file was saved under a generic "{tolid}_pri/alt" name with no
    # accession (that name can't be trusted -- see module docstring).
    if acc:
        return local_suffix_token(genome_path)
    sci_name = latin_lookup.get(species)
    info = vgp_by_species.get(sci_name) if sci_name else None
    if not info:
        return local_suffix_token(genome_path)
    is_pri_row = hap_label.lower().endswith(("_pri", ".pri"))
    return info["main"] if is_pri_row else info["second"]


def main():
    rows = list(csv.DictReader(open(GENOME_PATHS)))
    recs = []
    for r in rows:
        parts = r["Haplotype"].split("/")
        order, species, hap = parts[0], parts[1], "/".join(parts[2:])
        m = ACC_RE.search(r["Genome"])
        recs.append({"Order": order, "Species": species, "Haplotype": hap,
                     "Accession": m.group(1) if m else "", "Genome": r["Genome"]})

    accs = sorted({r["Accession"] for r in recs if r["Accession"]})
    cache = load_cache()
    cache = fetch_lineage(accs, cache)
    vgp_by_species, vgp_by_accession = load_vgp_details()

    for r in recs:
        r["Source"], r["SourceDetail"] = classify_source(
            r["Accession"], r["Genome"], cache, r["Species"], vgp_by_species, vgp_by_accession, MANUAL_LATIN)

    # sibling-pair groups, for complementary inference (pri<->alt, hap1<->hap2, mat<->pat)
    groups = defaultdict(list)
    for r in recs:
        groups[(r["Order"], r["Species"], pair_key(r["Haplotype"]))].append(r)

    for r in recs:
        r["_token"] = classify_haplotype_type(r["Order"], r["Species"], r["Haplotype"],
                                               r["Accession"], r["Genome"], r["Source"],
                                               vgp_by_species, MANUAL_LATIN)

    for key, members in groups.items():
        if len(members) != 2:
            continue
        toks = [m["_token"] for m in members]
        if toks[0] and not toks[1] and toks[0] in COMPLEMENT:
            members[1]["_token"] = COMPLEMENT[toks[0]]
        elif toks[1] and not toks[0] and toks[1] in COMPLEMENT:
            members[0]["_token"] = COMPLEMENT[toks[1]]

    for key, members in groups.items():
        for r in members:
            if r["_token"]:
                r["HaplotypeType"] = SUFFIX_LABEL.get(r["_token"], r["_token"])
            elif len(members) == 1:
                r["HaplotypeType"] = "Single assembly (no haplotype pair)"
            else:
                r["HaplotypeType"] = "Unresolved -- check Genome path manually"

    with open(OUT_PATH, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Order", "Species", "Haplotype", "Accession", "Source", "SourceDetail",
                    "HaplotypeType", "Published"])
        for r in recs:
            published = "FALSE" if r["Source"].startswith("Unpublished") else "TRUE"
            w.writerow([r["Order"], r["Species"], r["Haplotype"], r["Accession"], r["Source"],
                        r["SourceDetail"], r["HaplotypeType"], published])

    print(f"wrote {OUT_PATH} ({len(recs)} haplotypes)")
    flagged_source = [r for r in recs if r["Source"].startswith(("UNCLASSIFIED", "UNRESOLVED"))]
    flagged_type = [r for r in recs if r["HaplotypeType"] == "Unresolved -- check Genome path manually"]
    if flagged_source:
        print("SOURCE NEEDS MANUAL REVIEW:")
        for r in flagged_source:
            print(" ", r["Order"], r["Species"], r["Haplotype"], r["Accession"])
    if flagged_type:
        print("HAPLOTYPE TYPE NEEDS MANUAL REVIEW:")
        for r in flagged_type:
            print(" ", r["Order"], r["Species"], r["Haplotype"], r["Genome"])


if __name__ == "__main__":
    main()
