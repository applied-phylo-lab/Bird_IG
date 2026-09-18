"""
Build the publication-facing summary tables for the manuscript.

Combines summary_features.csv, gene_list.csv, inversion_stats.tsv,
D_inversions.tsv, palindromes.tsv, and haplotype_sources.csv (see
assign_haplotype_source.py) into a small set of clean, documented tables
under {INPUT_DIR}/summary_tables/.

Inversion/D-gene/palindrome analyses are IGH-only in this pipeline (see
README.md), so those columns are populated only for IGH rows.

A haplotype's IGH (or IGL) locus is sometimes split across multiple contigs
in summary_features.csv (one row per contig). This script aggregates to one
row per Order/Species/Haplotype/Locus, summing gene counts across contigs and
recomputing fractions from the summed numerator/denominator (never averaging
fractions directly).

Usage:
    python build_summary_tables.py [INPUT_DIR]

Writes, each as an "all haplotypes" version and a "_published" version that
drops the unpublished Darwin's finches (EXCLUDE_FROM_PUBLISHED, above):
    {INPUT_DIR}/summary_tables/main_table[_published].csv        -- one row per haplotype x locus
    {INPUT_DIR}/summary_tables/species_summary[_published].csv  -- one row per species (mean/range over haplotypes)
    {INPUT_DIR}/summary_tables/v_gene_table[_published].csv     -- one row per called V gene
    {INPUT_DIR}/summary_tables/d_gene_table[_published].csv     -- one row per called D gene (IGH only)
    {INPUT_DIR}/summary_tables/README.md                        -- column documentation
"""
import csv, json, os, sys
from collections import defaultdict

INPUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "/local/storage/kav67/clean_birds"
OUT_DIR = os.path.join(INPUT_DIR, "summary_tables")
INV_MINLEN = "250"  # the only threshold summarize_inversions.py emits; matches D_inversions.tsv default

# Sources excluded from the *_published.csv variants. Only Darwin's finches so
# far (explicitly asked for) -- note the one Red-winged Blackbird assembly is
# also Source-tagged "Unpublished" but is deliberately NOT in this set unless
# asked, since it wasn't part of the request that created this filter.
EXCLUDE_FROM_PUBLISHED = {"Unpublished Darwin finch"}

# Species with no LatinName in IGH_VGP_table.tsv / species_traits_avonet.csv / NCBI cache
# (mainly VGP doves whose local filenames don't embed a GCA accession).
MANUAL_LATIN = {
    "Eurasian_CollaredDove": "Streptopelia decaocto",
    "Pink_Pigeon": "Nesoenas mayeri",
    "Nicobar_Pigeon": "Caloenas nicobarica",
    "BandTailed_Pigeon": "Patagioenas fasciata",
    "Rock_Dove": "Columba livia",
    "European_Turtle_Dove": "Streptopelia turtur",
    # unpublished Darwin's finches, no NCBI accession to pull from
    "Common_cactus_finch": "Geospiza scandens",
    "Large_ground_finch": "Geospiza magnirostris",
    "Medium_ground_finch": "Geospiza fortis",
    "Small_ground_finch": "Geospiza fuliginosa",
}


def path(*p):
    return os.path.join(INPUT_DIR, *p)


def write_table(rows, cols, name, source_field="Source"):
    """Writes {name}.csv (all rows) and {name}_published.csv (rows whose
    Source isn't in EXCLUDE_FROM_PUBLISHED -- for species_summary.csv rows,
    source_field is "Sources" and holds a ";"-joined list instead of one value)."""
    full_path = os.path.join(OUT_DIR, f"{name}.csv")
    with open(full_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {full_path} ({len(rows)} rows)")

    def is_published(row):
        val = row.get(source_field, "")
        srcs = val.split(";") if source_field == "Sources" else [val]
        return not any(s in EXCLUDE_FROM_PUBLISHED for s in srcs)

    pub_rows = [r for r in rows if is_published(r)]
    pub_path = os.path.join(OUT_DIR, f"{name}_published.csv")
    with open(pub_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, restval="")
        w.writeheader()
        w.writerows(pub_rows)
    print(f"wrote {pub_path} ({len(pub_rows)} rows, {len(rows) - len(pub_rows)} excluded)")


def load_latin_names():
    latin = {}
    with open(path("IGH_VGP_table.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["LatinName"]:
                latin[row["Species"]] = row["LatinName"]
    with open(path("species_traits_avonet.csv")) as f:
        for row in csv.DictReader(f):
            latin.setdefault(row["Species"], row["LatinName"])

    cache_path = path("ncbi_bioproject_cache.json")
    if os.path.exists(cache_path):
        cache = json.load(open(cache_path))
        with open(path("haplotype_sources.csv")) as f:
            for row in csv.DictReader(f):
                sp = row["Species"]
                if sp in latin:
                    continue
                info = cache.get(row["Accession"])
                if info and info.get("organism"):
                    latin[sp] = info["organism"]

    for sp, name in MANUAL_LATIN.items():
        latin.setdefault(sp, name)
    return latin


def load_sources():
    src = {}
    with open(path("haplotype_sources.csv")) as f:
        for row in csv.DictReader(f):
            key = (row["Order"], row["Species"], row["Haplotype"])
            src[key] = (row["Source"], row["SourceDetail"], row["HaplotypeType"], row["Published"])
    return src


def load_gene_counts():
    """(Order/Species/Haplotype, Locus) -> dict of gene-level counts, from gene_list.csv."""
    counts = defaultdict(lambda: {"NumV": 0, "NumV_productive": 0, "NumV_with_RSS": 0})
    with open(path("gene_list.csv")) as f:
        for row in csv.DictReader(f):
            if row["Locus"] not in ("IGH", "IGL"):
                continue
            if row["Passes Filtering"] != "True":
                continue
            parts = row["Source"].split("/")
            key = (parts[0], parts[1], "/".join(parts[2:]), row["Locus"])
            c = counts[key]
            c["NumV"] += 1
            if row["Productive"] == "True":
                c["NumV_productive"] += 1
            if row["Heptamer"].strip() or row["Nonamer"].strip():
                c["NumV_with_RSS"] += 1
    return counts


def load_inversion_stats():
    """(Order,Species,Haplotype) -> summed IGH inversion stats at INV_MINLEN, across contigs."""
    agg = defaultdict(lambda: {"NumInversions": 0, "NumInversions_diag": 0,
                                "TotalSeqLength_bp": 0, "InvCoverage_bp": 0,
                                "TotalGenes": 0, "GenesOnInv": 0})
    with open(path("inversion_stats.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["minlen"] != INV_MINLEN:
                continue
            key = (row["order"], row["species"], row["haplotype"])
            a = agg[key]
            a["NumInversions"] += int(row["num_inversions"])
            a["NumInversions_diag"] += int(row["num_inversions_diag"])
            a["TotalSeqLength_bp"] += int(row["total_seq_length"])
            a["InvCoverage_bp"] += int(row["inv_cov_len"])
            a["TotalGenes"] += int(row["total_genes"])
            a["GenesOnInv"] += int(row["genes_on_inv"])
    return agg


def load_d_inversions():
    """D gene / inversion overlap, if the table is present.

    D_inversions.tsv was retired to {INPUT_DIR}/old_unused/ -- the upstream IGHD.csv
    files became a threshold-swept candidate list and the unfiltered counts are no
    longer meaningful (see that folder's README). Absent means the NumD columns are
    simply left blank rather than the build failing."""
    d = {}
    d_path = path("D_inversions.tsv")
    if not os.path.exists(d_path):
        print("note: D_inversions.tsv not present -- NumD columns left blank")
        return d
    with open(d_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = (row["Order"], row["Species"], row["Haplotype"])
            d[key] = {"NumD": int(row["n_d_genes"]), "NumD_on_inv": int(row["n_on_inversion"])}
    return d


def load_palindromes():
    agg = defaultdict(lambda: {"n": 0, "whole": 0.0, "mid20": 0.0, "rand20": 0.0, "delta": 0.0})
    with open(path("palindromes.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = (row["Order"], row["Species"], row["Haplotype"])
            a = agg[key]
            a["n"] += 1
            a["whole"] += float(row["WholeIdentity"])
            a["mid20"] += float(row["MiddleIdentity20bp"])
            a["rand20"] += float(row["RandomIdentity20bp"])
            a["delta"] += float(row["MiddleMinusWhole"])
    return agg


V_GENE_COLS = ["Order", "Species", "Haplotype", "Locus", "LatinName", "Source", "HaplotypeType",
               "Published", "Contig", "Pos", "Strand", "Productive", "PassesFiltering",
               "Heptamer", "Nonamer", "RSS_DownstreamBp", "Sequence"]
D_GENE_COLS = ["Order", "Species", "Haplotype", "Locus", "LatinName", "Source", "HaplotypeType",
               "Published", "Contig", "Pos", "Strand", "Productive", "UpstreamHeptamer",
               "UpstreamNonamer", "DownstreamHeptamer", "DownstreamNonamer",
               "LocationRelativeToVCluster", "Sequence"]


def _hap_context(source_field, sources, latin):
    parts = source_field.split("/")
    order, species, hap = parts[0], parts[1], "/".join(parts[2:])
    skey = (order, species, hap)
    source, _detail, hap_type, published = sources.get(skey, ("", "", "", ""))
    return {"Order": order, "Species": species, "Haplotype": hap,
            "LatinName": latin.get(species, ""), "Source": source,
            "HaplotypeType": hap_type, "Published": published}


def build_v_gene_table(sources, latin, indexed):
    """One row per called V gene from gene_list.csv, with haplotype-level
    context (species, source, haplotype type) joined in so the table is
    usable on its own. TRA/TRB/TRG/TRD rows in gene_list.csv are dropped --
    this project covers IGH/IGL only (see CLAUDE.md); everything else in
    there is leftover annotation from a shared upstream pipeline."""
    rows = []
    skipped = 0
    with open(path("gene_list.csv")) as f:
        for row in csv.DictReader(f):
            if row["Locus"] not in ("IGH", "IGL"):
                continue
            parts = row["Source"].split("/")
            if (parts[0], parts[1], "/".join(parts[2:])) not in indexed:
                skipped += 1
                continue
            rows.append({
                **_hap_context(row["Source"], sources, latin),
                "Locus": row["Locus"], "Contig": row["Contig"], "Pos": row["Pos"],
                "Strand": row["Strand"], "Productive": row["Productive"],
                "PassesFiltering": row["Passes Filtering"], "Heptamer": row["Heptamer"],
                "Nonamer": row["Nonamer"], "RSS_DownstreamBp": row["Number of bp Downstream (RSS)"],
                "Sequence": row["Sequence"],
            })
    print(f"v_gene_table: skipped {skipped} genes on haplotypes not in summary_features.csv")
    write_table(rows, V_GENE_COLS, "v_gene_table")


def build_d_gene_table(sources, latin, indexed):
    """One row per called D gene from bird_d_genes.csv (IGH only in this
    pipeline). Separate from v_gene_table.csv because V and D genes carry
    different RSS annotation (V: one downstream RSS; D: RSS on both sides)."""
    rows = []
    skipped = 0
    with open(path("bird_d_genes.csv")) as f:
        for row in csv.DictReader(f):
            parts = row["Source"].split("/")
            if (parts[0], parts[1], "/".join(parts[2:])) not in indexed:
                skipped += 1
                continue
            rows.append({
                **_hap_context(row["Source"], sources, latin),
                "Locus": row["Locus"], "Contig": row["Contig"], "Pos": row["Pos"],
                "Strand": row["Strand"], "Productive": row["Productive"],
                "UpstreamHeptamer": row["Upstream Heptamer"], "UpstreamNonamer": row["Upstream Nonamer"],
                "DownstreamHeptamer": row["Downstream Heptamer"], "DownstreamNonamer": row["Downstream Nonamer"],
                "LocationRelativeToVCluster": row["Location Relative to V-Cluster"],
                "Sequence": row["Sequence"],
            })
    print(f"d_gene_table: skipped {skipped} genes on haplotypes not in summary_features.csv")
    write_table(rows, D_GENE_COLS, "d_gene_table")


def load_indexed_haplotypes():
    """(Order, Species, Haplotype) present in summary_features.csv.

    gene_list.csv and bird_d_genes.csv come from the upstream pipeline and cover
    every haplotype ever annotated, including ones deliberately kept out of the
    cross-species analysis (config/excluded_haplotypes.csv). Without this filter
    main_table.csv honours the exclusions while v_gene_table.csv does not, so the
    two disagree about which assemblies exist."""
    haps = set()
    with open(path("summary_features.csv")) as f:
        for row in csv.DictReader(f):
            haps.add((row["Order"], row["Species"], row["Haplotype"]))
    return haps


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    latin = load_latin_names()
    sources = load_sources()
    gene_counts = load_gene_counts()
    inv_stats = load_inversion_stats()
    d_inv = load_d_inversions()
    palin = load_palindromes()

    # summary_features.csv is one row per contig; aggregate to one row per
    # Order/Species/Haplotype/Locus, summing NumV and collecting contig list.
    hap_locus = defaultdict(lambda: {"NumV_sf": 0, "Contigs": []})
    with open(path("summary_features.csv")) as f:
        for row in csv.DictReader(f):
            key = (row["Order"], row["Species"], row["Haplotype"], row["Locus"])
            hap_locus[key]["NumV_sf"] += int(row["NumV"])
            hap_locus[key]["Contigs"].append(row["Contig"])

    main_rows = []
    for (order, species, hap, locus), sf in sorted(hap_locus.items()):
        gkey = (order, species, hap, locus)
        gc = gene_counts.get(gkey, {"NumV": sf["NumV_sf"], "NumV_productive": "", "NumV_with_RSS": ""})
        skey = (order, species, hap)
        source, source_detail, hap_type, published = sources.get(skey, ("", "", "", ""))

        row = {
            "Order": order, "Species": species, "Haplotype": hap, "Locus": locus,
            "LatinName": latin.get(species, ""),
            "Source": source, "SourceDetail": source_detail, "HaplotypeType": hap_type,
            "Published": published,
            "NumContigs": len(sf["Contigs"]), "Contigs": ";".join(sf["Contigs"]),
            "NumV": gc["NumV"],
            "NumV_productive": gc["NumV_productive"],
            "NumV_with_RSS": gc["NumV_with_RSS"],
            "FracV_with_RSS": (round(gc["NumV_with_RSS"] / gc["NumV"], 4)
                               if gc["NumV_with_RSS"] != "" and gc["NumV"] else ""),
        }

        if locus == "IGH":
            inv = inv_stats.get(skey)
            if inv:
                row.update({
                    f"NumInversions_min{INV_MINLEN}bp": inv["NumInversions"],
                    f"NumInversionsDiag_min{INV_MINLEN}bp": inv["NumInversions_diag"],
                    "LocusLength_bp": inv["TotalSeqLength_bp"],
                    "InvCoverage_bp": inv["InvCoverage_bp"],
                    "FracGenesOnInv": (round(inv["GenesOnInv"] / inv["TotalGenes"], 4)
                                       if inv["TotalGenes"] else ""),
                })
            d = d_inv.get(skey)
            if d:
                row.update({
                    "NumD": d["NumD"], "NumD_on_inv": d["NumD_on_inv"],
                    "FracD_on_inv": round(d["NumD_on_inv"] / d["NumD"], 4) if d["NumD"] else "",
                })
            p = palin.get(skey)
            if p and p["n"]:
                row.update({
                    "NumCandidatePalindromes": p["n"],
                    "MeanWholeIdentity": round(p["whole"] / p["n"], 2),
                    "MeanMiddleIdentity20bp": round(p["mid20"] / p["n"], 2),
                    "MeanRandomIdentity20bp": round(p["rand20"] / p["n"], 2),
                    "MeanMiddleMinusWhole": round(p["delta"] / p["n"], 2),
                })

        main_rows.append(row)

    all_cols = ["Order", "Species", "Haplotype", "Locus", "LatinName", "Source", "SourceDetail",
                "HaplotypeType", "Published", "NumContigs", "Contigs", "NumV", "NumV_productive", "NumV_with_RSS",
                "FracV_with_RSS", f"NumInversions_min{INV_MINLEN}bp",
                f"NumInversionsDiag_min{INV_MINLEN}bp", "LocusLength_bp", "InvCoverage_bp",
                "FracGenesOnInv", "NumD", "NumD_on_inv", "FracD_on_inv",
                "NumCandidatePalindromes", "MeanWholeIdentity", "MeanMiddleIdentity20bp",
                "MeanRandomIdentity20bp", "MeanMiddleMinusWhole"]

    write_table(main_rows, all_cols, "main_table")

    # --- species-level rollup ---
    by_species = defaultdict(list)
    for row in main_rows:
        by_species[(row["Order"], row["Species"])].append(row)

    def stat(vals):
        vals = [v for v in vals if v != ""]
        if not vals:
            return "", "", ""
        return (round(sum(vals) / len(vals), 3), min(vals), max(vals))

    species_rows = []
    for (order, species), rows in sorted(by_species.items()):
        haps = sorted({r["Haplotype"] for r in rows})
        srcs = sorted({r["Source"] for r in rows if r["Source"]})
        igh = [r for r in rows if r["Locus"] == "IGH"]
        igl = [r for r in rows if r["Locus"] == "IGL"]

        numv_igh_mean, numv_igh_min, numv_igh_max = stat([r["NumV"] for r in igh])
        numv_igl_mean, numv_igl_min, numv_igl_max = stat([r["NumV"] for r in igl])
        frac_inv_mean, frac_inv_min, frac_inv_max = stat([r.get("FracGenesOnInv", "") for r in igh])

        species_rows.append({
            "Order": order, "Species": species, "LatinName": latin.get(species, ""),
            "NumHaplotypes": len(haps), "Sources": ";".join(srcs),
            "AnyUnpublished": any(s.startswith("Unpublished") for s in srcs),
            "NumV_IGH_mean": numv_igh_mean, "NumV_IGH_min": numv_igh_min, "NumV_IGH_max": numv_igh_max,
            "NumV_IGL_mean": numv_igl_mean, "NumV_IGL_min": numv_igl_min, "NumV_IGL_max": numv_igl_max,
            "FracGenesOnInv_IGH_mean": frac_inv_mean, "FracGenesOnInv_IGH_min": frac_inv_min,
            "FracGenesOnInv_IGH_max": frac_inv_max,
        })

    species_cols = ["Order", "Species", "LatinName", "NumHaplotypes", "Sources", "AnyUnpublished",
                     "NumV_IGH_mean", "NumV_IGH_min", "NumV_IGH_max",
                     "NumV_IGL_mean", "NumV_IGL_min", "NumV_IGL_max",
                     "FracGenesOnInv_IGH_mean", "FracGenesOnInv_IGH_min", "FracGenesOnInv_IGH_max"]
    write_table(species_rows, species_cols, "species_summary", source_field="Sources")

    missing_latin = sorted({r["Species"] for r in main_rows if not r["LatinName"]})
    if missing_latin:
        print(f"NOTE: {len(missing_latin)} species have no LatinName:", ", ".join(missing_latin))

    indexed = load_indexed_haplotypes()
    build_v_gene_table(sources, latin, indexed)
    build_d_gene_table(sources, latin, indexed)
    write_readme()


def write_readme():
    readme = """# summary_tables

Generated by `annotation_tables/build_summary_tables.py` (run after
`assign_haplotype_source.py`) from the per-haplotype pipeline outputs in the
parent directory. Do not hand-edit -- regenerate instead.

Every table below comes in two files: the plain name (all haplotypes) and a
`_published` version that drops the unpublished Darwin's finches. The one
Red-winged Blackbird assembly that is also unpublished (Source =
"Unpublished / Pennell Lab assembly") is *not* dropped by `_published` --
only Darwin's finches are excluded, since that's the only exclusion asked
for so far. Check `Published`/`Source` directly if you need to filter that
one out too.

## main_table.csv

One row per Order/Species/Haplotype/Locus (IGH and IGL). A haplotype's locus
can span multiple contigs in `summary_features.csv`; those are summed here
(`NumContigs`, `Contigs` list the detail).

| Column | Meaning |
|---|---|
| LatinName | From IGH_VGP_table.tsv / species_traits_avonet.csv / NCBI organism name / manual lookup, in that priority. Blank if unresolved (see build script output). |
| Source, SourceDetail, Published | From `haplotype_sources.csv` -- VGP, CCGP, house finch/jay/seedeater pangenome, unpublished, or one of the other sequencing initiatives the VGP master sheet also tracks (Darwin Tree of Life, AmaZoomics, Sanger 25G, etc.), cross-referenced by accession against `/local/storage/kav67/VGP_details.csv`. SourceDetail is blank except for a few one-off assemblies outside all of the above. |
| HaplotypeType | From `haplotype_sources.csv` -- Maternal/Paternal (trio-phased), Hap1/Hap2 (Hi-C-phased, no parent-of-origin call), Primary/Alternate (solo pseudohaplotype, no parent-of-origin call), Single assembly (no haplotype pair released), or Merged. These are not interchangeable -- only Maternal/Paternal reflects an actual parent-of-origin assignment. Resolved from the real assembly filename (cross-checked against `/local/storage/kav67/VGP_details.csv` for a few VGP haplotypes saved locally under a generic name that doesn't preserve the original suffix). |
| NumV, NumV_productive, NumV_with_RSS, FracV_with_RSS | From `gene_list.csv`, filtered to `Passes Filtering == True`. "With RSS" = has a called heptamer and/or nonamer. |
| NumInversions_min250bp, NumInversionsDiag_min250bp | From `inversion_stats.tsv` at the 250 bp threshold (the only one emitted; matches the default used in `d_genes_on_inversions.py`). IGH only -- inversion detection in this pipeline is not run on IGL. |
| LocusLength_bp, InvCoverage_bp, FracGenesOnInv | Total self-alignment length, bp covered by inversions, and fraction of V genes falling in an inverted region (summed numerator/denominator across contigs, not averaged). |
| NumD, NumD_on_inv, FracD_on_inv | From `D_inversions.tsv`. IGH only. **Currently blank** -- that table was retired to `old_unused/` because the upstream D gene calls became a threshold-swept candidate list; see that folder's README. |
| NumCandidatePalindromes, MeanWholeIdentity, MeanMiddleIdentity20bp, MeanRandomIdentity20bp, MeanMiddleMinusWhole | From `palindromes.tsv` (`hairpin.py` output). `MeanMiddleMinusWhole` > 0 is the diagnostic for hairpin-like structure: the inversion's center is more self-similar than the whole alignment. IGH only. |

## species_summary.csv

One row per Order/Species, collapsing across haplotypes (mean/min/max --
several species have many individuals from pangenome projects, so a single
"the" value would hide real within-species variation).

`AnyUnpublished = TRUE` flags species with at least one unpublished haplotype
(currently: several Darwin's finches, plus one Red-winged Blackbird assembly)
-- check before including in any figure/table meant for public release.

## v_gene_table.csv

One row per called V gene from `gene_list.csv`, IGH + IGL only -- TRA/TRB/
TRG/TRD rows in that file are dropped, they're leftover annotation from a
shared upstream pipeline and out of scope for this project. Haplotype-level
context (`LatinName`, `Source`, `HaplotypeType`, `Published`) is joined in
from `haplotype_sources.csv` so the table stands alone. `PassesFiltering` is
included but *not* pre-filtered on -- this table has every called V gene
(real and filtered-out), unlike `NumV` in `main_table.csv` which counts only
`PassesFiltering == True` genes.

## d_gene_table.csv

One row per called D gene from `bird_d_genes.csv` (IGH only in this
pipeline). Kept as a separate file from `v_gene_table.csv` rather than one
combined gene table because V and D genes carry different RSS annotation:
V has a single downstream RSS (`Heptamer`/`Nonamer`/`RSS_DownstreamBp`), D
has RSS on both sides (`UpstreamHeptamer`/`UpstreamNonamer`/
`DownstreamHeptamer`/`DownstreamNonamer`/`LocationRelativeToVCluster`) --
forcing them into one table meant half the RSS columns were always blank
depending on gene type.

## Known gaps

- LatinName is unresolved for a handful of species with no NCBI accession and
  no AVONET match (see script stdout when re-run).
- No phylogenetic-regression-results table yet (e.g. phylolm fits of NumV ~
  n_inversions) -- those live in `plots/phylolm_tree.R` and friends and would
  need a dedicated export step in the `bird_ig` R environment.
- Ecological trait data (`species_traits_avonet.csv`) and pathogen exposure
  data are intentionally not joined in here -- keep genomic/provenance facts
  separate from ecological covariates used in specific analyses.
"""
    with open(os.path.join(OUT_DIR, "README.md"), "w") as f:
        f.write(readme)


if __name__ == "__main__":
    main()
