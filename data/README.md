# Data files

## Tool resources

| File | Purpose |
|---|---|
| `CYP2D6_SNP_38.txt` | Paralogous CYP2D6/CYP2D7 discriminating sites (GRCh38) |
| `CYP2D6_target_variant_38.txt` | Target variants outside the homology blocks |
| `CYP2D6_target_variant_homology_region_38.txt` | Target variants within the homology blocks |
| `CYP2D6_haplotype_38.txt` | Sites used by the haplotype-resolved callers |
| `star_table.txt` | Star-allele definitions (BCyrius) |
| `PGxProbe_region_hg38.bed` | Panel target intervals; replace with your own panel |

## Reference genotypes

`CYP2D6_deCYPher.csv` — CYP2D6 diplotypes called by **deCYPher** from HPRC
haplotype-resolved assemblies, with the derived metabolizer phenotype and population
labels. 233 samples.

These are the reference genotypes against which CyriPanel was evaluated. They are
assembly-derived and therefore independent of the short-read evidence CyriPanel uses, so
they do not share its failure modes.

| Column | Notes |
|---|---|
| `Sample_ID` | Identifier as published by deCYPher |
| `Coriell_ID` | Coriell identifier; identical to `Sample_ID` except for the two GIAB samples, which deCYPher lists under their GIAB names (HG002 = NA24385, HG005 = NA24631). Join on this column |
| `Genotype` | Diplotype, unmodified from the source |
| `Phenotype` | Metabolizer phenotype |
| `Pop`, `Superpop` | 1000 Genomes phase 3 codes where applicable; `MKK` is HapMap 3 (Maasai in Kinyawa, Kenya); HG06807 is an HPRC sample not in either project and carries a descriptive label |

Changes made to the source file, none of which touch a genotype:

- Rows sorted by `Coriell_ID`; the `Coriell_ID` column added.
- `HG03492` had no population labels; filled as `PJL` / `SAS` from IGSR
  (https://www.internationalgenome.org/data-portal/sample/HG03492).
- `NA21309` carried the population as `MAASAI in KINYAWA, KENYA`, the only field in the
  file requiring quoting; recoded to `MKK`.

Reproduced with the authors' permission. Please cite:

> Chang T-Y, Liu Y-S, Lai H-S, Hung T-K, Lin H-F, Lin Y-H, Hsu C-L, Yang Y-C, Chen C-Y,
> Chen P-L, Hsu JS-J. deCYPher: star allele-resolution computational framework of
> pharmacogenes for haplotype-resolved long-read assemblies. *bioRxiv* 2025.10.13.681303.
> https://www.biorxiv.org/content/10.1101/2025.10.13.681303v2.full
