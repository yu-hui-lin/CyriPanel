# Sensitivity analyses from archived calls

`sensitivity.csv` uses panel size 5 and equal weighting of sample-level rates.
All evaluated observations, including no-calls, remain in denominators.
Primary scoring accepts any matching candidate. Single-candidate scoring requires
exactly one candidate and a match. Exact rational rates are in `per_sample_rates.csv`.
Aldy scores refer only to the retained major calls in the archived CSV; this audit
does not verify candidate selection against original Aldy solution files.

The *49 scenario forces every evaluation of HG00597, HG00609 and NA24631 to be
incorrect while preserving denominators. This is hypothetical rescoring, not a
filter-removal experiment or a demonstrated worst-case bound on filter effects.
Across all panel sizes and both aligners, no *49 call was observed in 15,960
archived evaluations of the 112 gold-non-*49 samples (no-calls included).

Across all 16,680 CyriPanel evaluations, candidate-inclusive discordance comprises:

- 106 no-calls.
- 430 phased calls with different per-haplotype gene-copy counts.
- 172 unphased calls with a different total gene-copy count.
- 270 phased calls with equal per-haplotype copy counts but different normal/hybrid counts.
- 8 unphased HG00597 calls whose total copy count does not resolve haplotype structure.

The structure signature is [phased, sorted [normal_count, hybrid_count] groups].
*5 contributes zero copies; *13, *36, *68 and *83 are collapsed to one hybrid class.
Other star identities are ignored. These coarse categories do not establish
biological equivalence, a causal mechanism, or validated structural truth.
The eight unresolved calls cannot support a claim that all discordance is structural.
This audit re-scores archived results; it does not rerun genotyping or CNV estimation.
In INPUT_SHA256SUMS, CSV paths are relative to --benchmark; script paths are
relative to the directory containing sensitivity_analyses.py.
