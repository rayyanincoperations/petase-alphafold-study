# Referee report: "How well does AlphaFold know a plastic-eating enzyme?" (v4)

Reviewed 21 Nov 2026, written as a journal referee would. Every major point was checked against the project's own data.

## Summary of the submission
The study compares the AlphaFold DB model of *Ideonella sakaiensis* PETase with an ensemble of crystal structures. It also correlates pLDDT with crystallographic B-factors and with variation across the ensemble, looks at the Trp185 (W156) side-chain conformation, compares 11 PET hydrolase relatives, and tests AlphaFold on 6 relatives solved after its training cutoff.

## Overall assessment
The question is clear and well chosen, the analysis is reproducible in principle, and the limitations section is unusually honest for work at this level. The main conclusions are probably right. But one error in how the dataset was built (M1) affects every ensemble number, and some statements claim more than the data shows (M2). **Recommendation: major revision.** M1 and M2 are required, M3–M6 would strengthen the paper, and the minor points are optional.

---

## Major issues

**M1. The "wild-type" set includes engineered and non-native proteins.** *(Required)*
Wild-type entries were selected with the RCSB `rcsb_mutation_count = 0` field. That field is unreliable. Checking each sequence against UniProt over residues 30–290 shows:
- **8CRU** is an ancestral sequence reconstruction ("PETase ancestral sequence reconstruction 008"), with ≥29 substitutions and 88% identity.
- **8D1D** is a computationally stabilised PROSS design, with ≥17 substitutions and 93% identity.
- **6KY5** is titled "hydrolase mutant" and has 10 substitutions (96% identity).
These three entries give 6 of the 32 chains. They are different proteins, and 8CRU was the chain furthest from all the others (average 0.38 Å).
*Fix:* define wild-type as zero substitutions against UniProt in the mature region (expression tags allowed). That leaves 13 entries. A preliminary re-run on that set gives lab vs lab median 0.253 Å (n = 279, range 0.06–0.36) and AlphaFold vs lab median 0.266 Å (n = 25, range 0.16–0.35). **The conclusion survives**, but every ensemble figure, the rank, the RMSF, the correlations and the Trp185 plot need recomputing.

**M2. Some statements claim more than the data shows.** *(Required)*
"As close to the experiments as they are to each other" and "as good as a second experiment" are slightly stronger than the evidence. AlphaFold's median (0.265 Å) sits a little above the lab vs lab median, and it ranks 23rd of 33 by average agreement. The accurate claim is: *AlphaFold falls within the range of experimental variation, slightly above its median.*
The lab vs lab distribution also mixes two kinds of pair. Pairs from the same crystal (n = 25, median 0.18 Å) are much closer than pairs from different crystals (median 0.26 Å). Same-crystal pairs should be shown separately or left out, because the fair comparison for AlphaFold is different crystals.

**M3. The ensemble isn't independent.** *(Recommended)*
13 of the 32 chains come from five entries by one group under similar conditions (6EQD–6EQH). The lab vs lab distribution and the RMSF are weighted towards that group's crystal form.
*Fix:* add a sensitivity analysis using one chain per entry, or weight pairs by entry, and report whether the medians and ρ change.

**M4. The "core" definition favours ordered regions.** *(Recommended)*
Requiring a residue to be present in every chain removes 43 of 257 residues (17%). Those are mostly loops that are disordered in some crystals, which are exactly the regions where prediction and flexibility are most interesting. Leaving them out probably lowers every RMSD and may change the correlations.
*Fix:* for each pair, superpose on the residues those two structures share, and report how much the results move.

**M5. The blind test needs tighter framing.** *(Recommended)*
(a) AlphaFold DB models are built with templates, so the test is blind only to the protein's own structure. The page says this, which is good. It can go further, because the model files carry dates:

| Enzyme (lab structure, released) | Model file date | Built before its structure existed? |
|---|---|---|
| 7QJO (Dec 2022) | Jun 2022 | **Yes** |
| 8GZD (Jul 2023) | Jun 2022 | **Yes** |
| 9O9W (Oct 2025) | Jun 2022 | **Yes** |
| 7QJT (Dec 2022) | Mar 2025 | No |
| 8BRA (Apr 2023) | Mar 2025 | No |
| 7QJS (Dec 2022) | Aug 2025 | No |

So only **three** predictions are strictly blind: those models existed before their own lab structures did. Their core RMSDs are 0.35, 0.34 and 0.42 Å, so the conclusion holds even on this strict subset. Split the blind test into these two tiers. One caution: a header date may be when the file was regenerated rather than when the prediction was run, so say "dated" rather than "predicted". (The PETase model is dated Aug 2025, so it isn't blind in any sense.)
(b) "Core within 2 Å" is defined by the comparison itself (iterative pruning), which flatters the result. Report the fraction of residues kept (95–100% here), so readers can see it's mild.
(c) The sample is small: n = 6, one structure each. Avoid statements like "6/6", which suggest a benchmark.

**M6. The family alignment method is basic.** *(Recommended)*
Identity values come from a simple Needleman–Wunsch alignment that scores only identical/non-identical residues (+5/−3/−6, free end gaps), not a substitution matrix (e.g. BLOSUM62) or a structure-based alignment. At 44–56% identity this is probably fine, but identity values depend on the method. Reporting a structure-based score (e.g. TM-score) would make the "shape is more conserved than sequence" claim more robust.

## Minor issues

- **m1. Ser193.** The 0.79 Å error is isolated: neighbours are at 0.16 and 0.34 Å, and the crystals agree closely there (RMSF 0.19 Å). Before interpreting it, check the backbone angles (φ/ψ) at 192–194 for a peptide-plane flip, and confirm that AlphaFold is off from *every* chain and not just from the average. As it stands, the text slightly over-reads one residue.
- **m2. 8BRA alanine.** Both the 8BRA structure *and* its UniProt record (A0AA82WPD4) have Ala at the catalytic-serine position. This is probably an inactive variant used to capture the bound product (the entry is titled "PHL7 bound to terephthalic acid"), but the page shouldn't imply it either way. Check the original 8BRA paper and say so explicitly.
- **m3. Version labelling.** The model file's header says "AlphaFold Monomer v2.0", while the page says "model v6". These are the database release and the network version. State both, plus the download date.
- **m4. No CI on the AI-error correlations.** ρ = 0.51 (AI error vs RMSF) and 0.48 have no confidence intervals, unlike the others. Add them.
- **m5. Signal peptide.** Low pLDDT is *consistent with* the peptide being unstructured, but it's also simply absent from all the training structures. Say "consistent with" rather than "matches the biology".
- **m6. B-factors.** Refinement methods (e.g. TLS) vary between entries and affect B-factors. The per-chain z-scoring helps but doesn't remove this. One sentence would cover it.
- **m7. Figure 3** should colour same-crystal pairs differently (see M2).
- **m8. Reproducibility.** The code exists only on a local disk. Publishing `lib.cs`, `ensemble.ps1`, `family.ps1` and a download script on GitHub, with the list of PDB IDs and retrieval dates, would make the study fully reproducible. That's a strong signal to any lab reading it.

## Presentation
The writing is clear and well pitched for a general reader, and the figures are informative and honest (the Ser193 exception is shown rather than hidden). Defining RMSD, pLDDT and B-factors up front works well. The abstract will need updating after M1/M2.

## What the author did well
- Tested against an ensemble rather than a single structure, which is the right control and one many published comparisons skip.
- Used both B-factors and variation between structures as independent flexibility measures.
- Reported negative and awkward results (Ser193, AlphaFold ranking 23rd of 33).
- Stated the training-set and template caveats openly.
