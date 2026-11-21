# How well does AlphaFold know a plastic-eating enzyme? Code and data

By Rayyan Niyazi, November 2026. Write-up: https://rayyanincoperations.github.io/petase-alphafold-study/

Everything runs on Windows PowerShell 5.1 with no extra installs. The structure maths is in C# and compiled on the fly with `Add-Type`.

## Reproduce from scratch
```powershell
.\fetch-data.ps1     # downloads every input (UniProt, AlphaFold DB, RCSB PDB) into data\
.\compare.ps1        # single comparison: AlphaFold vs 5XJH, and 6EQE vs 5XJH
.\ensemble.ps1       # all wild-type PETase chains: RMSD, sensitivity checks, RMSF, B-factors, correlations, Ser193, Trp185
.\family.ps1         # 11 relatives (EC 3.1.1.101): identity, TM-score, active site, blind test with model dates
.\part2.ps1          # ESMFold (via the ESM Atlas API) vs the crystals, and FAST-PETase (7SH6) real vs predicted change
```
Outputs go to `results\` (CSV) and `page\*.json` (data for the write-up page).

## Files
| File | What it does |
|---|---|
| `lib.cs` | PDB reading, Needleman–Wunsch alignment, quaternion (Horn) superposition, RMSD, iterative core fit, TM-score, dihedrals, Pearson and Spearman |
| `fetch-data.ps1` | Downloads all inputs, and records the RCSB search results and metadata as CSV |
| `compare.ps1` | The first, single-structure comparison |
| `ensemble.ps1` | The main analysis |
| `family.ps1` | Relatives and the blind test |
| `part2.ps1` | ESMFold comparison and the FAST-PETase variant |
| `referee-report.md` | Self-review that led to version 5 |

## Key choices
- **Wild-type** means zero amino-acid differences from UniProt A0A0K8P6T7 over residues 34–290. The RCSB "mutation count" field is *not* used, because it missed 8CRU, 8D1D and 6KY5.
- Residues are matched by sequence alignment, not by residue number, because structures use different numbering and carry tags.
- RMSDs use C-alpha atoms. The core is the residues present in every wild-type chain (34–289).
- Correlations are Spearman, with 2,000 bootstrap resamples (seed 42) for 95% intervals.
- "Strictly blind" means the AlphaFold model file's HEADER date is earlier than the lab structure's initial release date.

## Data versions
AlphaFold DB model v6 (AlphaFold Monomer v2.0 network). PDB and UniProt were as of 21 Nov 2026.
