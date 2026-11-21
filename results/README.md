# Results

Output tables from the analysis scripts. The file name prefix says which script and step made each file:
- `single_*`: compare.ps1 (AlphaFold vs 5XJH)
- - `baseline_6EQE_vs_5XJH_*`: compare.ps1 (lab vs lab baseline)
  - - `ensemble_*`: ensemble.ps1 (26 wild-type chains)
    - - `family_*`: family.ps1 (11 relatives and the blind test)
      - - `part2_*`: part2.ps1 (ESMFold and FAST-PETase)
       
        - Running the scripts writes these into subfolders of results/ instead.
