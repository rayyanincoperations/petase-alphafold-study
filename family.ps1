# Part 3: PETase's relatives. For each related PET-degrading enzyme (EC 3.1.1.101):
#  - sequence identity and 3D similarity to PETase
#  - which amino acids sit where PETase has its active site and Trp185
#  - blind test: AlphaFold DB model vs the enzyme's own lab structure
$Root = $PSScriptRoot
Add-Type -Path (Join-Path $Root "lib.cs")
$D = Join-Path $Root "data\family"; New-Item -ItemType Directory -Force $D | Out-Null
$Out = Join-Path $Root "results\family"; New-Item -ItemType Directory -Force $Out | Out-Null
$AF_TRAINING_CUTOFF = "2018-04-30"   # AlphaFold 2 was trained on PDB structures released before this date

$uniPet = ((Get-Content (Join-Path $Root "data\PETase_A0A0K8P6T7.fasta") | Select-Object -Skip 1) -join "").Trim()
$petRes = [Lib]::Residues((Join-Path $Root "data\6EQE_experimental.pdb"), "A")
$petSeq = [Lib]::Seq($petRes)
# UniProt position of each 6EQE residue
$pmap = [Lib]::Align($petSeq, $uniPet); $petPos = @{}; for ($i = 0; $i -lt $petSeq.Length; $i++) { if ($pmap[$i] -ge 0) { $petPos[$pmap[$i] + 1] = $i } }

$fam = Import-Csv (Join-Path $Root "data\pet_hydrolase_family.csv") | Where-Object { $_.uniprot -and $_.uniprot -ne "A0A0K8P6T7" -and $_.uniprot -notmatch "," }
$results = foreach ($g in ($fam | Group-Object uniprot)) {
  $uni = $g.Name
  $first = ($g.Group | Sort-Object released | Select-Object -First 1).released
  $best = $g.Group | Where-Object muts -eq "0" | Sort-Object { [double]$_.res } | Select-Object -First 1
  if (-not $best) { continue }
  $pdb = $best.entity.Split("_")[0]; $chain = $best.chains.Substring(0, 1)
  $pf = Join-Path $D "$pdb.pdb"
  if (-not (Test-Path $pf)) { try { Invoke-WebRequest "https://files.rcsb.org/download/$pdb.pdb" -OutFile $pf -UseBasicParsing -ErrorAction Stop } catch { "No PDB-format file for $pdb, skipped"; continue } }
  $lab = [Lib]::Residues($pf, $chain); $labSeq = [Lib]::Seq($lab)

  # PETase (6EQE) vs this relative's lab structure
  $m = [Lib]::Align($petSeq, $labSeq)
  $ident = [Lib]::Identity($petSeq, $labSeq, $m)
  $mov = New-Object 'System.Collections.Generic.List[double[]]'; $tgt = New-Object 'System.Collections.Generic.List[double[]]'
  for ($i = 0; $i -lt $petSeq.Length; $i++) { if ($m[$i] -ge 0) { $mov.Add([Lib]::P($petRes[$i].CA)); $tgt.Add([Lib]::P($lab[$m[$i]].CA)) } }
  $kept = 0; $coreR = 0.0; [void][Lib]::FitCore($mov, $tgt, 2.0, [ref]$kept, [ref]$coreR)
  $tm = [math]::Round([Lib]::TmScore($mov, $tgt, $petRes.Count), 3)
  $equiv = foreach ($pos in 160, 185, 206, 237) { $i = $petPos[$pos]; if ($i -ne $null -and $m[$i] -ge 0) { [Lib]::One($lab[$m[$i]].Res) + $lab[$m[$i]].Num } else { "-" } }

  # Blind test: AlphaFold DB model for this UniProt vs its own lab structure
  $afR = $null; $afPl = $null; $afVer = $null; $afN = $null; $afCore = $null; $afK = $null; $afOffN = $null; $afOffPl = $null; $afDate = $null
  try {
    $api = Invoke-RestMethod "https://alphafold.ebi.ac.uk/api/prediction/$uni" -ErrorAction Stop
    $af = Join-Path $D "AF-$uni.pdb"; if (-not (Test-Path $af)) { Invoke-WebRequest $api[0].pdbUrl -OutFile $af -UseBasicParsing }
    $afVer = $api[0].latestVersion; $afPl = $api[0].globalMetricValue
    # Date in the model file's HEADER line (when the file was made; the prediction can't be later)
    $hdr = (Get-Content $af -TotalCount 1).Substring(50).Trim()
    $afDate = [datetime]::ParseExact($hdr, "dd-MMM-yy", [Globalization.CultureInfo]::InvariantCulture).ToString("yyyy-MM-dd")
    $afRes = [Lib]::Residues($af, "A"); $afSeq = [Lib]::Seq($afRes); $am = [Lib]::Align($labSeq, $afSeq)
    $mv = New-Object 'System.Collections.Generic.List[double[]]'; $tg = New-Object 'System.Collections.Generic.List[double[]]'
    for ($i = 0; $i -lt $labSeq.Length; $i++) { $j = $am[$i]; if ($j -ge 0 -and $labSeq[$i] -eq $afSeq[$j]) { $mv.Add([Lib]::P($afRes[$j].CA)); $tg.Add([Lib]::P($lab[$i].CA)) } }
    $T = [Lib]::Fit($mv, $tg); $afR = [math]::Round([Lib]::Rmsd($mv, $tg, $T), 3); $afN = $mv.Count
    # Core after dropping residues >2 A off, and how confident AlphaFold was about the ones it got wrong
    $afK = 0; $afC = 0.0; [void][Lib]::FitCore($mv, $tg, 2.0, [ref]$afK, [ref]$afC); $afCore = [math]::Round($afC, 3)
    $plOff = @(); for ($i = 0; $i -lt $labSeq.Length; $i++) { $j = $am[$i]; if ($j -ge 0 -and $labSeq[$i] -eq $afSeq[$j]) {
      if ([Lib]::Dist([Lib]::Apply($T, [Lib]::P($afRes[$j].CA)), [Lib]::P($lab[$i].CA)) -gt 2) { $plOff += $afRes[$j].CA.B } } }
    $afOffN = $plOff.Count; $afOffPl = if ($plOff.Count) { [math]::Round(($plOff | Measure-Object -Average).Average, 1) } else { $null }
  } catch { "No AlphaFold DB model for $uni" }

  [pscustomobject]@{
    uniprot = $uni; organism = ($g.Group[0].org); name = ($g.Group[0].desc); pdb = "$pdb-$chain"; resolution = [double]$best.res
    first_structure = $first; blind = ($first -gt $AF_TRAINING_CUTOFF)
    seq_identity_pct = [math]::Round($ident, 1); core_rmsd_to_petase = [math]::Round($coreR, 2); core_residues = $kept; length = $lab.Count; tm_score = $tm
    at_S160 = $equiv[0]; at_W185 = $equiv[1]; at_D206 = $equiv[2]; at_H237 = $equiv[3]
    af_vs_lab_rmsd = $afR; af_residues = $afN; af_mean_plddt = $afPl; af_version = $afVer
    af_core_rmsd = $afCore; af_core_residues = $afK; af_off_by_2A = $afOffN; af_off_mean_plddt = $afOffPl
    af_core_pct = $(if ($afN) { [math]::Round(100 * $afK / $afN, 1) } else { $null })
    af_model_date = $afDate; strict_blind = ($afDate -and $afDate -lt $first)
  }
}
$results = @($results | Where-Object { $_ -is [pscustomobject] })
$results | Export-Csv (Join-Path $Out "family.csv") -NoTypeInformation -Encoding UTF8
$results | ConvertTo-Json -Depth 4 -Compress | Set-Content (Join-Path $Root "page\family.json") -Encoding UTF8
$results | Sort-Object first_structure | Format-Table uniprot, organism, pdb, first_structure, blind, seq_identity_pct, core_rmsd_to_petase, core_residues, at_S160, at_W185, at_D206, at_H237, tm_score, af_vs_lab_rmsd, af_core_rmsd, af_core_pct, af_model_date, strict_blind -AutoSize | Out-String -Width 260

