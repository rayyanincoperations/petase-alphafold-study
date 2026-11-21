# Part 2: a second predictor (ESMFold) and an engineered variant (FAST-PETase).
#  A. ESMFold's prediction of wild-type PETase vs the same 26 wild-type crystal chains used in ensemble.ps1
#  B. Do predictors capture what FAST-PETase's five mutations (S121E/D186H/R224Q/N233K/R280A) do to the structure?
#     Real change: FAST-PETase crystal 7SH6 vs the wild-type crystal ensemble. Predicted change: ESMFold(FAST) vs ESMFold(WT).
# ESMFold is run through the public ESM Metagenomic Atlas API.
$Root = $PSScriptRoot
Add-Type -Path (Join-Path $Root "lib.cs")
$D = Join-Path $Root "data"; $Out = Join-Path $Root "results\part2"; New-Item -ItemType Directory -Force $Out | Out-Null
$uni = ((Get-Content (Join-Path $D "PETase_A0A0K8P6T7.fasta") | Select-Object -Skip 1) -join "").Trim()
$rng = New-Object Random 42
$MUT = @{ 121 = "E"; 186 = "H"; 224 = "Q"; 233 = "K"; 280 = "A" }

function Map-Chain($residues, [switch]$AllowMut) {
  $seq = [Lib]::Seq($residues); $map = [Lib]::Align($seq, $uni); $h = @{}; $diff = 0
  for ($i = 0; $i -lt $seq.Length; $i++) { $j = $map[$i]; if ($j -lt 0) { continue }
    if ($seq[$i] -eq $uni[$j] -or ($AllowMut -and $MUT.ContainsKey($j + 1))) { $h[$j + 1] = $residues[$i] } elseif ($j -ge 33) { $diff++ } }
  [pscustomobject]@{ map = $h; diff = $diff }
}
function Med($xs) { $s = @($xs | Sort-Object); $s[[int][math]::Floor(($s.Count - 1) / 2)] }
function Boot([double[]]$x, [double[]]$y) { $n = $x.Length; $v = New-Object 'System.Collections.Generic.List[double]'
  for ($b = 0; $b -lt 2000; $b++) { $xi = New-Object double[] $n; $yi = New-Object double[] $n
    for ($i = 0; $i -lt $n; $i++) { $k = $rng.Next($n); $xi[$i] = $x[$k]; $yi[$i] = $y[$k] }; $v.Add([Lib]::Spearman($xi, $yi)) }
  $s = $v.ToArray(); [Array]::Sort($s); [ordered]@{ rho = [math]::Round([Lib]::Spearman($x, $y), 2); lo = [math]::Round($s[49], 2); hi = [math]::Round($s[1949], 2) } }
function CA-List($map, $positions) { $l = New-Object 'System.Collections.Generic.List[double[]]'; foreach ($p in $positions) { $l.Add([Lib]::P($map[$p].CA)) }; ,$l }

# ---- ESMFold predictions (mature enzyme, residues 30-290)
$wtSeq = $uni.Substring(29)
$fastSeq = ($wtSeq.ToCharArray() | ForEach-Object -Begin { $i = 30 } -Process { if ($MUT.ContainsKey($i)) { $MUT[$i] } else { $_ }; $i++ }) -join ""
foreach ($v in @(@{ n = "WT"; s = $wtSeq }, @{ n = "FAST"; s = $fastSeq })) {
  $f = Join-Path $D "esmfold_$($v.n).pdb"
  if (-not (Test-Path $f)) { $r = Invoke-WebRequest -Method Post -Uri "https://api.esmatlas.com/foldSequence/v1/pdb/" -Body $v.s -UseBasicParsing -TimeoutSec 300
    [IO.File]::WriteAllBytes($f, $r.Content) }
}
$esmWT = (Map-Chain ([Lib]::Residues((Join-Path $D "esmfold_WT.pdb"), "A"))).map
$esmFAST = (Map-Chain ([Lib]::Residues((Join-Path $D "esmfold_FAST.pdb"), "A")) -AllowMut).map
$af = (Map-Chain ([Lib]::Residues((Join-Path $D "AF-A0A0K8P6T7_predicted.pdb"), "A"))).map
# ESMFold writes pLDDT on a 0-1 scale in some versions; put it on 0-100
$scale = if ((@($esmWT.Values | ForEach-Object { $_.CA.B }) | Measure-Object -Maximum).Maximum -le 1) { 100 } else { 1 }

# ---- Wild-type crystal chains (same rule as ensemble.ps1)
$chains = @()
foreach ($f in Get-ChildItem (Join-Path $D "wt") -Filter *.pdb) { foreach ($ch in [Lib]::Chains($f.FullName)) {
  $m = Map-Chain ([Lib]::Residues($f.FullName, $ch)); if ($m.map.Count -ge 200 -and $m.diff -eq 0) { $chains += [pscustomobject]@{ name = "$($f.BaseName)-$ch"; map = $m.map } } } }
$fastX = (Map-Chain ([Lib]::Residues((Join-Path $D "7SH6_fastpetase.pdb"), "A")) -AllowMut).map
$core = @(34..289 | Where-Object { $p = $_; $esmWT.ContainsKey($p) -and $esmFAST.ContainsKey($p) -and $fastX.ContainsKey($p) -and $af.ContainsKey($p) -and -not ($chains | Where-Object { -not $_.map.ContainsKey($p) }) })
"Wild-type chains: $($chains.Count); core residues: $($core.Count)"

# ---- A. ESMFold vs crystals, next to AlphaFold
$ca = @{}; foreach ($c in $chains) { $ca[$c.name] = CA-List $c.map $core }
$caE = CA-List $esmWT $core; $caA = CA-List $af $core
$esmR = @(foreach ($c in $chains) { [Lib]::Rmsd($caE, $ca[$c.name], [Lib]::Fit($caE, $ca[$c.name])) })
$afR  = @(foreach ($c in $chains) { [Lib]::Rmsd($caA, $ca[$c.name], [Lib]::Fit($caA, $ca[$c.name])) })
$labR = @(); for ($i = 0; $i -lt $chains.Count; $i++) { for ($j = $i + 1; $j -lt $chains.Count; $j++) {
  if ($chains[$i].name.Split("-")[0] -ne $chains[$j].name.Split("-")[0]) { $labR += [Lib]::Rmsd($ca[$chains[$i].name], $ca[$chains[$j].name], [Lib]::Fit($ca[$chains[$i].name], $ca[$chains[$j].name])) } } }
$esmVsAf = [Lib]::Rmsd($caE, $caA, [Lib]::Fit($caE, $caA))
"A. Median RMSD to wild-type crystals: ESMFold {0:N3} (range {1:N3}-{2:N3}), AlphaFold {3:N3}, lab-lab {4:N3}; ESMFold vs AlphaFold {5:N3}" -f (Med $esmR), ($esmR | Measure-Object -Minimum).Minimum, ($esmR | Measure-Object -Maximum).Maximum, (Med $afR), (Med $labR), $esmVsAf

$ens = @{}; Import-Csv (Join-Path $Root "results\ensemble\per_residue_ensemble.csv") | ForEach-Object { $ens[[int]$_.pos] = $_ }
$both = @($core | Where-Object { $ens.ContainsKey($_) })
$uE = [double[]]@($both | ForEach-Object { 100 - $esmWT[$_].CA.B * $scale }); $zb = [double[]]@($both | ForEach-Object { [double]$ens[$_].zb })
$rm = [double[]]@($both | ForEach-Object { [double]$ens[$_].rmsf })
$plE = @($both | ForEach-Object { $esmWT[$_].CA.B * $scale })
$corrE = [ordered]@{ esm_uncertainty_vs_bfactor = Boot $uE $zb; esm_uncertainty_vs_rmsf = Boot $uE $rm }
"   ESMFold mean pLDDT over core {0:N1} (range {1:N1}-{2:N1}); uncertainty vs B-factor rho {3} [{4},{5}]; vs RMSF rho {6} [{7},{8}]" -f ($plE | Measure-Object -Average).Average, ($plE | Measure-Object -Minimum).Minimum, ($plE | Measure-Object -Maximum).Maximum, $corrE.esm_uncertainty_vs_bfactor.rho, $corrE.esm_uncertainty_vs_bfactor.lo, $corrE.esm_uncertainty_vs_bfactor.hi, $corrE.esm_uncertainty_vs_rmsf.rho, $corrE.esm_uncertainty_vs_rmsf.lo, $corrE.esm_uncertainty_vs_rmsf.hi
$k193 = [array]::IndexOf($core, 193)
$e193 = @(foreach ($c in $chains) { $T = [Lib]::Fit($caE, $ca[$c.name]); [Lib]::Dist([Lib]::Apply($T, $caE[$k193]), $ca[$c.name][$k193]) })
"   Ser193: ESMFold median distance to crystals {0:N2} A" -f (Med $e193)

# ---- B. FAST-PETase: real vs predicted structural change
$caF = CA-List $fastX $core; $caEF = CA-List $esmFAST $core
$fastR = @(foreach ($c in $chains) { [Lib]::Rmsd($caF, $ca[$c.name], [Lib]::Fit($caF, $ca[$c.name])) })
$predGlobal = [Lib]::Rmsd($caEF, $caE, [Lib]::Fit($caEF, $caE))
"B. FAST-PETase crystal vs wild-type crystals: median {0:N3} (range {1:N3}-{2:N3}); lab-lab median {3:N3}" -f (Med $fastR), ($fastR | Measure-Object -Minimum).Minimum, ($fastR | Measure-Object -Maximum).Maximum, (Med $labR)
"   ESMFold FAST vs ESMFold WT: {0:N3} A" -f $predGlobal

# Per-residue: real change = 7SH6 distance from the wild-type average, in units of wild-type wobble (RMSF)
$mean = $ca[$chains[0].name]
for ($it = 0; $it -lt 3; $it++) { $acc = New-Object 'double[,]' $core.Count, 3
  foreach ($c in $chains) { $T = [Lib]::Fit($ca[$c.name], $mean); for ($k = 0; $k -lt $core.Count; $k++) { $p = [Lib]::Apply($T, $ca[$c.name][$k]); for ($q = 0; $q -lt 3; $q++) { $acc[$k, $q] += $p[$q] / $chains.Count } } }
  $mean = New-Object 'System.Collections.Generic.List[double[]]'; for ($k = 0; $k -lt $core.Count; $k++) { $mean.Add([double[]]@($acc[$k,0], $acc[$k,1], $acc[$k,2])) } }
$rmsf = New-Object double[] $core.Count
foreach ($c in $chains) { $T = [Lib]::Fit($ca[$c.name], $mean); for ($k = 0; $k -lt $core.Count; $k++) { $d = [Lib]::Dist([Lib]::Apply($T, $ca[$c.name][$k]), $mean[$k]); $rmsf[$k] += $d * $d / $chains.Count } }
$TF = [Lib]::Fit($caF, $mean); $TEF = [Lib]::Fit($caEF, $caE)
$rows = for ($k = 0; $k -lt $core.Count; $k++) { $pos = $core[$k]; $w = [math]::Sqrt($rmsf[$k])
  $real = [Lib]::Dist([Lib]::Apply($TF, $caF[$k]), $mean[$k]); $pred = [Lib]::Dist([Lib]::Apply($TEF, $caEF[$k]), $caE[$k])
  [pscustomobject]@{ pos = $pos; wt = $uni[$pos - 1]; mutated = $MUT.ContainsKey($pos); rmsf = [math]::Round($w, 3); real_change = [math]::Round($real, 3)
    real_z = [math]::Round(($real / [math]::Max($w, 0.05)), 2); predicted_change = [math]::Round($pred, 3) } }
$rows | Export-Csv (Join-Path $Out "fast_per_residue.csv") -NoTypeInformation -Encoding UTF8
"   Mutation sites (real = 7SH6 from wild-type average; predicted = ESMFold FAST vs WT), A:"
$rows | Where-Object mutated | Format-Table pos, wt, rmsf, real_change, real_z, predicted_change
"   Residues the crystal says moved most (real change > 3x wild-type wobble):"
$rows | Where-Object { $_.real_z -gt 3 } | Sort-Object real_change -Descending | Select-Object -First 10 | Format-Table pos, wt, mutated, rmsf, real_change, real_z, predicted_change
$corrF = Boot ([double[]]@($rows | ForEach-Object predicted_change)) ([double[]]@($rows | ForEach-Object real_change))
"   Predicted vs real per-residue change: rho {0} [{1},{2}]" -f $corrF.rho, $corrF.lo, $corrF.hi

[ordered]@{
  esm = [ordered]@{ median_to_crystals = [math]::Round((Med $esmR), 3); min = [math]::Round(($esmR | Measure-Object -Minimum).Minimum, 3); max = [math]::Round(($esmR | Measure-Object -Maximum).Maximum, 3)
    af_median = [math]::Round((Med $afR), 3); lab_median = [math]::Round((Med $labR), 3); esm_vs_af = [math]::Round($esmVsAf, 3)
    mean_plddt = [math]::Round(($plE | Measure-Object -Average).Average, 1); ser193 = [math]::Round((Med $e193), 2); correlations = $corrE }
  fast = [ordered]@{ crystal_vs_wt_median = [math]::Round((Med $fastR), 3); crystal_vs_wt_min = [math]::Round(($fastR | Measure-Object -Minimum).Minimum, 3); crystal_vs_wt_max = [math]::Round(($fastR | Measure-Object -Maximum).Maximum, 3)
    esm_fast_vs_wt = [math]::Round($predGlobal, 3); sites = @($rows | Where-Object mutated); pred_vs_real = $corrF }
  n_core = $core.Count; n_chains = $chains.Count
} | ConvertTo-Json -Depth 5 | Set-Content (Join-Path $Out "part2_summary.json") -Encoding UTF8
"Wrote results\part2\part2_summary.json"
