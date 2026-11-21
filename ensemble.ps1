# Part 2 (revised after referee review): AlphaFold vs every truly wild-type PETase crystal structure.
# Wild-type = zero amino-acid differences from UniProt A0A0K8P6T7 over residues 34-290 (the part every crystal shows);
# residues 30-33 are skipped because several constructs carry expression-tag residues there. The RCSB "mutation count" field is NOT used (it missed
# engineered entries such as 8CRU, 8D1D and 6KY5).
# Outputs: results/ensemble/*.csv and page/ensemble.json
$Root = $PSScriptRoot
Add-Type -Path (Join-Path $Root "lib.cs")
$Out = Join-Path $Root "results\ensemble"; New-Item -ItemType Directory -Force $Out | Out-Null
$meta = @{}; Import-Csv (Join-Path $Root "data\petase_structures.csv") | ForEach-Object { $meta[$_.id] = $_ }
$uni = ((Get-Content (Join-Path $Root "data\PETase_A0A0K8P6T7.fasta") | Select-Object -Skip 1) -join "").Trim()
$rng = New-Object Random 42

function Map-Chain($residues) {
  $seq = [Lib]::Seq($residues); $map = [Lib]::Align($seq, $uni); $h = @{}; $diff = @()
  for ($i = 0; $i -lt $seq.Length; $i++) { $j = $map[$i]; if ($j -lt 0) { continue }
    if ($seq[$i] -eq $uni[$j]) { $h[$j + 1] = $residues[$i] } elseif ($j -ge 33) { $diff += "$($uni[$j])$($j + 1)$($seq[$i])" } }
  [pscustomobject]@{ map = $h; diff = $diff }
}
function Stats($xs) { $s = @($xs | Sort-Object); if (-not $s.Count) { return $null }
  [pscustomobject]@{ n = $s.Count; min = $s[0]; median = $s[[int][math]::Floor(($s.Count - 1) / 2)]; mean = [math]::Round(($s | Measure-Object -Average).Average, 3); max = $s[-1] } }
function Boot([double[]]$x, [double[]]$y) {
  $n = $x.Length; $v = New-Object 'System.Collections.Generic.List[double]'
  for ($b = 0; $b -lt 2000; $b++) { $xi = New-Object double[] $n; $yi = New-Object double[] $n
    for ($i = 0; $i -lt $n; $i++) { $k = $rng.Next($n); $xi[$i] = $x[$k]; $yi[$i] = $y[$k] }; $v.Add([Lib]::Spearman($xi, $yi)) }
  $s = $v.ToArray(); [Array]::Sort($s)
  [ordered]@{ rho = [math]::Round([Lib]::Spearman($x, $y), 2); lo = [math]::Round($s[49], 2); hi = [math]::Round($s[1949], 2) }
}
function Lst() { ,(New-Object 'System.Collections.Generic.List[double[]]') }

# ---- Load structures, keep only true wild-type chains
$af = (Map-Chain ([Lib]::Residues((Join-Path $Root "data\AF-A0A0K8P6T7_predicted.pdb"), "A"))).map
$entries = New-Object System.Collections.ArrayList; $excluded = @()
foreach ($f in Get-ChildItem (Join-Path $Root "data\wt") -Filter *.pdb | Sort-Object Name) {
  $id = $f.BaseName
  foreach ($ch in [Lib]::Chains($f.FullName)) {
    $m = Map-Chain ([Lib]::Residues($f.FullName, $ch))
    if ($m.map.Count -lt 200) { continue }
    if ($m.diff.Count -gt 0) { $excluded += [pscustomobject]@{ chain = "$id-$ch"; differences = $m.diff.Count; examples = (($m.diff | Select-Object -First 5) -join " ") }; continue }
    [void]$entries.Add([pscustomobject]@{ name = "$id-$ch"; id = $id; chain = $ch; map = $m.map; res = [double]$meta[$id].res; released = $meta[$id].released })
  }
}
$nEntries = @($entries | Select-Object -ExpandProperty id -Unique).Count
"Wild-type chains: $($entries.Count) from $nEntries entries"
"Excluded (not wild-type):"; $excluded | Group-Object { $_.chain.Split("-")[0] } | ForEach-Object { "  $($_.Name): $($_.Group[0].differences) differences, e.g. $($_.Group[0].examples)" }

$core = @(34..290 | Where-Object { $pos = $_; $af.ContainsKey($pos) -and -not ($entries | Where-Object { -not $_.map.ContainsKey($pos) }) })
"Core residues: $($core.Count) of 257 ($($core[0])-$($core[-1]))"

function CA-List($map, $positions) { $l = New-Object 'System.Collections.Generic.List[double[]]'; foreach ($p in $positions) { $l.Add([Lib]::P($map[$p].CA)) }; ,$l }
$all = @([pscustomobject]@{ name = "AlphaFold"; id = "AF"; map = $af }) + @($entries)
$ca = @{}; foreach ($e in $all) { $ca[$e.name] = CA-List $e.map $core }

# ---- 1. Pairwise RMSD on the shared core, plus a sensitivity version on each pair's own shared residues
$pairs = New-Object System.Collections.ArrayList
for ($i = 0; $i -lt $all.Count; $i++) { for ($j = $i + 1; $j -lt $all.Count; $j++) {
  $a = $all[$i]; $b = $all[$j]
  $T = [Lib]::Fit($ca[$a.name], $ca[$b.name]); $r = [Lib]::Rmsd($ca[$a.name], $ca[$b.name], $T)
  $own = @(34..290 | Where-Object { $a.map.ContainsKey($_) -and $b.map.ContainsKey($_) })
  $ma = CA-List $a.map $own; $mb = CA-List $b.map $own; $T2 = [Lib]::Fit($ma, $mb); $r2 = [Lib]::Rmsd($ma, $mb, $T2)
  $kind = if ($a.id -eq "AF") { "ai-lab" } elseif ($a.id -eq $b.id) { "same-crystal" } else { "lab-lab" }
  [void]$pairs.Add([pscustomobject]@{ a = $a.name; b = $b.name; kind = $kind; rmsd = [math]::Round($r, 3); rmsd_own = [math]::Round($r2, 3); n_own = $own.Count }) } }
$pairs | Export-Csv (Join-Path $Out "pairwise_rmsd.csv") -NoTypeInformation -Encoding UTF8

$labLab = @($pairs | Where-Object kind -eq "lab-lab" | ForEach-Object rmsd)
$same   = @($pairs | Where-Object kind -eq "same-crystal" | ForEach-Object rmsd)
$aiLab  = @($pairs | Where-Object kind -eq "ai-lab" | ForEach-Object rmsd)
"Different crystals:"; Stats $labLab | Format-List; "Same crystal:"; Stats $same | Format-List; "AI vs lab:"; Stats $aiLab | Format-List

# Where AlphaFold sits: share of different-crystal lab pairs that agree better than AlphaFold's median
$pctBelow = [math]::Round(100 * @($labLab | Where-Object { $_ -lt (Stats $aiLab).median }).Count / $labLab.Count)
"Different-crystal pairs closer than AlphaFold's median: $pctBelow%"

# Each structure's average RMSD to lab chains from OTHER crystals
$avgToLab = foreach ($e in $all) {
  if ($e.name -eq "AlphaFold") { $v = $aiLab }
  else { $v = @($pairs | Where-Object { ($_.a -eq $e.name -or $_.b -eq $e.name) -and $_.kind -eq "lab-lab" } | ForEach-Object rmsd) }
  [pscustomobject]@{ name = $e.name; mean_rmsd_to_lab = [math]::Round(($v | Measure-Object -Average).Average, 3) } }
$avgToLab = @($avgToLab | Sort-Object mean_rmsd_to_lab)
$afRank = [array]::IndexOf(@($avgToLab | ForEach-Object name), "AlphaFold") + 1
"AlphaFold rank: $afRank of $($avgToLab.Count)"

# Sensitivity A: one chain per crystal (first chain of each entry)
$firsts = @($entries | Group-Object id | ForEach-Object { $_.Group[0].name })
$s1Lab = @($pairs | Where-Object { $_.kind -eq "lab-lab" -and $firsts -contains $_.a -and $firsts -contains $_.b } | ForEach-Object rmsd)
$s1AI  = @($pairs | Where-Object { $_.kind -eq "ai-lab" -and $firsts -contains $_.b } | ForEach-Object rmsd)
# Sensitivity B: one chain per research group's crystal form; entries 6EQD-6EQH (one study) count once
$grp = @($firsts | Where-Object { $_ -notmatch "^6EQ[DFGH]" })
$s2Lab = @($pairs | Where-Object { $_.kind -eq "lab-lab" -and $grp -contains $_.a -and $grp -contains $_.b } | ForEach-Object rmsd)
$s2AI  = @($pairs | Where-Object { $_.kind -eq "ai-lab" -and $grp -contains $_.b } | ForEach-Object rmsd)
# Sensitivity C: each pair superposed on all residues both structures have (not just the shared core)
$s3Lab = @($pairs | Where-Object kind -eq "lab-lab" | ForEach-Object rmsd_own)
$s3AI  = @($pairs | Where-Object kind -eq "ai-lab" | ForEach-Object rmsd_own)
$sens = @(
  [ordered]@{ name = "Main analysis (shared core, all chains)"; lab = Stats $labLab; ai = Stats $aiLab }
  [ordered]@{ name = "One chain per crystal"; lab = Stats $s1Lab; ai = Stats $s1AI }
  [ordered]@{ name = "One crystal per study (6EQD-6EQH counted once)"; lab = Stats $s2Lab; ai = Stats $s2AI }
  [ordered]@{ name = "Each pair on all residues it shares"; lab = Stats $s3Lab; ai = Stats $s3AI }
)
"Sensitivity:"; $sens | ForEach-Object { "  {0,-48} lab median {1:N3} (n={2})  AI median {3:N3} (n={4})" -f $_.name, $_.lab.median, $_.lab.n, $_.ai.median, $_.ai.n }

# ---- 2. Mean lab structure, per-residue wobble (RMSF), B-factors
$mean = $ca[($entries | Sort-Object res | Select-Object -First 1).name]
for ($iter = 0; $iter -lt 3; $iter++) {
  $acc = New-Object 'double[,]' $core.Count, 3
  foreach ($e in $entries) { $T = [Lib]::Fit($ca[$e.name], $mean); for ($k = 0; $k -lt $core.Count; $k++) { $p = [Lib]::Apply($T, $ca[$e.name][$k]); for ($c = 0; $c -lt 3; $c++) { $acc[$k, $c] += $p[$c] / $entries.Count } } }
  $mean = New-Object 'System.Collections.Generic.List[double[]]'; for ($k = 0; $k -lt $core.Count; $k++) { $mean.Add([double[]]@($acc[$k,0], $acc[$k,1], $acc[$k,2])) }
}
$sq = New-Object 'double[]' $core.Count
foreach ($e in $entries) { $T = [Lib]::Fit($ca[$e.name], $mean)
  for ($k = 0; $k -lt $core.Count; $k++) { $d = [Lib]::Dist([Lib]::Apply($T, $ca[$e.name][$k]), $mean[$k]); $sq[$k] += $d * $d / $entries.Count } }
$Taf = [Lib]::Fit($ca["AlphaFold"], $mean)
$toMean = @(foreach ($e in $all) { $T = [Lib]::Fit($ca[$e.name], $mean); [pscustomobject]@{ name = $e.name; rmsd_to_mean = [math]::Round([Lib]::Rmsd($ca[$e.name], $mean, $T), 3) } }) | Sort-Object rmsd_to_mean
$afMeanRank = [array]::IndexOf(@($toMean | ForEach-Object name), "AlphaFold") + 1

$zb = New-Object 'double[]' $core.Count
foreach ($e in $entries) { $bs = @($core | ForEach-Object { $e.map[$_].CA.B }); $mu = ($bs | Measure-Object -Average).Average
  $sd = [math]::Sqrt((($bs | ForEach-Object { ($_ - $mu) * ($_ - $mu) }) | Measure-Object -Sum).Sum / $bs.Count)
  for ($k = 0; $k -lt $core.Count; $k++) { $zb[$k] += (($bs[$k] - $mu) / $sd) / $entries.Count } }

$rows = for ($k = 0; $k -lt $core.Count; $k++) { $pos = $core[$k]
  [pscustomobject]@{ pos = $pos; aa = $af[$pos].Res; plddt = $af[$pos].CA.B; rmsf = [math]::Round([math]::Sqrt($sq[$k]), 3)
    zb = [math]::Round($zb[$k], 3); af_dev = [math]::Round([Lib]::Dist([Lib]::Apply($Taf, $ca["AlphaFold"][$k]), $mean[$k]), 3) } }
$rows | Export-Csv (Join-Path $Out "per_residue_ensemble.csv") -NoTypeInformation -Encoding UTF8

$unc = [double[]]@($rows | ForEach-Object { 100 - $_.plddt }); $rmsf = [double[]]@($rows | ForEach-Object rmsf)
$zbA = [double[]]@($rows | ForEach-Object zb); $dev = [double[]]@($rows | ForEach-Object af_dev)
$corr = [ordered]@{
  uncertainty_vs_bfactor = Boot $unc $zbA; uncertainty_vs_rmsf = Boot $unc $rmsf; rmsf_vs_bfactor = Boot $rmsf $zbA
  afdev_vs_rmsf = Boot $dev $rmsf; afdev_vs_uncertainty = Boot $dev $unc }
"Spearman (95% CI):"; $corr.GetEnumerator() | ForEach-Object { "  {0,-24} {1} [{2}, {3}]" -f $_.Key, $_.Value.rho, $_.Value.lo, $_.Value.hi }
"pLDDT range in core: $(($rows | Measure-Object plddt -Minimum -Maximum | ForEach-Object { "$($_.Minimum)-$($_.Maximum)" }))"
"Top AI errors:"; $rows | Sort-Object af_dev -Descending | Select-Object -First 5 | Format-Table pos, aa, plddt, rmsf, af_dev

# ---- 3. Ser193 check: is AlphaFold off from EVERY chain there, and is the backbone flipped?
function Dh($r1, $n1, $r2, $n2, $r3, $n3, $r4, $n4) { [math]::Round([Lib]::Dihedral([Lib]::P($r1.Atoms[$n1]), [Lib]::P($r2.Atoms[$n2]), [Lib]::P($r3.Atoms[$n3]), [Lib]::P($r4.Atoms[$n4])), 0) }
function Angles($m) { [ordered]@{ psi192 = Dh $m[192] "N" $m[192] "CA" $m[192] "C" $m[193] "N"; phi193 = Dh $m[192] "C" $m[193] "N" $m[193] "CA" $m[193] "C"
  psi193 = Dh $m[193] "N" $m[193] "CA" $m[193] "C" $m[194] "N"; phi194 = Dh $m[193] "C" $m[194] "N" $m[194] "CA" $m[194] "C" } }
$k193 = [array]::IndexOf($core, 193)
$d193 = @(foreach ($e in $entries) { $T = [Lib]::Fit($ca["AlphaFold"], $ca[$e.name]); [Lib]::Dist([Lib]::Apply($T, $ca["AlphaFold"][$k193]), $ca[$e.name][$k193]) })
$labAng = @($entries | ForEach-Object { Angles $_.map })
$ser = [ordered]@{
  af_to_each_chain = [ordered]@{ min = [math]::Round(($d193 | Measure-Object -Minimum).Minimum, 2); max = [math]::Round(($d193 | Measure-Object -Maximum).Maximum, 2) }
  af_angles = Angles $af
  lab_angles_median = [ordered]@{}
}
foreach ($key in "psi192", "phi193", "psi193", "phi194") { $v = @($labAng | ForEach-Object { $_[$key] } | Sort-Object); $ser.lab_angles_median[$key] = $v[[int][math]::Floor(($v.Count - 1) / 2)]
  $ser["lab_" + $key + "_range"] = "$($v[0]) to $($v[-1])" }
"Ser193:"; $ser | ConvertTo-Json -Depth 4

# ---- 4. Trp185 side-chain angles and the active site
function Chi($r) { if (-not ($r -and $r.Atoms.ContainsKey("CB") -and $r.Atoms.ContainsKey("CG") -and $r.Atoms.ContainsKey("CD1"))) { return $null }
  [pscustomobject]@{ chi1 = [math]::Round([Lib]::Dihedral([Lib]::P($r.Atoms["N"]), [Lib]::P($r.CA), [Lib]::P($r.Atoms["CB"]), [Lib]::P($r.Atoms["CG"])), 1)
                     chi2 = [math]::Round([Lib]::Dihedral([Lib]::P($r.CA), [Lib]::P($r.Atoms["CB"]), [Lib]::P($r.Atoms["CG"]), [Lib]::P($r.Atoms["CD1"])), 1) } }
$w185 = @(foreach ($e in @([pscustomobject]@{ name = "AlphaFold"; map = $af; released = "" }) + @($entries)) {
  $c = Chi $e.map[185]; if ($c) { [pscustomobject]@{ name = $e.name; released = $e.released; chi1 = $c.chi1; chi2 = $c.chi2 } } })
$w185 | Export-Csv (Join-Path $Out "trp185.csv") -NoTypeInformation -Encoding UTF8
$site = foreach ($e in $entries) { $T = [Lib]::Fit($ca["AlphaFold"], $ca[$e.name]); $s = 0.0; $n = 0
  foreach ($pos in 160, 206, 237) { foreach ($kv in $af[$pos].Atoms.GetEnumerator()) { if ($e.map[$pos].Atoms.ContainsKey($kv.Key)) {
    $d = [Lib]::Dist([Lib]::Apply($T, [Lib]::P($kv.Value)), [Lib]::P($e.map[$pos].Atoms[$kv.Key])); $s += $d * $d; $n++ } } }
  [math]::Round([math]::Sqrt($s / $n), 3) }
"Active site:"; Stats $site | Format-List

# ---- JSON for the page
$json = [ordered]@{
  n_chains = $entries.Count; n_entries = $nEntries; core = "$($core[0])-$($core[-1])"; n_core = $core.Count
  excluded = @($excluded | Group-Object { $_.chain.Split("-")[0] } | ForEach-Object { [ordered]@{ entry = $_.Name; chains = $_.Count; differences = $_.Group[0].differences } })
  chains = @($entries | ForEach-Object { [ordered]@{ name = $_.name; res = $_.res; released = $_.released } })
  lab_lab = Stats $labLab; same_crystal = Stats $same; ai_lab = Stats $aiLab
  lab_lab_values = $labLab; same_crystal_values = $same; ai_lab_values = $aiLab; pct_lab_below_ai_median = $pctBelow
  avg_to_lab = $avgToLab; af_rank = $afRank; to_mean = $toMean; af_mean_rank = $afMeanRank
  af_to_mean = [math]::Round([Lib]::Rmsd($ca['AlphaFold'], $mean, $Taf), 3)
  sensitivity = $sens; correlations = $corr; residues = @($rows); trp185 = $w185; active_site = Stats $site; ser193 = $ser
}
$json | ConvertTo-Json -Depth 6 -Compress | Set-Content (Join-Path $Root "page\ensemble.json") -Encoding UTF8
"Wrote page\ensemble.json"
