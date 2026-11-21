# Compares the AlphaFold prediction of PETase with the X-ray structure 5XJH.
# Superposes C-alpha atoms (Horn quaternion method), reports RMSD, per-residue
# deviation and AlphaFold confidence (pLDDT), and writes results.csv + a superposed model.
param(
  [string]$Data = (Join-Path $PSScriptRoot "data"),
  [string]$Out  = (Join-Path $PSScriptRoot "results"),
  [string]$ModelFile = "AF-A0A0K8P6T7_predicted.pdb",
  [string]$ModelChain = "",
  [int]$MaxResidue = 290   # PETase ends at 290; higher numbers are expression tags
)
New-Item -ItemType Directory -Force $Out | Out-Null

Add-Type -TypeDefinition @'
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

public class Atom { public string Name, Res, Chain, Line; public int Num; public double X, Y, Z, B; }

public static class Prot {
  static CultureInfo C = CultureInfo.InvariantCulture;

  public static List<Atom> Read(string path, string chain) {
    var list = new List<Atom>();
    foreach (var l in File.ReadLines(path)) {
      if (l.StartsWith("ENDMDL")) break;
      if (!l.StartsWith("ATOM  ") || l.Length < 66) continue;
      char alt = l[16];
      if (alt != ' ' && alt != 'A') continue;
      string ch = l.Substring(21, 1);
      if (!string.IsNullOrEmpty(chain) && ch != chain) continue;
      list.Add(new Atom {
        Name = l.Substring(12, 4).Trim(), Res = l.Substring(17, 3), Chain = ch, Line = l,
        Num = int.Parse(l.Substring(22, 4).Trim()),
        X = double.Parse(l.Substring(30, 8), C), Y = double.Parse(l.Substring(38, 8), C),
        Z = double.Parse(l.Substring(46, 8), C), B = double.Parse(l.Substring(60, 6), C) });
    }
    return list;
  }

  // Largest eigenvector of a symmetric 4x4 matrix (Jacobi rotations).
  static double[] TopEigen(double[,] a) {
    int n = 4; var v = new double[n, n];
    for (int i = 0; i < n; i++) v[i, i] = 1;
    for (int sweep = 0; sweep < 100; sweep++) {
      double off = 0;
      for (int p = 0; p < n; p++) for (int q = p + 1; q < n; q++) off += a[p, q] * a[p, q];
      if (off < 1e-20) break;
      for (int p = 0; p < n; p++) for (int q = p + 1; q < n; q++) {
        if (Math.Abs(a[p, q]) < 1e-30) continue;
        double th = (a[q, q] - a[p, p]) / (2 * a[p, q]);
        double t = Math.Sign(th) / (Math.Abs(th) + Math.Sqrt(th * th + 1)); if (th == 0) t = 1;
        double c = 1 / Math.Sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < n; k++) { double akp = a[k, p], akq = a[k, q]; a[k, p] = c * akp - s * akq; a[k, q] = s * akp + c * akq; }
        for (int k = 0; k < n; k++) { double apk = a[p, k], aqk = a[q, k]; a[p, k] = c * apk - s * aqk; a[q, k] = s * apk + c * aqk; }
        for (int k = 0; k < n; k++) { double vkp = v[k, p], vkq = v[k, q]; v[k, p] = c * vkp - s * vkq; v[k, q] = s * vkp + c * vkq; }
      }
    }
    int best = 0; for (int i = 1; i < n; i++) if (a[i, i] > a[best, best]) best = i;
    return new[] { v[0, best], v[1, best], v[2, best], v[3, best] };
  }

  // Returns rotation R (3x3) and translation so that target ~= R*(moving - cm) + ct.
  public static double[] Fit(double[][] mov, double[][] tgt) {
    int n = mov.Length; var cm = new double[3]; var ct = new double[3];
    for (int i = 0; i < n; i++) for (int k = 0; k < 3; k++) { cm[k] += mov[i][k] / n; ct[k] += tgt[i][k] / n; }
    var S = new double[3, 3];
    for (int i = 0; i < n; i++) for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++)
      S[r, c] += (mov[i][r] - cm[r]) * (tgt[i][c] - ct[c]);
    double Sxx=S[0,0],Sxy=S[0,1],Sxz=S[0,2],Syx=S[1,0],Syy=S[1,1],Syz=S[1,2],Szx=S[2,0],Szy=S[2,1],Szz=S[2,2];
    var N = new double[,] {
      { Sxx+Syy+Szz, Syz-Szy,      Szx-Sxz,      Sxy-Syx },
      { Syz-Szy,     Sxx-Syy-Szz,  Sxy+Syx,      Szx+Sxz },
      { Szx-Sxz,     Sxy+Syx,     -Sxx+Syy-Szz,  Syz+Szy },
      { Sxy-Syx,     Szx+Sxz,      Syz+Szy,     -Sxx-Syy+Szz } };
    var q = TopEigen(N); double a=q[0], b=q[1], c2=q[2], d=q[3];
    return new[] {
      a*a+b*b-c2*c2-d*d, 2*(b*c2-a*d),       2*(b*d+a*c2),
      2*(b*c2+a*d),      a*a-b*b+c2*c2-d*d,  2*(c2*d-a*b),
      2*(b*d-a*c2),      2*(c2*d+a*b),       a*a-b*b-c2*c2+d*d,
      cm[0], cm[1], cm[2], ct[0], ct[1], ct[2] };
  }

  public static double[] Apply(double[] T, double[] p) {
    double x = p[0]-T[9], y = p[1]-T[10], z = p[2]-T[11];
    return new[] { T[0]*x+T[1]*y+T[2]*z+T[12], T[3]*x+T[4]*y+T[5]*z+T[13], T[6]*x+T[7]*y+T[8]*z+T[14] };
  }

  public static double Dist(double[] p, double[] q) {
    double dx=p[0]-q[0], dy=p[1]-q[1], dz=p[2]-q[2]; return Math.Sqrt(dx*dx+dy*dy+dz*dz);
  }
  public static double[] P(Atom a) { return new[] { a.X, a.Y, a.Z }; }
}
'@

$model = [Prot]::Read((Join-Path $Data $ModelFile), $ModelChain)
$xray  = [Prot]::Read((Join-Path $Data "5XJH_experimental.pdb"), "A")

$mCA = @{}; foreach ($a in $model) { if ($a.Name -eq "CA" -and $a.Num -le $MaxResidue) { $mCA[$a.Num] = $a } }
$xCA = @{}; foreach ($a in $xray)  { if ($a.Name -eq "CA" -and $a.Num -le $MaxResidue) { $xCA[$a.Num] = $a } }

$shared = @($mCA.Keys | Where-Object { $xCA.ContainsKey($_) } | Sort-Object)
# Crystal residues whose amino acid differs from PETase (expression tag) are left out
$mismatch = @($shared | Where-Object { $mCA[$_].Res -ne $xCA[$_].Res })
$common = @($shared | Where-Object { $mCA[$_].Res -eq $xCA[$_].Res })
foreach ($n in $mismatch) { $xCA.Remove($n) }

function Get-Fit($nums) {
  $mov = [double[][]]@($nums | ForEach-Object { ,[Prot]::P($mCA[$_]) })
  $tgt = [double[][]]@($nums | ForEach-Object { ,[Prot]::P($xCA[$_]) })
  [Prot]::Fit($mov, $tgt)
}
function Get-Rmsd($T, $nums) {
  $s = 0.0
  foreach ($n in $nums) { $d = [Prot]::Dist([Prot]::Apply($T, [Prot]::P($mCA[$n])), [Prot]::P($xCA[$n])); $s += $d * $d }
  [math]::Sqrt($s / $nums.Count)
}

$T = Get-Fit $common
$rmsdAll = Get-Rmsd $T $common

# Per-residue table over the full model
$rows = foreach ($n in ($mCA.Keys | Sort-Object)) {
  $dev = $null
  if ($xCA.ContainsKey($n)) { $dev = [math]::Round([Prot]::Dist([Prot]::Apply($T, [Prot]::P($mCA[$n])), [Prot]::P($xCA[$n])), 2) }
  [pscustomobject]@{ residue = $n; aa = $mCA[$n].Res; plddt = $mCA[$n].B; in_xray = $xCA.ContainsKey($n); deviation_A = $dev }
}
$rows | Export-Csv (Join-Path $Out "per_residue.csv") -NoTypeInformation -Encoding UTF8

$sig   = @($rows | Where-Object { $_.residue -le 29 })
$rest  = @($rows | Where-Object { $_.residue -gt 29 })
$conf  = @($common | Where-Object { $mCA[$_].B -ge 90 })
$low   = @($common | Where-Object { $mCA[$_].B -lt 90 })
$worst = @($rows | Where-Object { $_.deviation_A -ne $null } | Sort-Object deviation_A -Descending | Select-Object -First 8)

# Catalytic triad: all shared heavy atoms after the C-alpha superposition
$triad = foreach ($n in 160, 206, 237) {
  $ma = @($model | Where-Object { $_.Num -eq $n }); $xa = @($xray | Where-Object { $_.Num -eq $n })
  $s = 0.0; $k = 0
  foreach ($a in $ma) { $b = $xa | Where-Object { $_.Name -eq $a.Name } | Select-Object -First 1
    if ($b) { $d = [Prot]::Dist([Prot]::Apply($T, [Prot]::P($a)), [Prot]::P($b)); $s += $d*$d; $k++ } }
  [pscustomobject]@{ residue = $n; model_aa = $ma[0].Res; xray_aa = $xa[0].Res; atoms = $k; rmsd_A = [math]::Round([math]::Sqrt($s / [math]::Max($k,1)), 2) }
}

# Superposed model for viewing alongside 5XJH
$sb = New-Object System.Text.StringBuilder
foreach ($a in $model) {
  $p = [Prot]::Apply($T, [Prot]::P($a))
  $line = $a.Line.Substring(0,30) + ("{0,8:F3}{1,8:F3}{2,8:F3}" -f $p[0],$p[1],$p[2]) + $a.Line.Substring(54)
  [void]$sb.AppendLine($line)
}
[void]$sb.AppendLine("END")
[IO.File]::WriteAllText((Join-Path $Out "AF_superposed_on_5XJH.pdb"), $sb.ToString())

$avg = { param($xs) if ($xs.Count) { [math]::Round((($xs | Measure-Object -Average).Average), 1) } else { "n/a" } }
$summary = [ordered]@{
  model_residues        = $mCA.Count
  xray_residues_chainA  = $xCA.Count
  common_residues       = $common.Count
  xray_range            = "$($common[0])-$($common[-1])"
  residue_name_mismatch = ($mismatch -join ",")
  rmsd_all_CA_A         = [math]::Round($rmsdAll, 2)
  rmsd_plddt90plus_A    = [math]::Round((Get-Rmsd $T $conf), 2)
  n_plddt90plus         = $conf.Count
  rmsd_plddt_below90_A  = $(if ($low.Count) { [math]::Round((Get-Rmsd $T $low), 2) } else { "n/a" })
  n_plddt_below90       = $low.Count
  mean_plddt_signal_1_29 = & $avg @($sig | ForEach-Object plddt)
  mean_plddt_30_290      = & $avg @($rest | ForEach-Object plddt)
  residues_within_1A    = @($rows | Where-Object { $_.deviation_A -ne $null -and $_.deviation_A -le 1 }).Count
}
$summary | ConvertTo-Json | Set-Content (Join-Path $Out "summary.json") -Encoding UTF8
"== Summary =="; $summary | Format-List
"== Biggest differences =="; $worst | Format-Table residue, aa, plddt, deviation_A
"== Catalytic triad (S160, D206, H237) =="; $triad | Format-Table
