// Shared structure-analysis helpers for the PETase study (loaded with Add-Type -Path).
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

public class PAtom { public string Name, Res, Chain; public int Num; public string Ins; public double X, Y, Z, B; }

public class PResidue {
  public string Res; public int Num; public string Key;
  public Dictionary<string, PAtom> Atoms = new Dictionary<string, PAtom>();
  public PAtom CA { get { PAtom a; return Atoms.TryGetValue("CA", out a) ? a : null; } }
}

public static class Lib {
  static readonly CultureInfo C = CultureInfo.InvariantCulture;
  static readonly Dictionary<string, char> Three = new Dictionary<string, char> {
    {"ALA",'A'},{"ARG",'R'},{"ASN",'N'},{"ASP",'D'},{"CYS",'C'},{"GLN",'Q'},{"GLU",'E'},{"GLY",'G'},{"HIS",'H'},{"ILE",'I'},
    {"LEU",'L'},{"LYS",'K'},{"MET",'M'},{"PHE",'F'},{"PRO",'P'},{"SER",'S'},{"THR",'T'},{"TRP",'W'},{"TYR",'Y'},{"VAL",'V'},{"MSE",'M'} };

  public static char One(string r) { char c; return Three.TryGetValue(r, out c) ? c : 'X'; }

  public static List<string> Chains(string path) {
    var s = new List<string>();
    foreach (var l in File.ReadLines(path)) {
      if (l.StartsWith("ENDMDL")) break;
      if (l.StartsWith("ATOM  ") && l.Length > 21) { var c = l.Substring(21, 1); if (!s.Contains(c)) s.Add(c); }
    }
    return s;
  }

  // Residues of one chain, in file order, protein atoms only (first model, altloc blank or A).
  public static List<PResidue> Residues(string path, string chain) {
    var list = new List<PResidue>(); PResidue cur = null;
    foreach (var l in File.ReadLines(path)) {
      if (l.StartsWith("ENDMDL")) break;
      bool atom = l.StartsWith("ATOM  ") || (l.StartsWith("HETATM") && l.Substring(17, 3) == "MSE");
      if (!atom || l.Length < 66) continue;
      char alt = l[16]; if (alt != ' ' && alt != 'A') continue;
      string ch = l.Substring(21, 1); if (!string.IsNullOrEmpty(chain) && ch != chain) continue;
      var a = new PAtom { Name = l.Substring(12, 4).Trim(), Res = l.Substring(17, 3), Chain = ch,
        Num = int.Parse(l.Substring(22, 4).Trim()), Ins = l.Substring(26, 1),
        X = double.Parse(l.Substring(30, 8), C), Y = double.Parse(l.Substring(38, 8), C),
        Z = double.Parse(l.Substring(46, 8), C), B = double.Parse(l.Substring(60, 6), C) };
      string key = ch + a.Num + a.Ins;
      if (cur == null || cur.Key != key) { cur = new PResidue { Res = a.Res, Num = a.Num, Key = key }; list.Add(cur); }
      if (!cur.Atoms.ContainsKey(a.Name)) cur.Atoms[a.Name] = a;
    }
    return list.Where(r => r.CA != null).ToList();
  }

  public static string Seq(List<PResidue> rs) { return new string(rs.Select(r => One(r.Res)).ToArray()); }

  // Global alignment (Needleman-Wunsch, BLOSUM-free: match +5, mismatch -3, gap -6 open-free linear).
  // Returns for each position i in a the index in b it aligns to, or -1.
  public static int[] Align(string a, string b, double match = 5, double mismatch = -3, double gap = -6) {
    int n = a.Length, m = b.Length; var S = new double[n + 1, m + 1]; var P = new byte[n + 1, m + 1];
    for (int i = 1; i <= n; i++) { S[i, 0] = 0; P[i, 0] = 1; }           // free end gaps
    for (int j = 1; j <= m; j++) { S[0, j] = 0; P[0, j] = 2; }
    for (int i = 1; i <= n; i++) for (int j = 1; j <= m; j++) {
      double d = S[i - 1, j - 1] + (a[i - 1] == b[j - 1] ? match : mismatch);
      double u = S[i - 1, j] + (j == m ? 0 : gap);
      double l = S[i, j - 1] + (i == n ? 0 : gap);
      if (d >= u && d >= l) { S[i, j] = d; P[i, j] = 0; } else if (u >= l) { S[i, j] = u; P[i, j] = 1; } else { S[i, j] = l; P[i, j] = 2; }
    }
    var map = Enumerable.Repeat(-1, n).ToArray(); int x = n, y = m;
    while (x > 0 && y > 0) {
      if (P[x, y] == 0) { map[x - 1] = y - 1; x--; y--; } else if (P[x, y] == 1) x--; else y--;
    }
    return map;
  }

  public static double Identity(string a, string b, int[] map) {
    int same = 0, aligned = 0;
    for (int i = 0; i < a.Length; i++) if (map[i] >= 0) { aligned++; if (a[i] == b[map[i]]) same++; }
    return aligned == 0 ? 0 : 100.0 * same / aligned;
  }

  static double[] TopEigen(double[,] a) {
    int n = 4; var v = new double[n, n]; for (int i = 0; i < n; i++) v[i, i] = 1;
    for (int sweep = 0; sweep < 100; sweep++) {
      double off = 0; for (int p = 0; p < n; p++) for (int q = p + 1; q < n; q++) off += a[p, q] * a[p, q];
      if (off < 1e-20) break;
      for (int p = 0; p < n; p++) for (int q = p + 1; q < n; q++) {
        if (Math.Abs(a[p, q]) < 1e-30) continue;
        double th = (a[q, q] - a[p, p]) / (2 * a[p, q]);
        double t = th == 0 ? 1 : Math.Sign(th) / (Math.Abs(th) + Math.Sqrt(th * th + 1));
        double c = 1 / Math.Sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < n; k++) { double x1 = a[k, p], x2 = a[k, q]; a[k, p] = c * x1 - s * x2; a[k, q] = s * x1 + c * x2; }
        for (int k = 0; k < n; k++) { double x1 = a[p, k], x2 = a[q, k]; a[p, k] = c * x1 - s * x2; a[q, k] = s * x1 + c * x2; }
        for (int k = 0; k < n; k++) { double x1 = v[k, p], x2 = v[k, q]; v[k, p] = c * x1 - s * x2; v[k, q] = s * x1 + c * x2; }
      }
    }
    int best = 0; for (int i = 1; i < n; i++) if (a[i, i] > a[best, best]) best = i;
    return new[] { v[0, best], v[1, best], v[2, best], v[3, best] };
  }

  // Transform T maps "mov" coordinates onto "tgt".
  public static double[] Fit(List<double[]> mov, List<double[]> tgt) {
    int n = mov.Count; var cm = new double[3]; var ct = new double[3];
    for (int i = 0; i < n; i++) for (int k = 0; k < 3; k++) { cm[k] += mov[i][k] / n; ct[k] += tgt[i][k] / n; }
    var S = new double[3, 3];
    for (int i = 0; i < n; i++) for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) S[r, c] += (mov[i][r] - cm[r]) * (tgt[i][c] - ct[c]);
    double Sxx=S[0,0],Sxy=S[0,1],Sxz=S[0,2],Syx=S[1,0],Syy=S[1,1],Syz=S[1,2],Szx=S[2,0],Szy=S[2,1],Szz=S[2,2];
    var N = new double[,] {
      { Sxx+Syy+Szz, Syz-Szy, Szx-Sxz, Sxy-Syx }, { Syz-Szy, Sxx-Syy-Szz, Sxy+Syx, Szx+Sxz },
      { Szx-Sxz, Sxy+Syx, -Sxx+Syy-Szz, Syz+Szy }, { Sxy-Syx, Szx+Sxz, Syz+Szy, -Sxx-Syy+Szz } };
    var q = TopEigen(N); double a=q[0], b=q[1], c2=q[2], d=q[3];
    return new[] { a*a+b*b-c2*c2-d*d, 2*(b*c2-a*d), 2*(b*d+a*c2), 2*(b*c2+a*d), a*a-b*b+c2*c2-d*d, 2*(c2*d-a*b),
                   2*(b*d-a*c2), 2*(c2*d+a*b), a*a-b*b-c2*c2+d*d, cm[0], cm[1], cm[2], ct[0], ct[1], ct[2] };
  }
  public static double[] Apply(double[] T, double[] p) {
    double x = p[0]-T[9], y = p[1]-T[10], z = p[2]-T[11];
    return new[] { T[0]*x+T[1]*y+T[2]*z+T[12], T[3]*x+T[4]*y+T[5]*z+T[13], T[6]*x+T[7]*y+T[8]*z+T[14] };
  }
  public static double[] P(PAtom a) { return new[] { a.X, a.Y, a.Z }; }
  public static double Dist(double[] p, double[] q) { double dx=p[0]-q[0], dy=p[1]-q[1], dz=p[2]-q[2]; return Math.Sqrt(dx*dx+dy*dy+dz*dz); }

  public static double Rmsd(List<double[]> mov, List<double[]> tgt, double[] T) {
    double s = 0; for (int i = 0; i < mov.Count; i++) { var d = Dist(Apply(T, mov[i]), tgt[i]); s += d * d; }
    return Math.Sqrt(s / mov.Count);
  }

  // Superpose, then repeatedly drop pairs > cutoff and refit (for comparing different proteins' cores).
  public static double[] FitCore(List<double[]> mov, List<double[]> tgt, double cutoff, out int kept, out double rmsdCore) {
    var idx = Enumerable.Range(0, mov.Count).ToList(); double[] T = Fit(mov, tgt);
    for (int it = 0; it < 10; it++) {
      var keep = idx.Where(i => Dist(Apply(T, mov[i]), tgt[i]) <= cutoff).ToList();
      if (keep.Count < 3 || keep.Count == idx.Count) { idx = keep.Count < 3 ? idx : keep; break; }
      idx = keep; T = Fit(idx.Select(i => mov[i]).ToList(), idx.Select(i => tgt[i]).ToList());
    }
    kept = idx.Count; rmsdCore = Rmsd(idx.Select(i => mov[i]).ToList(), idx.Select(i => tgt[i]).ToList(), T);
    return T;
  }

  // TM-score (Zhang & Skolnick 2004) for an existing residue pairing, normalised by lengthNorm.
  // The superposition is searched over several pruning cutoffs and the best score kept; the residue
  // pairing itself comes from the sequence alignment, so this is a lower bound on what TM-align would report.
  public static double TmScore(List<double[]> mov, List<double[]> tgt, int lengthNorm) {
    double d0 = 1.24 * Math.Pow(Math.Max(lengthNorm, 19) - 15, 1.0 / 3) - 1.8; double best = 0;
    foreach (var cut in new[] { 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 100.0 }) {
      int k; double r; var T = FitCore(mov, tgt, cut, out k, out r);
      double s = 0; for (int i = 0; i < mov.Count; i++) { double d = Dist(Apply(T, mov[i]), tgt[i]); s += 1 / (1 + (d / d0) * (d / d0)); }
      best = Math.Max(best, s / lengthNorm);
    }
    return best;
  }

  public static double Dihedral(double[] p0, double[] p1, double[] p2, double[] p3) {
    double[] b0 = Sub(p0, p1), b1 = Sub(p2, p1), b2 = Sub(p3, p2);
    double n1 = Norm(b1); var b1n = new[] { b1[0]/n1, b1[1]/n1, b1[2]/n1 };
    var v = Sub(b0, Mul(b1n, DotP(b0, b1n))); var w = Sub(b2, Mul(b1n, DotP(b2, b1n)));
    double x = DotP(v, w); double y = DotP(Cross(b1n, v), w);
    return Math.Atan2(y, x) * 180 / Math.PI;
  }
  static double[] Sub(double[] a, double[] b) { return new[] { a[0]-b[0], a[1]-b[1], a[2]-b[2] }; }
  static double[] Mul(double[] a, double k) { return new[] { a[0]*k, a[1]*k, a[2]*k }; }
  static double DotP(double[] a, double[] b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
  static double Norm(double[] a) { return Math.Sqrt(DotP(a, a)); }
  static double[] Cross(double[] a, double[] b) { return new[] { a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] }; }

  public static double Pearson(double[] x, double[] y) {
    double mx = x.Average(), my = y.Average(), sxy = 0, sxx = 0, syy = 0;
    for (int i = 0; i < x.Length; i++) { sxy += (x[i]-mx)*(y[i]-my); sxx += (x[i]-mx)*(x[i]-mx); syy += (y[i]-my)*(y[i]-my); }
    return sxy / Math.Sqrt(sxx * syy);
  }
  static double[] Ranks(double[] x) {
    var order = Enumerable.Range(0, x.Length).OrderBy(i => x[i]).ToArray(); var r = new double[x.Length];
    for (int i = 0; i < order.Length; ) { int j = i; while (j + 1 < order.Length && x[order[j + 1]] == x[order[i]]) j++;
      double avg = (i + j) / 2.0 + 1; for (int k = i; k <= j; k++) r[order[k]] = avg; i = j + 1; }
    return r;
  }
  public static double Spearman(double[] x, double[] y) { return Pearson(Ranks(x), Ranks(y)); }
}
