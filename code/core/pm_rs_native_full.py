#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pm_rs_native_full.py — RS-native mass tables with universal μ*, top threshold (nf=6),
global uncertainty, consolidated CSV/TeX, and non-circular residuals on demand.


New in this version:
  • --resid-at-mu-star : TeX residuals compare RS to PDG masses transported to the SAME μ* (PDG→μ⋆).
    CSV stays unchanged (no breaking changes).


Everything else remains as in b2:
  • RS-native fixed point with ONE universal μ* (no per-flavor μ* solves)
  • 4L QCD, 2L QED, nf stepping at mc, mb, and mt (nf=6 above mt)
  • α_em policy switch (frozen | leptonic1L) and optional policy band
  • Global uncertainty via MC over {αs(MZ), thresholds, μ*, α policy}
  • Classical ablation via --emit-classical
  • Consolidated outputs (quarks+leptons+bosons+ν)


Usage:
  FAST_RG=1 PM_VERBOSE=1 python3 pm_rs_native_full.py --resid-at-mu-star
  FAST_RG=1 PM_VERBOSE=1 python3 pm_rs_native_full.py --mu-star 182.201 --emit-classical --alpha-policy-band --resid-at-mu-star
"""


import os, csv, math, sys, argparse, random
from typing import Optional


# -------------------- verbosity --------------------
FAST    = os.environ.get("FAST_RG", "0") != "0"
VERBOSE = os.environ.get("PM_VERBOSE", "1") != "0"
def log(msg: str):
    if VERBOSE: print(msg)


# -------------------- math & constants --------------------
PHI   = (1.0 + 5.0**0.5)/2.0
LNPHI = math.log(PHI)


# Anchors / PDG pulls (central)
ALPHA_S_MZ_DEFAULT = 0.1179
ALPHA_mZ_DEFAULT   = 1.0 / 127.955
MZ_GeV             = 91.1876
MC_THR_CENTRAL     = 1.27
MB_THR_CENTRAL     = 4.18
MT_THR_CENTRAL     = 172.76   # used as nf switch to 6 (see comments)


# Envelope (rough, for global bands)
SIG_ALPHAS   = 0.0009
SIG_MC_THR   = 0.02
SIG_MB_THR   = 0.03
SIG_MT_THR   = 0.30


# Bridge μ* relative uncertainty (dominated by G): u_rel(μ*) ≈ 0.5 u_rel(G)
U_REL_G      = 2.2e-5
U_REL_MUSTAR = 0.5 * U_REL_G


# -------------------- RS sector baselines (provenance inline) --------------------
ECOH_eV  = PHI**(-5.0)
ECOH_GeV = ECOH_eV * 1e-9


# Baselines (closed-form; not tuned)
A_U = (2.0 ** (-1)) * ECOH_GeV * (PHI ** 35)
A_D = (2.0 ** (23)) * ECOH_GeV * (PHI ** (-5))


# Integer rungs (from RS constructor; pinned for reproducibility)
R_UP = {"u": 4, "c": 15, "t": 21}
R_DN = {"d": 4, "s": 15, "b": 21}


# Electric charges
Q_UP = {"u": +2.0/3.0, "c": +2.0/3.0, "t": +2.0/3.0}
Q_DN = {"d": -1.0/3.0, "s": -1.0/3.0, "b": -1.0/3.0}


# -------------------- consolidated reference rows --------------------
LEPTON_REFS_GeV = {'e': 0.00051099895, 'mu': 0.1056583755, 'tau': 1.77686}
BOSON_REFS_GeV  = {'W': 80.379, 'Z': 91.1876, 'H': 125.2}


# PDG quark references (for classical control and/or PDG→μ⋆ audit)
PDG_QUARKS = {
    'u': {'mu_ref': 2.0,  'm_ref': 0.00216,  'Q': +2/3, 'scheme': 'MSbar@2GeV'},
    'd': {'mu_ref': 2.0,  'm_ref': 0.00467,  'Q': -1/3, 'scheme': 'MSbar@2GeV'},
    's': {'mu_ref': 2.0,  'm_ref': 0.0930,   'Q': -1/3, 'scheme': 'MSbar@2GeV'},
    'c': {'mu_ref': 1.27, 'm_ref': 1.27,     'Q': +2/3, 'scheme': 'MSbar@$m_c$'},
    'b': {'mu_ref': 4.18, 'm_ref': 4.18,     'Q': -1/3, 'scheme': 'MSbar@$m_b$'},
}


# Dirac NO neutrinos (paper’s RS predictions; Σ≈0.0605 eV)
NU_RS_eV = { 'nu1': 2.0832e-3, 'nu2': 9.0225e-3, 'nu3': 4.9427e-2 }


# -------------------- universal μ* from the bridge --------------------
def universal_mu_star_from_bridge_GeV(alpha_mZ: float) -> float:
    HBAR = 1.054_571_817e-34  # J·s
    C    = 2.997_924_58e8     # m/s
    G    = 6.674_30e-11       # m^3 kg^-1 s^-2
    J_per_GeV = 1.602_176_634e-10
    lP = (HBAR * G / (C**3))**0.5
    tau0    = lP / C
    tau_rec = (2.0 * math.pi / (8.0 * LNPHI)) * tau0
    mu_J    = HBAR / (tau_rec * (PHI**8))
    return mu_J / J_per_GeV


# -------------------- QCD 4L β and 4L γ_m --------------------
ZETA3 = 1.2020569031595942
def beta_coeffs_4L(nf: int):
    b0 = 11.0 - 2.0/3.0 * nf
    b1 = 102.0 - 38.0/3.0 * nf
    b2 = 2857.0/2.0 - 5033.0/18.0 * nf + 325.0/54.0 * nf*nf
    b3 = (149753.0/6.0 + 3564.0*ZETA3) \
         - (1078361.0/162.0 + 6508.0/27.0*ZETA3)*nf \
         + (50065.0/162.0 + 6472.0/81.0*ZETA3)*nf*nf \
         + (1093.0/729.0)*nf*nf*nf
    return b0, b1, b2, b3



def gamma_m_qcd_4L(a_s_over_pi: float, nf: int) -> float:
    a = a_s_over_pi
    g0 = 1.0
    g1 = (101.0/24.0) - (5.0/36.0)*nf
    g2 = (1249.0/64.0) - (277.0/216.0)*nf - (35.0/1296.0)*nf*nf
    g3 = (4603055.0/41472.0 + (41.0/18.0)*ZETA3) \
         - ((91723.0/20736.0) + (11.0/12.0)*ZETA3)*nf \
         + (151.0/162.0)*nf*nf + (5.0/972.0)*nf*nf*nf
    return - (g0*a + g1*a*a + g2*a*a*a + g3*a*a*a*a)


# -------------------- QED 2L γ_m and α_em policies --------------------
def gamma_m_qed_2L(alpha_em: float, Qq: float) -> float:
    a = alpha_em / math.pi
    return - (3.0*(Qq**2)*a + 1.5*(Qq**4)*a*a)



def alpha_em_running(mu: float, policy: str, alpha_mZ: float) -> float:
    if policy == 'frozen':
        return alpha_mZ
    # 1L leptonic running (e, μ, τ)
    m_e, m_mu, m_tau = 0.000510999, 0.105658, 1.77686
    alpha0 = alpha_mZ
    S = 0.0
    for mℓ in (m_e, m_mu, m_tau):
        if mu > mℓ:
            S += math.log(mu / mℓ)
    denom = 1.0 - (2.0 * alpha0 / (3.0 * math.pi)) * S
    return alpha0 / denom if denom > 0 else alpha0 * 1e6


# -------------------- nf stepping and αs(μ) with mt threshold --------------------
def run_alpha_s_4L(mu_lo: float, mu_hi: float, alpha_s_lo: float, nf: int) -> float:
    def step_a4(a4, dln):
        b0,b1,b2,b3 = beta_coeffs_4L(nf)
        def beta(a): return - (b0*a*a + b1*a*a*a + b2*a*a*a*a + b3*a*a*a*a*a)
        k1 = beta(a4)
        k2 = beta(a4 + 0.5*dln*k1)
        k3 = beta(a4 + 0.5*dln*k2)
        k4 = beta(a4 + dln*k3)
        return a4 + (dln/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    if mu_lo <= 0: mu_lo = 1e-12
    if mu_hi <= 0: mu_hi = 1e-12
    if mu_hi == mu_lo: return alpha_s_lo
    a4 = alpha_s_lo / (4.0*math.pi)
    step_mult = 180 if FAST else 900
    base_steps = 90 if FAST else 300
    n = max(base_steps, int(abs(math.log(mu_hi/mu_lo))*step_mult))
    sign = 1.0 if mu_hi>mu_lo else -1.0
    dln  = sign*abs(math.log(mu_hi/mu_lo))/n
    for _ in range(n):
        a4 = step_a4(a4, dln)
        a4 = max(min(a4, 0.5), 1e-8)
    return (4.0*math.pi)*a4



def alpha_s_MSbar(mu: float, alpha_s_MZ: float, mc: float, mb: float, mt: float) -> float:
    if mu == MZ_GeV:
        return alpha_s_MZ
    as_cur = alpha_s_MZ
    legs = []
    if mu < MZ_GeV:
        legs.append((MZ_GeV, max(mu, mb), 5))
        if mu < mb:
            legs.append((mb, max(mu, mc), 4))
            if mu < mc:
                legs.append((mc, mu, 3))
    else:
        if mu <= mt:
            legs.append((MZ_GeV, mu, 5))
        else:
            legs.append((MZ_GeV, mt, 5))
            # αs continuity at mt (tree-level matching). TODO: higher-order matching if desired.
            legs.append((mt, mu, 6))
    for lo, hi, nf in legs:
        as_cur = run_alpha_s_4L(lo, hi, as_cur, nf)
    return as_cur


# -------------------- RG integrals & fixed point --------------------
def ln_R_qcd(mu_lo: float, mu_hi: float, alpha_s_MZ: float, mc: float, mb: float, mt: float) -> float:
    mu_lo = max(mu_lo, 1e-12); mu_hi = max(mu_hi, 1e-12)
    if mu_lo == mu_hi: return 0.0
    x0, x1 = math.log(mu_lo), math.log(mu_hi)
    base = 200 if FAST else 900
    mult = 500 if FAST else 1400
    n = max(base, int(abs(x1-x0)*mult))
    s = 0.0
    dx = (x1 - x0)/n
    for k in range(n):
        xm  = x0 + (k + 0.5)*dx
        mu  = math.exp(xm)
        a_s = alpha_s_MSbar(mu, alpha_s_MZ, mc, mb, mt) / math.pi
        # nf appears only inside γ_m here (captured by αs profile); included for clarity
        nf = 3 if mu < mc else (4 if mu < mb else (5 if mu < mt else 6))
        s += gamma_m_qcd_4L(a_s, nf) * dx
    return s



def ln_R_qed(mu_lo: float, mu_hi: float, Qq: float, alpha_policy: str, alpha_mZ: float) -> float:
    if mu_lo == mu_hi: return 0.0
    x0, x1 = math.log(max(mu_lo,1e-12)), math.log(max(mu_hi,1e-12))
    n = 60 if FAST else 200
    s = 0.0
    dx = (x1 - x0)/n
    for k in range(n):
        xm  = x0 + (k + 0.5)*dx
        mu  = math.exp(xm)
        alpha = alpha_em_running(mu, alpha_policy, alpha_mZ)
        s += gamma_m_qed_2L(alpha, Qq) * dx
    return s



def residue_total(mu_star: float, mu: float, Qq: float, alpha_s_MZ: float, mc: float, mb: float, mt: float,
                  alpha_policy: str, alpha_mZ: float) -> float:
    return ( ln_R_qcd(mu_star, mu, alpha_s_MZ, mc, mb, mt)
           + ln_R_qed(mu_star, mu, Qq, alpha_policy, alpha_mZ) ) / LNPHI



def fixed_point_mass(A_B: float, r_i: int, mu_star: float, Qq: float, alpha_s_MZ: float, mc: float, mb: float, mt: float,
                     alpha_policy: str, alpha_mZ: float, tol: float=1e-9, itmax: int=300) -> float:
    m = max(A_B * (PHI ** r_i), 1e-9)
    for _ in range(itmax):
        mu_eval = max(m, 1e-12)
        f = residue_total(mu_star, mu_eval, Qq, alpha_s_MZ, mc, mb, mt, alpha_policy, alpha_mZ)
        expo = max(min(r_i + f, 80.0), -80.0)
        m_new = A_B * (PHI ** expo)
        if abs(m_new - m)/max(m,1e-16) < tol:
            return m_new
        m = 0.5*m + 0.5*m_new
    return m


# -------------------- predictions & classical control --------------------
def quark_predictions_rs(mu_star: float, alpha_s_MZ: float, mc: float, mb: float, mt: float,
                         alpha_policy: str, alpha_mZ: float) -> dict:
    out = {}
    out['d'] = fixed_point_mass(A_D, R_DN['d'], mu_star, Q_DN['d'], alpha_s_MZ, mc, mb, mt, alpha_policy, alpha_mZ)
    out['s'] = fixed_point_mass(A_D, R_DN['s'], mu_star, Q_DN['s'], alpha_s_MZ, mc, mb, mt, alpha_policy, alpha_mZ)
    out['u'] = fixed_point_mass(A_U, R_UP['u'], mu_star, Q_UP['u'], alpha_s_MZ, mc, mb, mt, alpha_policy, alpha_mZ)
    out['c'] = fixed_point_mass(A_U, R_UP['c'], mu_star, Q_UP['c'], alpha_s_MZ, mc, mb, mt, alpha_policy, alpha_mZ)
    out['b'] = fixed_point_mass(A_D, R_DN['b'], mu_star, Q_DN['b'], alpha_s_MZ, mc, mb, mt, alpha_policy, alpha_mZ)
    return out



def quark_predictions_classical(mu_star: float, alpha_s_MZ: float, mc: float, mb: float, mt: float,
                                alpha_policy: str, alpha_mZ: float) -> dict:
    out = {}
    for sp, d in PDG_QUARKS.items():
        mu_ref, m_ref, Qq = d['mu_ref'], d['m_ref'], d['Q']
        lr_qcd = ln_R_qcd(mu_ref, mu_star, alpha_s_MZ, mc, mb, mt)
        lr_qed = ln_R_qed(mu_ref, mu_star, Qq, alpha_policy, alpha_mZ)
        out[sp] = m_ref * math.exp(lr_qcd + lr_qed)
    return out


# -------------------- uncertainty MC --------------------
def mc_uncertainty(mu_star_c: float, alpha_s_c: float, alpha_policy: str, alpha_mZ: float,
                   N: int=200, policy_band: bool=False) -> dict:
    vals = { k: [] for k in ('d','s','u','c','b') }
    # α policy band proxy (central half-diff)
    delta_policy = {k:0.0 for k in vals}
    if policy_band:
        pf = quark_predictions_rs(mu_star_c, alpha_s_c, MC_THR_CENTRAL, MB_THR_CENTRAL, MT_THR_CENTRAL, 'frozen',      alpha_mZ)
        pl = quark_predictions_rs(mu_star_c, alpha_s_c, MC_THR_CENTRAL, MB_THR_CENTRAL, MT_THR_CENTRAL, 'leptonic1L', alpha_mZ)
        delta_policy = {k:0.5*abs(pl[k]-pf[k]) for k in vals}
    for _ in range(max(1,N)):
        alphas = random.gauss(alpha_s_c, SIG_ALPHAS)
        mc_thr = max(0.3, random.gauss(MC_THR_CENTRAL, SIG_MC_THR))
        mb_thr = max(mc_thr+0.3, random.gauss(MB_THR_CENTRAL, SIG_MB_THR))
        mt_thr = max(mb_thr+20.0, random.gauss(MT_THR_CENTRAL, SIG_MT_THR))
        mu_fac = random.gauss(1.0, U_REL_MUSTAR)
        mu_star= max(1e-6, mu_star_c * mu_fac)
        p = quark_predictions_rs(mu_star, alphas, mc_thr, mb_thr, mt_thr, alpha_policy, alpha_mZ)
        for k in vals:
            jitter = delta_policy[k] * (1 if random.random()<0.5 else -1)
            vals[k].append(max(0.0, p[k] + jitter))
    sig = {}
    for k, arr in vals.items():
        m = sum(arr)/len(arr)
        v = sum((x-m)**2 for x in arr)/(len(arr)-1 if len(arr)>1 else 1)
        sig[k] = math.sqrt(max(0.0, v))
    return sig


# -------------------- writers --------------------
def ensure_dirs():
    os.makedirs(os.path.join("out","csv"), exist_ok=True)
    os.makedirs(os.path.join("out","tex"), exist_ok=True)



def fmt(x: float, nd: int=6) -> str:
    if x is None: return "--"
    ax = abs(x)
    if (ax != 0 and ax < 1e-6) or ax >= 1e6: return f"{x:.3e}"
    return f"{x:.{nd}f}"



def write_consolidated_csv(path: str, mu_star: float, quarks: dict, sig: Optional[dict], mode: str):
    ensure_dirs()
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Species","Predicted[GeV]","±(global)","Reference[GeV]","Scheme/Scale","Mode","μ*[GeV]"])
        for sp in ['e','mu','tau']:
            w.writerow([sp,"--","--",fmt(LEPTON_REFS_GeV[sp]),"pole",mode,fmt(mu_star)])
        for sp in ['nu1','nu2','nu3']:
            w.writerow([sp,fmt(NU_RS_eV[sp]*1e-9),"--","--","",mode,fmt(mu_star)])
        for sp in ['d','s','u','c','b']:
            ref = PDG_QUARKS[sp]['m_ref']
            sch = PDG_QUARKS[sp]['scheme']
            w.writerow([sp,fmt(quarks[sp]),fmt((sig or {}).get(sp,float('nan'))),fmt(ref),sch,mode,fmt(mu_star)])
        for sp in ['W','Z','H']:
            w.writerow([sp,"--","--",fmt(BOSON_REFS_GeV[sp]),"pole",mode,fmt(mu_star)])



def write_consolidated_tex(path: str, mu_star: float, quarks: dict, sig: Optional[dict], mode: str,
                           include_nu: bool=True, resid_at_mu_star: bool=False, pdg_to_mu_star: Optional[dict]=None):
    """
    If resid_at_mu_star=True and pdg_to_mu_star provided, show 'Reference' column as PDG→μ⋆ and
    compute Residual vs those values (non-circular). Else use PDG scheme refs as before.
    """
    ensure_dirs()
    with open(path, "w", encoding="utf-8") as f:
        f.write("\\begin{center}\n\\begin{tabular}{lrrlc}\n\\toprule\n")
        if resid_at_mu_star:
            f.write("Species & Predicted [GeV] & PDG$\\to\\,\\mu_\\star$ [GeV] & Scheme/Scale & Residual \\\\\n")
        else:
            f.write("Species & Predicted [GeV] & Reference [GeV] & Scheme/Scale & Residual \\\\\n")
        f.write("\\midrule\n")
        # leptons (refs only)
        for sp in ['e','mu','tau']:
            f.write(f"{sp} & -- & {fmt(LEPTON_REFS_GeV[sp])} & pole & -- \\\\\n")
        # neutrinos (optional)
        if include_nu:
            for sp in ['nu1','nu2','nu3']:
                f.write(f"{sp} & {fmt(NU_RS_eV[sp]*1e-9)} & -- &  & -- \\\\\n")
        # quarks
        for sp in ['d','s','u','c','b']:
            if resid_at_mu_star and pdg_to_mu_star is not None:
                ref_val = pdg_to_mu_star[sp]
                sch     = "PDG$\\to\\,\\mu_\\star$"
            else:
                ref_val = PDG_QUARKS[sp]['m_ref']
                sch     = PDG_QUARKS[sp]['scheme']
            resid = (quarks[sp]-ref_val)/ref_val if ref_val else float('nan')
            pred_s = fmt(quarks[sp])
            if sig is not None and sp in sig:
                pred_s = f"{pred_s} $\\pm$ {fmt(sig[sp])}"
            f.write(f"{sp} & {pred_s} & {fmt(ref_val)} & {sch} & {resid:.3e} \\\\\n")
        # bosons (refs only)
        for sp in ['W','Z','H']:
            f.write(f"{sp} & -- & {fmt(BOSON_REFS_GeV[sp])} & pole & -- \\\\\n")
        f.write("\\bottomrule\n\\par\\vspace{2pt}\n")
        f.write(f"\\emph{{Mode:}} {mode}; $\\;\\mu_\\star={fmt(mu_star)}$ GeV.\\end{tabular}\n\\end{center}\n")


# -------------------- CLI --------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="RS-native mass tables with universal μ*, nf=6 top threshold, consolidated outputs.")
    ap.add_argument("--alpha-s", type=float, default=ALPHA_S_MZ_DEFAULT, help="αs(MZ) central value (default 0.1179)")
    ap.add_argument("--alpha-policy", choices=("frozen","leptonic1L"), default="frozen", help="α_em(μ) policy (default: frozen)")
    ap.add_argument("--alpha-policy-band", action="store_true", help="Add α_em policy delta into global uncertainty band")
    ap.add_argument("--mu-star", type=float, default=None, help="Override universal μ* in GeV (else bridge)")
    ap.add_argument("--mc-samples", type=int, default=200, help="MC samples for global band (default 200)")
    ap.add_argument("--no-nu", action="store_true", help="Omit neutrino rows in consolidated outputs")
    ap.add_argument("--emit-classical", action="store_true", help="Also emit classical control CSV/TeX")
    ap.add_argument("--resid-at-mu-star", action="store_true", help="TeX residuals vs PDG masses transported to the SAME μ* (non-circular)")
    args = ap.parse_args(argv)


    alpha_s_MZ  = float(args.alpha_s)
    alpha_mZ    = ALPHA_mZ_DEFAULT
    alpha_policy= args.alpha_policy


    # μ* central
    if args.mu_star is not None:
        mu_star = max(1e-6, float(args.mu_star)); mu_src = "user"
    else:
        mu_star = universal_mu_star_from_bridge_GeV(alpha_mZ); mu_src = "bridge"


    # thresholds central
    mc, mb, mt = MC_THR_CENTRAL, MB_THR_CENTRAL, MT_THR_CENTRAL
    log(f"[INFO] FAST_RG={'1' if FAST else '0'}  αs(MZ)={alpha_s_MZ:.6f}  α_policy={alpha_policy}  μ*({mu_src})={mu_star:.6f} GeV  thresholds=(mc={mc}, mb={mb}, mt={mt})")


    # RS predictions (central)
    quarks_rs = quark_predictions_rs(mu_star, alpha_s_MZ, mc, mb, mt, alpha_policy, alpha_mZ)
    # Global band via MC
    sig_rs    = mc_uncertainty(mu_star, alpha_s_MZ, alpha_policy, alpha_mZ, N=max(1,args.mc_samples), policy_band=args.alpha_policy_band)


    # PDG→μ⋆ map (only if needed for non-circular residuals)
    pdg_muStar = None
    if args.resid_at_mu_star:
        pdg_muStar = quark_predictions_classical(mu_star, alpha_s_MZ, mc, mb, mt, alpha_policy, alpha_mZ)


    # RS consolidated
    ensure_dirs()
    write_consolidated_csv(os.path.join("out","csv","all_masses_rs.csv"), mu_star, quarks_rs, sig_rs, "RS-native")
    write_consolidated_tex(os.path.join("out","tex","all_masses_rs.tex"), mu_star, quarks_rs, sig_rs, "RS-native",
                           include_nu=(not args.no_nu), resid_at_mu_star=args.resid_at_mu_star, pdg_to_mu_star=pdg_muStar)
    log("[OK] wrote RS consolidated CSV/TeX")


    # Classical control (optional)
    if args.emit_classical:
        quarks_cl = quark_predictions_classical(mu_star, alpha_s_MZ, mc, mb, mt, alpha_policy, alpha_mZ)
        sig_cl    = mc_uncertainty(mu_star, alpha_s_MZ, alpha_policy, alpha_mZ, N=max(1,args.mc_samples), policy_band=args.alpha_policy_band)
        write_consolidated_csv(os.path.join("out","csv","all_masses_classical.csv"), mu_star, quarks_cl, sig_cl, "Classical control")
        write_consolidated_tex(os.path.join("out","tex","all_masses_classical.tex"), mu_star, quarks_cl, sig_cl, "Classical control",
                               include_nu=(not args.no_nu), resid_at_mu_star=args.resid_at_mu_star, pdg_to_mu_star=pdg_muStar)
        log("[OK] wrote Classical control CSV/TeX")


    # Console skim
    print("\n[RS-native @ universal μ*]")
    print(f"μ* = {mu_star:.6g} GeV   (source: {mu_src})")
    for sp in ['d','s','u','c','b']:
        print(f"{sp}: {quarks_rs[sp]:.9g} GeV   (± {sig_rs[sp]:.3g})")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr); sys.exit(1)
