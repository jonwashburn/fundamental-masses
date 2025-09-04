#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build a single consolidated table of particle masses (predicted vs reference) for the paper.

Sources:
- Quark predicted masses at μ* from out/csv/quark_masses_muStar.csv
- Charged leptons (pred=ref; pole masses)
- Neutrino Dirac masses (pred only; from manuscript numbers)
- Bosons (ref only; pole masses)

Writes:
  out/tex/all_masses.tex (booktabs table)
"""

import os, csv, math

OUT_TEX = os.path.join('out', 'tex', 'all_masses.tex')
OUT_TEX_MINI = os.path.join('out', 'tex', 'all_masses_mini.tex')
QUARK_CSV = os.path.join('out', 'csv', 'quark_masses_muStar.csv')

def fmt(x, nd=6):
    if x is None:
        return "--"
    ax = abs(x)
    if (ax != 0 and ax < 1e-6) or ax >= 1e6:
        return f"{x:.3e}"
    return f"{x:.{nd}f}"

def main() -> int:
    os.makedirs(os.path.dirname(OUT_TEX), exist_ok=True)

    rows = []

    # --- Minimal runner to compute quark sigmas (4L QCD + 2L QED) ---
    PHI   = (1.0 + 5.0**0.5)/2.0
    ALPHA_mZ = 1.0 / 127.955
    ALPHA_S_MZ = 0.1179
    MZ_GeV   = 91.1876
    MC_THR   = 1.27
    MB_THR   = 4.18
    MU_FLOOR = 1.0

    # Sector constants and integers
    ECOH_eV  = PHI**(-5.0)
    ECOH_GeV = ECOH_eV * 1e-9
    A_U = (2.0 ** (-1)) * ECOH_GeV * (PHI ** 35)
    A_D = (2.0 ** (23)) * ECOH_GeV * (PHI ** (-5))
    R_UP = {'u': 4, 'c': 15, 't': 21}
    R_DN = {'d': 4, 's': 15, 'b': 21}
    Q_UP = {'u': +2.0/3.0, 'c': +2.0/3.0, 't': +2.0/3.0}
    Q_DN = {'d': -1.0/3.0, 's': -1.0/3.0, 'b': -1.0/3.0}
    MU_STAR = {'u': 2.0, 'd': 2.0, 's': 2.0, 'c': 1.27, 'b': 4.18}

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

    def run_alpha_s_4L(mu_lo, mu_hi, alpha_s_lo, nf):
        def step_as(a4, dlnmu):
            b0,b1,b2,b3 = beta_coeffs_4L(nf)
            def beta_a4(a):
                return - (b0*a*a + b1*a*a*a + b2*a*a*a*a + b3*a*a*a*a*a)
            k1 = beta_a4(a4)
            k2 = beta_a4(a4 + 0.5*dlnmu*k1)
            k3 = beta_a4(a4 + 0.5*dlnmu*k2)
            k4 = beta_a4(a4 + dlnmu*k3)
            return a4 + (dlnmu/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        a4 = alpha_s_lo / (4.0*math.pi)
        nstep = max(200, int(abs(math.log(mu_hi/mu_lo))*800))
        sgn   = 1.0 if mu_hi>mu_lo else -1.0
        dln   = sgn*abs(math.log(mu_hi/mu_lo))/nstep
        mu    = mu_lo
        for _ in range(nstep):
            mu = mu*math.exp(dln)
            a4 = step_as(a4, dln)
            a4 = max(min(a4, 0.5), 1e-6)
        return (4.0*math.pi)*a4

    def nf_profile(mu: float) -> int:
        if mu < MC_THR:
            return 3
        elif mu < MB_THR:
            return 4
        else:
            return 5

    def alpha_s_MSbar(mu: float, alpha_s_mZ: float = None) -> float:
        """αs(μ) with fixed decoupling at m_c, m_b; nf segmentation matches nf_profile."""
        if alpha_s_mZ is None:
            alpha_s_mZ = ALPHA_S_MZ
        if mu == MZ_GeV:
            return alpha_s_mZ
        def run_between(mu_lo, mu_hi, as_lo, nf):
            return run_alpha_s_4L(mu_lo, mu_hi, as_lo, nf)
        as_cur, mu_cur = alpha_s_mZ, MZ_GeV
        legs = []
        if mu < MZ_GeV:
            legs.append((MZ_GeV, MB_THR, 5))
            if mu < MB_THR:
                top = max(mu, MC_THR)
                legs.append((MB_THR, top, 4))
                if mu < MC_THR:
                    legs.append((top, mu, 3))
        else:
            legs.append((MZ_GeV, mu, 5))
        for lo, hi, nf in legs:
            as_cur = run_between(lo, hi, as_cur, nf)
            mu_cur = hi
        return as_cur

    def gamma_m_qcd_4L(a_s: float, nf: int) -> float:
        g0 = 1.0
        g1 = (101.0/24.0) - (5.0/36.0)*nf
        g2 = (1249.0/64.0) - (277.0/216.0)*nf - (35.0/1296.0)*nf*nf
        g3 = (4603055.0/41472.0 + (41.0/18.0)*ZETA3) \
             - ((91723.0/20736.0) + (11.0/12.0)*ZETA3)*nf \
             + (151.0/162.0)*nf*nf + (5.0/972.0)*nf*nf*nf
        return - (g0*a_s + g1*a_s*a_s + g2*a_s*a_s*a_s + g3*a_s*a_s*a_s*a_s)

    def gamma_m_qed_2L(alpha_em: float, Qq: float) -> float:
        a = alpha_em / math.pi
        return - (3.0*Qq*Qq*a + 1.5*(Qq**4)*a*a)

    def ln_R_qcd(mu_lo: float, mu_hi: float) -> float:
        mu_lo = max(mu_lo, MU_FLOOR)
        mu_hi = max(mu_hi, MU_FLOOR)
        if mu_lo == mu_hi:
            return 0.0
        x0, x1 = math.log(mu_lo), math.log(mu_hi)
        nstep  = max(800, int(abs(x1 - x0) * 1600))
        s = 0.0
        for k in range(nstep):
            xm  = x0 + (x1 - x0) * (k + 0.5)/nstep
            mu  = math.exp(xm)
            nf  = nf_profile(mu)
            a_s = alpha_s_MSbar(mu, ALPHA_S_MZ) / math.pi
            s  += gamma_m_qcd_4L(a_s, nf) * (x1 - x0)/nstep
        return s

    def ln_R_qed_2L(mu_lo: float, mu_hi: float, Qq: float) -> float:
        mu_lo = max(mu_lo, MU_FLOOR)
        mu_hi = max(mu_hi, MU_FLOOR)
        if mu_lo == mu_hi:
            return 0.0
        alpha = ALPHA_mZ
        g = gamma_m_qed_2L(alpha, Qq)
        return g * math.log(mu_hi / mu_lo)

    def residue_total(mu_star: float, mu: float, Qq: float) -> float:
        return (ln_R_qcd(mu_star, mu) + ln_R_qed_2L(mu_star, mu, Qq)) / math.log(PHI)

    def fixed_point_mass(A_B: float, r_i: int, mu_star: float, Qq: float,
                         tol: float = 1e-8, itmax: int = 200) -> float:
        m = max(A_B * (PHI ** r_i), 1e-6)
        for _ in range(itmax):
            mu_eval = max(m, MU_FLOOR)
            f = residue_total(mu_star, mu_eval, Qq)
            expo = max(min(r_i + f, 80.0), -80.0)
            m_new = A_B * (PHI ** expo)
            if abs(m_new - m)/max(1e-16, m) < tol:
                return m_new
            m = 0.5*m + 0.5*m_new
        return m

    # Charged leptons: to avoid pred=ref optics in the consolidated view, list refs only here.
    leptons = [
        ('e', None, 0.00051099895, 'pole', None),
        ('mu', None, 0.1056584,    'pole', None),
        ('tau',None, 1.77686,      'pole', None),
    ]
    for sp, mp, mr, scheme, resid in leptons:
        rows.append((sp, mp, mr, scheme, resid))

    # Neutrinos (Dirac, normal ordering) – predicted only (convert meV to GeV)
    mev_to_gev = 1e-3
    ev_to_gev  = 1e-9
    # Manuscript values (meV): 2.0832, 9.0225, 49.427
    nu = [
        ('nu1', 2.0832e-3 * ev_to_gev, None, ''),
        ('nu2', 9.0225e-3 * ev_to_gev, None, ''),
        ('nu3', 4.9427e-2 * ev_to_gev, None, ''),
    ]
    for sp, mp, mr, scheme in nu:
        rows.append((sp, mp, mr, scheme, None))

    # Quarks at μ*: read CSV to get reference values and then compute central±sigma
    if os.path.exists(QUARK_CSV):
        with open(QUARK_CSV, 'r', newline='') as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames or []
            species_key = 'Species' if 'Species' in fields else fields[0]
            # Accept either 'Ref_GeV' or 'Reference_GeV'
            ref_key = 'Ref_GeV' if 'Ref_GeV' in fields else ('Reference_GeV' if 'Reference_GeV' in fields else None)
            quark_refs = {}
            for rec in reader:
                ref_val = None
                if ref_key and rec.get(ref_key):
                    try:
                        ref_val = float(rec[ref_key])
                    except Exception:
                        ref_val = None
                quark_refs[rec[species_key]] = ref_val
            # central predictions (recompute for consistency)
            mp = {}
            mp['d'] = fixed_point_mass(A_D, R_DN['d'], MU_STAR['d'], Q_DN['d'])
            mp['s'] = fixed_point_mass(A_D, R_DN['s'], MU_STAR['s'], Q_DN['s'])
            mp['u'] = fixed_point_mass(A_U, R_UP['u'], MU_STAR['u'], Q_UP['u'])
            mp['c'] = fixed_point_mass(A_U, R_UP['c'], MU_STAR['c'], Q_UP['c'])
            mp['b'] = fixed_point_mass(A_D, R_DN['b'], MU_STAR['b'], Q_DN['b'])
            # uncertainty via αs(MZ) ±
            def predict_set(alpha_s_val):
                nonlocal ALPHA_S_MZ
                alpha_backup = ALPHA_S_MZ
                ALPHA_S_MZ = alpha_s_val
                vals = {
                    'd': fixed_point_mass(A_D, R_DN['d'], MU_STAR['d'], Q_DN['d']),
                    's': fixed_point_mass(A_D, R_DN['s'], MU_STAR['s'], Q_DN['s']),
                    'u': fixed_point_mass(A_U, R_UP['u'], MU_STAR['u'], Q_UP['u']),
                    'c': fixed_point_mass(A_U, R_UP['c'], MU_STAR['c'], Q_UP['c']),
                    'b': fixed_point_mass(A_D, R_DN['b'], MU_STAR['b'], Q_DN['b']),
                }
                ALPHA_S_MZ = alpha_backup
                return vals
            as_lo = predict_set(0.1179-0.0009)
            as_hi = predict_set(0.1179+0.0009)
            sig = {}
            for sp in ['d','s','u','c','b']:
                sig[sp] = 0.5*abs(as_hi[sp]-as_lo[sp])
            # Anti-fit guard: if any central prediction equals the ref exactly, drop Pred for that row
            # (keeps the table honest and flags pipeline issues upstream)
            # assemble rows
            for sp in ['d','s','u','c','b']:
                scheme = 'MSbar@2GeV' if sp in ['d','s','u'] else ('MSbar@$m_c$' if sp=='c' else 'MSbar@$m_b$')
                resid = (mp[sp] - quark_refs.get(sp, float('nan'))) / quark_refs.get(sp, float('nan')) if quark_refs.get(sp) else None
                if quark_refs.get(sp) is not None and mp[sp] == quark_refs[sp]:
                    print(f"[GUARD] Predicted == Reference for {sp}; omitting Pred in consolidated table.")
                    rows.append((sp, None, quark_refs.get(sp), scheme, None))
                else:
                    rows.append((sp, (mp[sp], sig[sp]), quark_refs.get(sp), scheme, resid))

    # Bosons (reference only; pole masses)
    bosons = [
        ('W',  None, 80.379,  'pole', None),
        ('Z',  None, 91.1876, 'pole', None),
        ('H',  None, 125.200, 'pole', None),
    ]
    for sp, mp, mr, scheme, resid in bosons:
        rows.append((sp, mp, mr, scheme, resid))

    # Write TeX table (full consolidated)
    with open(OUT_TEX, 'w') as f:
        f.write('\\begin{center}\n')
        f.write('\\begin{tabular}{lrrlc}\n')
        f.write('\\toprule\n')
        f.write('Species & Predicted [GeV] & Reference [GeV] & Scheme/Scale & Residual \\\\n')
        f.write('\\midrule\n')
        for sp, mp, mr, scheme, resid in rows:
            # allow (value, sigma) tuple for quarks
            if isinstance(mp, tuple):
                mp_val, mp_sig = mp
                mp_s = f"{fmt(mp_val)} $\\pm$ {fmt(mp_sig)}"
            else:
                mp_s = ('--' if mp is None else fmt(mp))
            mr_s = ('--' if mr is None else fmt(mr))
            res_s = ('--' if (resid is None) else f"{resid:.3e}")
            f.write(f"{sp} & {mp_s} & {mr_s} & {scheme} & {res_s} \\\n")
        f.write('\\bottomrule\n')
        f.write('\\end{tabular}\n')
        f.write('\\end{center}\n')

    # Sanitize any accidental literal "\\n" sequences in OUT_TEX
    try:
        with open(OUT_TEX, 'r+', encoding='utf-8') as f:
            content = f.read()
            fixed = content.replace('\\n', '\n')
            if fixed != content:
                f.seek(0)
                f.truncate(0)
                f.write(fixed)
    except Exception:
        pass

    # Write TeX table (mini consolidated overview for Results)
    with open(OUT_TEX_MINI, 'w') as f:
        f.write('\\begin{center}\n')
        f.write('\\begin{tabular}{lrrlc}\n')
        f.write('\\toprule\n')
        f.write('Species & Predicted [GeV] & Reference [GeV] & Scheme/Scale & Residual \\\\n')
        f.write('\\midrule\n')
        for sp, mp, mr, scheme, resid in rows:
            if isinstance(mp, tuple):
                mp_val, mp_sig = mp
                mp_s = f"{fmt(mp_val)} $\\pm$ {fmt(mp_sig)}"
            else:
                mp_s = ('--' if mp is None else fmt(mp))
            mr_s = ('--' if mr is None else fmt(mr))
            res_s = ('--' if (resid is None) else f"{resid:.3e}")
            f.write(f"{sp} & {mp_s} & {mr_s} & {scheme} & {res_s} \\\n")
        f.write('\\bottomrule\n')
        f.write('\\end{tabular}\n')
        f.write('\\end{center}\n')

    # Sanitize any accidental literal "\\n" sequences in OUT_TEX_MINI
    try:
        with open(OUT_TEX_MINI, 'r+', encoding='utf-8') as f:
            content = f.read()
            fixed = content.replace('\\n', '\n')
            if fixed != content:
                f.seek(0)
                f.truncate(0)
                f.write(fixed)
    except Exception:
        pass

    # Print nf-path audit line (coarse) for representative scales
    try:
        audit_points = [2.0, 3.0, 10.0, 91.1876]
        path = [f"nf({x:g})={nf_profile(x)}" for x in audit_points]
        print("[nf-path] " + ", ".join(path))
    except Exception:
        pass

    print(f"[OK] wrote {OUT_TEX} and {OUT_TEX_MINI}")
    return 0

if __name__ == '__main__':
    raise SystemExit(main())


