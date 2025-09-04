#!/usr/bin/env python3
"""
Minimal quark RG runner scaffold exposing residue_total for QCD(4L)+QED(2L).

This file intentionally avoids hard-coding 4L coefficients and instead expects
either a mass_run_factor(mu_lo, mu_hi) or run_mass_MSbar(m0, mu0, mu1) that
already incorporate 4L QCD running with decoupling (e.g., via CRunDec).

If those are not yet wired, the functions will not crash table generation; the
writer script will skip emission until predictions succeed.
"""

from dataclasses import dataclass
import math

PHI: float = (1.0 + 5.0 ** 0.5) / 2.0
LNPHI: float = math.log(PHI)


@dataclass
class PDGInputs:
    # Light MSbar masses at 2 GeV (GeV)
    m_u_2: float = 0.00216
    m_d_2: float = 0.00467
    m_s_2: float = 0.093
    # Heavy MSbar reference scales (GeV)
    m_c_mc: float = 1.27
    m_b_mb: float = 4.18
    # Top pole (not used by default table)
    m_t_pole: float = 172.69


class QuarkRunner:
    def __init__(self, pdg: PDGInputs):
        self.pdg = pdg

    # --- QCD mass ratio ln R_QCD ---
    def _ln_R_qcd(self, mu_lo: float, mu_hi: float) -> float:
        # Prefer subclass to provide precise 4L mass_run_factor; fallback is neutral
        if hasattr(self, "mass_run_factor"):
            R = self.mass_run_factor(mu_lo, mu_hi)  # expected to include decoupling
            return math.log(max(R, 1e-300))
        if hasattr(self, "run_mass_MSbar"):
            m_seed = 1.0
            m_hi = self.run_mass_MSbar(m_seed, mu_lo, mu_hi)
            return math.log(max(m_hi / m_seed, 1e-300))
        # No backend: return 0 so residue is purely QED (effectively skip QCD)
        return 0.0

    # --- QED(2L) mass ratio ln R_QED ---
    def _ln_R_qed_2L(self, mu_lo: float, mu_hi: float, Qq: float) -> float:
        # Use constant alpha at mZ if no running alpha provided
        alpha = self.alpha_em(91.1876) if hasattr(self, "alpha_em") else 1.0 / 127.955
        a = alpha / math.pi
        g1 = -3.0 * (Qq ** 2) * a
        g2 = -1.5 * (Qq ** 4) * a * a
        return (g1 + g2) * math.log(mu_hi / mu_lo)

    # --- Public residue ---
    def residue_total(self, mu_star: float, mu: float, Qq: float) -> float:
        lnR_qcd = self._ln_R_qcd(mu_star, mu)
        lnR_qed = self._ln_R_qed_2L(mu_star, mu, Qq)
        return (lnR_qcd + lnR_qed) / LNPHI

#!/usr/bin/env python3
# Minimal quark RG utilities using CRunDec for MSbar running with threshold stepping.

from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple
import math

from rundec import CRunDec


@dataclass
class PDGInputs:
    alpha_s_MZ: float = 0.1179
    MZ: float = 91.1876
    sigma_alpha_s: float = 0.0009
    # MSbar reference masses (GeV)
    m_u_2: float = 0.00216  # at 2 GeV
    m_d_2: float = 0.00467  # at 2 GeV
    m_s_2: float = 0.0934   # at 2 GeV
    m_c_mc: float = 1.27    # at m_c
    m_b_mb: float = 4.18    # at m_b
    # 1-sigma uncertainties (GeV) — conservative PDG-style
    sigma_m_u_2: float = 0.00011
    sigma_m_d_2: float = 0.00020
    sigma_m_s_2: float = 0.00110
    sigma_m_c_mc: float = 0.02
    sigma_m_b_mb: float = 0.03
    # Thresholds (GeV). Use MSbar values for c,b; top not needed below MZ.
    mc_thr: float = 1.27
    mb_thr: float = 4.18


class QuarkRunner:
    def __init__(self, pdg: PDGInputs, loops: int = 4):
        self.pdg = pdg
        self.loops = loops
        self.rd = CRunDec()
        # Purely classical pipeline: no RS sector constants or rungs.
        # Compute Lambda_5 from alpha_s(MZ)
        # signature LamImpl(asmu, mu, nf, nloops)
        try:
            self.Lam5 = self.rd.LamImpl(pdg.alpha_s_MZ, pdg.MZ, 5, loops)
        except Exception:
            # fallback to 2-loop for robustness
            self.loops = 2
            self.Lam5 = self.rd.LamImpl(pdg.alpha_s_MZ, pdg.MZ, 5, 2)
        # Derive Lambda_4 and Lambda_3 by decoupling down at thresholds
        # Use DecLambdaDown(alh, Mth, nth) with appropriate inputs
        try:
            self.Lam4 = self.rd.DecLambdaDown(self.Lam5, pdg.mb_thr, 5, self.loops)
        except Exception:
            self.Lam4 = self.rd.DecLambdaDown(self.Lam5, pdg.mb_thr, 5, 2)
        try:
            self.Lam3 = self.rd.DecLambdaDown(self.Lam4, pdg.mc_thr, 4, self.loops)
        except Exception:
            self.Lam3 = self.rd.DecLambdaDown(self.Lam4, pdg.mc_thr, 4, 2)

    # --- Minimal helpers for EM running and nf profile (classical stubs) ---
    def alpha_em(self, mu: float) -> float:
        """Return alpha(μ). Stub: use PDG alpha(MZ) as constant until VP model is wired."""
        return 1.0 / 127.955

    def nf_profile(self, mu: float) -> int:
        """Active flavours for QCD running based on MSbar thresholds."""
        if mu < self.pdg.mc_thr:
            return 3
        if mu < self.pdg.mb_thr:
            return 4
        # below MZ we rarely need nf=6; cap at 5 for this pipeline
        return 5

    def alphas(self, mu: float, nf: int) -> float:
        Lam = {3: self.Lam3, 4: self.Lam4, 5: self.Lam5}[nf]
        try:
            return self.rd.AlphasLam(Lam, mu, nf, self.loops)
        except Exception:
            return self.rd.AlphasLam(Lam, mu, nf)

    def _run_factor_segment(self, m0: float, mu0: float, mu1: float, nf: int) -> float:
        """Return m(mu1) given m(mu0)=m0 in fixed nf using RGI transport."""
        if mu0 == mu1:
            return m0
        return self.run_mass_between(m0, mu0, mu1, nf)

    def mass_run_factor(self, mu_lo: float, mu_hi: float) -> float:
        """Compute factor m(mu_hi)/m(mu_lo), stepping thresholds at mc, mb.
        Uses nf=3 below mc, nf=4 between mc..mb, nf=5 above mb.
        """
        if mu_lo <= 0.0 or mu_hi <= 0.0:
            return 1.0
        if abs(mu_hi - mu_lo) < 1e-16:
            return 1.0
        # impose a safe infrared floor to avoid alpha_s blow-up in nf=3
        mu_floor = 2.0
        if mu_lo < mu_floor:
            mu_lo = mu_floor
        m = 1.0
        mb = self.pdg.mb_thr
        mc = self.pdg.mc_thr
        if mu_lo < mu_hi:
            # upward
            cur_mu = mu_lo
            # up to mc
            if cur_mu < mc:
                to_mu = min(mc, mu_hi)
                m = self._run_factor_segment(m, cur_mu, to_mu, 3)
                cur_mu = to_mu
                if mu_hi > mc:
                    # decouple 3->4 at mc
                    m = self.decouple_up_mass(m, mc, 3)
            # up to mb
            if cur_mu < mb and mu_hi > cur_mu:
                to_mu = min(mb, mu_hi)
                m = self._run_factor_segment(m, cur_mu, to_mu, 4)
                if mu_hi > mb:
                    m = self.decouple_up_mass(m, mb, 4)
                cur_mu = to_mu
            # above mb
            if mu_hi > cur_mu:
                m = self._run_factor_segment(m, cur_mu, mu_hi, 5)
        else:
            # downward
            cur_mu = mu_lo
            # from >mb down to <=mb
            if cur_mu > mb:
                to_mu = max(mb, mu_hi)
                m = self._run_factor_segment(m, cur_mu, to_mu, 5)
                if mu_hi < mb:
                    m = self.decouple_down_mass(m, mb, 5)
                cur_mu = to_mu
            # from between mc..mb down
            if cur_mu > mc and mu_hi < cur_mu:
                to_mu = max(mc, mu_hi)
                m = self._run_factor_segment(m, cur_mu, to_mu, 4)
                if mu_hi < mc:
                    m = self.decouple_down_mass(m, mc, 4)
                cur_mu = to_mu
            # below mc
            if mu_hi < cur_mu:
                m = self._run_factor_segment(m, cur_mu, mu_hi, 3)
        return m

    # No RS fixed-point mass prediction for quarks in the classical pipeline.

    # (removed) — no RS-derived quark "predicted masses" here

    def predicted_with_uncertainty(self, flavors: list[str], sigma_alpha_s: float = 0.0009) -> dict:
        """Return central and ±1σ from varying alpha_s(MZ) globally.
        """
        # central
        cen = {f: self.fixed_point_mass(f) for f in flavors}
        # minus
        pdg_minus = PDGInputs(alpha_s_MZ=self.pdg.alpha_s_MZ - sigma_alpha_s,
                              MZ=self.pdg.MZ,
                              m_u_2=self.pdg.m_u_2, m_d_2=self.pdg.m_d_2, m_s_2=self.pdg.m_s_2,
                              m_c_mc=self.pdg.m_c_mc, m_b_mb=self.pdg.m_b_mb,
                              mc_thr=self.pdg.mc_thr, mb_thr=self.pdg.mb_thr)
        qm_minus = QuarkRunner(pdg_minus, loops=self.loops)
        lo = {f: qm_minus.fixed_point_mass(f) for f in flavors}
        # plus
        pdg_plus = PDGInputs(alpha_s_MZ=self.pdg.alpha_s_MZ + sigma_alpha_s,
                             MZ=self.pdg.MZ,
                             m_u_2=self.pdg.m_u_2, m_d_2=self.pdg.m_d_2, m_s_2=self.pdg.m_s_2,
                             m_c_mc=self.pdg.m_c_mc, m_b_mb=self.pdg.m_b_mb,
                             mc_thr=self.pdg.mc_thr, mb_thr=self.pdg.mb_thr)
        qm_plus = QuarkRunner(pdg_plus, loops=self.loops)
        hi = {f: qm_plus.fixed_point_mass(f) for f in flavors}
        out = {}
        for f in flavors:
            sigma = 0.5 * abs(hi[f] - lo[f])
            out[f] = (cen[f], sigma)
        return out

    def run_mass_between(self, m_ref: float, mu_ref: float, mu_tgt: float, nf: int) -> float:
        # Use RGI mass to transport between scales in fixed nf
        as_ref = self.alphas(mu_ref, nf)
        as_tgt = self.alphas(mu_tgt, nf)
        m_rgi = self.rd.mMS2mRGI(m_ref, as_ref, nf)
        m_tgt = self.rd.mRGI2mMS(m_rgi, as_tgt, nf)
        return m_tgt

    def decouple_up_mass(self, m_light: float, Mth: float, nf_from: int) -> float:
        # DecMqUpMS(mL, asL, asH, Mth, nfL)
        # Build alpha_s at threshold on both sides
        # get Lambda for nf_from and nf_from+1
        if nf_from == 3:
            LamL = self.Lam3
            LamH = self.Lam4
        elif nf_from == 4:
            LamL = self.Lam4
            LamH = self.Lam5
        else:
            raise ValueError('nf_from must be 3 or 4')
        try:
            asL = self.rd.AlphasLam(LamL, Mth, nf_from, self.loops)
            asH = self.rd.AlphasLam(LamH, Mth, nf_from+1, self.loops)
            return self.rd.DecMqUpMS(m_light, asL, asH, Mth, nf_from, self.loops)
        except Exception:
            asL = self.rd.AlphasLam(LamL, Mth, nf_from, 2)
            asH = self.rd.AlphasLam(LamH, Mth, nf_from+1, 2)
            return self.rd.DecMqUpMS(m_light, asL, asH, Mth, nf_from)

    def decouple_down_mass(self, m_light: float, Mth: float, nf_from: int) -> float:
        as_thr = self.alphas(Mth, nf_from)
        try:
            return self.rd.DecMqDownMS(m_light, as_thr, Mth, nf_from, self.loops)
        except Exception:
            return self.rd.DecMqDownMS(m_light, as_thr, Mth, nf_from)

    def run_to_MZ_light(self, m2: float) -> float:
        # Run mass defined at 2 GeV (nf=3) to MZ, stepping nf=3->4 at mc, 4->5 at mb
        m = m2
        # 3 -> run to mc
        m = self.run_mass_between(m, 2.0, self.pdg.mc_thr, 3)
        # decouple up at mc: 3->4
        m = self.decouple_up_mass(m, self.pdg.mc_thr, 3)
        # 4 -> run to mb
        m = self.run_mass_between(m, self.pdg.mc_thr, self.pdg.mb_thr, 4)
        # decouple up at mb: 4->5
        m = self.decouple_up_mass(m, self.pdg.mb_thr, 4)
        # 5 -> run to MZ
        m = self.run_mass_between(m, self.pdg.mb_thr, self.pdg.MZ, 5)
        return m

    def run_to_MZ_heavy(self, mQ_mQ: float, mQ: float) -> float:
        # treat in nf=5 up to MZ (adequate for c/b to MZ)
        return self.run_mass_between(mQ_mQ, mQ, self.pdg.MZ, 5)

    def compute_MZ_masses(self) -> Tuple[float, float, float, float, float]:
        mu_MZ = self.run_to_MZ_light(self.pdg.m_u_2)
        md_MZ = self.run_to_MZ_light(self.pdg.m_d_2)
        ms_MZ = self.run_to_MZ_light(self.pdg.m_s_2)
        mc_MZ = self.run_to_MZ_heavy(self.pdg.m_c_mc, self.pdg.m_c_mc)
        mb_MZ = self.run_to_MZ_heavy(self.pdg.m_b_mb, self.pdg.m_b_mb)
        return mu_MZ, md_MZ, ms_MZ, mc_MZ, mb_MZ

    def compute_mu_star(self, name: str) -> float:
        """Self-consistent μ* such that μ* = mbar_i(μ*) in MSbar (classical convention).
        Light flavours use the standard 2 GeV convention (nf=3) for transparency."""
        if name == 'u':
            return 2.0
            # Use nf stepping and light mass at 2 GeV
            mref = self.pdg.m_u_2
            refmu = 2.0
            nf_at = lambda x: 3 if x < self.pdg.mc_thr else (4 if x < self.pdg.mb_thr else 5)
        elif name == 'd':
            return 2.0
            mref = self.pdg.m_d_2
            refmu = 2.0
            nf_at = lambda x: 3 if x < self.pdg.mc_thr else (4 if x < self.pdg.mb_thr else 5)
        elif name == 's':
            return 2.0
            mref = self.pdg.m_s_2
            refmu = 2.0
            nf_at = lambda x: 3 if x < self.pdg.mc_thr else (4 if x < self.pdg.mb_thr else 5)
        elif name == 'c':
            mu = self.pdg.m_c_mc
            mref = self.pdg.m_c_mc
            refmu = self.pdg.m_c_mc
            nf_at = lambda x: 4 if x < self.pdg.mb_thr else 5
        elif name == 'b':
            mu = self.pdg.m_b_mb
            mref = self.pdg.m_b_mb
            refmu = self.pdg.m_b_mb
            nf_at = lambda x: 5
        else:
            raise ValueError(name)

        for _ in range(100):
            nf = nf_at(mu)
            # run m from reference to μ (fixed-nf segments with decoupling at thresholds)
            if refmu == 2.0 and mu <= self.pdg.mc_thr:
                mu_eff = max(mu, 1.0)
                m_mu = self.run_mass_between(mref, refmu, mu_eff, 3)
            else:
                m_temp = mref
                cur_mu = refmu
                cur_nf = 3 if refmu <= self.pdg.mc_thr else (4 if refmu <= self.pdg.mb_thr else 5)
                if refmu == 2.0:
                    m_temp = self.run_mass_between(m_temp, cur_mu, self.pdg.mc_thr, 3)
                    m_temp = self.decouple_up_mass(m_temp, self.pdg.mc_thr, 3)
                    cur_mu = self.pdg.mc_thr
                    cur_nf = 4
                if mu > self.pdg.mb_thr:
                    if cur_nf == 4:
                        m_temp = self.run_mass_between(m_temp, cur_mu, self.pdg.mb_thr, 4)
                        m_temp = self.decouple_up_mass(m_temp, self.pdg.mb_thr, 4)
                        cur_mu = self.pdg.mb_thr
                        cur_nf = 5
                    m_temp = self.run_mass_between(m_temp, cur_mu, mu, 5)
                else:
                    m_temp = self.run_mass_between(m_temp, cur_mu, mu, 4)
                m_mu = m_temp
            new_mu = m_mu
            if abs(new_mu - mu) / max(1e-9, mu) < 1e-10:
                mu = new_mu
                break
            mu = new_mu
        return mu

    # --- Global residue building blocks (stubs compile; fill coefficients later) ---
    def gamma_m_qcd_4L(self, alphas: float, nf: int) -> float:
        """MSbar quark mass anomalous dimension through 4 loops.
        Convention: a_s = α_s/π,  γ_m(a_s) = - (γ0 a + γ1 a^2 + γ2 a^3 + γ3 a^4)
        Coefficients specialized to SU(3) with nf active flavours.
        """
        ZETA3 = 1.2020569031595942
        a = alphas / math.pi
        # γ0, γ1
        g0 = 1.0
        g1 = (101.0/24.0) - (5.0/36.0)*nf
        # γ2
        g2 = (1249.0/64.0) - (277.0/216.0)*nf - (35.0/1296.0)*nf*nf
        # γ3 with ζ3
        g3 = (4603055.0/41472.0 + (41.0/18.0)*ZETA3) \
             - ((91723.0/20736.0) + (11.0/12.0)*ZETA3)*nf \
             + (151.0/162.0)*nf*nf + (5.0/972.0)*nf*nf*nf
        return -(g0*a + g1*a*a + g2*a*a*a + g3*a*a*a*a)

    def gamma_m_qed_2L(self, alpha: float, Qq: float) -> float:
        """QED mass anomalous dimension up to 2 loops; returns negative by convention."""
        a = alpha / math.pi
        return -(3.0 * (Qq ** 2) * a + 1.5 * (Qq ** 4) * a * a)

    def beta_qed_2L(self, alpha: float, nf_em: int) -> float:
        """QED beta up to 2 loops for running alpha(μ) if desired.
        Not used if alpha(μ) is supplied externally via dispersion.
        """
        a = alpha / math.pi
        b1 = (2.0/3.0) * nf_em
        b2 = 0.25 * nf_em
        return b1*a*a + b2*a*a*a

    def residue_integral(self, mu_star: float, mu: float, Qq: float, nf_profile=None) -> float:
        """Compute f_q(μ) = (1/ln φ) ∫_{ln μ*}^{ln μ} [γ_m^QCD + γ_m^QED] d ln μ'.
        Integrates on a uniform log grid; uses self.alphas and self.alpha_em with fixed decoupling.
        """
        lnphi = math.log((1.0 + math.sqrt(5.0)) / 2.0)
        # handle degenerate case
        if mu == mu_star:
            return 0.0
        if nf_profile is None:
            nf_profile = self.nf_profile
        # integrate with simple midpoint rule
        nstep = 1000
        log_lo, log_hi = math.log(min(mu_star, mu)), math.log(max(mu_star, mu))
        s = 0.0
        for k in range(nstep):
            x0 = log_lo + (log_hi - log_lo) * k / nstep
            x1 = log_lo + (log_hi - log_lo) * (k + 1) / nstep
            xm = 0.5 * (x0 + x1)
            mum = math.exp(xm)
            nf = nf_profile(mum)
            as_m = self.alphas(mum, nf)
            a_em = self.alpha_em(mum)
            gm = self.gamma_m_qcd_4L(as_m, nf) + self.gamma_m_qed_2L(a_em, Qq)
            s += gm * (x1 - x0)
        # direction (mu above/below mu_star)
        if mu < mu_star:
            s = -s
        return s / lnphi

    # Convenience: PDG-driven ratios at mu* (light flavours use 2 GeV convention)
    def pdg_ratios_mu_star(self) -> dict:
        mu_u = 2.0; mu_d = 2.0; mu_s = 2.0
        mu_c = self.pdg.m_c_mc; mu_b = self.pdg.m_b_mb
        # Masses at their conventional reference points
        m_u = self.pdg.m_u_2
        m_d = self.pdg.m_d_2
        m_s = self.pdg.m_s_2
        m_c = self.pdg.m_c_mc
        m_b = self.pdg.m_b_mb
        return {
            's/d': m_s/m_d,
            'b/s': m_b/m_s,
            'b/d': m_b/m_d,
            'c/u': m_c/m_u,
        }

    def pdg_ratios_mu_star_with_unc(self) -> dict:
        # central
        cen = self.pdg_ratios_mu_star()
        # fractional priors
        fu = self.pdg.sigma_m_u_2 / max(1e-16, self.pdg.m_u_2)
        fd = self.pdg.sigma_m_d_2 / max(1e-16, self.pdg.m_d_2)
        fs = self.pdg.sigma_m_s_2 / max(1e-16, self.pdg.m_s_2)
        fc = self.pdg.sigma_m_c_mc / max(1e-16, self.pdg.m_c_mc)
        fb = self.pdg.sigma_m_b_mb / max(1e-16, self.pdg.m_b_mb)
        # uncorrelated propagation for ratios
        sig = {}
        sig['s/d'] = cen['s/d'] * ( (fs*fs + fd*fd) ** 0.5 )
        sig['b/s'] = cen['b/s'] * ( (fb*fb + fs*fs) ** 0.5 )
        sig['b/d'] = cen['b/d'] * ( (fb*fb + fd*fd) ** 0.5 )
        sig['c/u'] = cen['c/u'] * ( (fc*fc + fu*fu) ** 0.5 )
        return {k: (cen[k], sig[k]) for k in cen}

    def ratios_MZ_with_unc(self) -> dict:
        # central ratios at MZ
        mu_MZ, md_MZ, ms_MZ, mc_MZ, mb_MZ = self.compute_MZ_masses()
        cen = {
            's/d': ms_MZ / md_MZ,
            'b/s': mb_MZ / ms_MZ,
            'c/u': mc_MZ / mu_MZ,
        }
        # alpha_s contribution via +/- sigma
        pdg_lo = PDGInputs(alpha_s_MZ=self.pdg.alpha_s_MZ - self.pdg.sigma_alpha_s,
                           MZ=self.pdg.MZ,
                           m_u_2=self.pdg.m_u_2, m_d_2=self.pdg.m_d_2, m_s_2=self.pdg.m_s_2,
                           m_c_mc=self.pdg.m_c_mc, m_b_mb=self.pdg.m_b_mb)
        pdg_hi = PDGInputs(alpha_s_MZ=self.pdg.alpha_s_MZ + self.pdg.sigma_alpha_s,
                           MZ=self.pdg.MZ,
                           m_u_2=self.pdg.m_u_2, m_d_2=self.pdg.m_d_2, m_s_2=self.pdg.m_s_2,
                           m_c_mc=self.pdg.m_c_mc, m_b_mb=self.pdg.m_b_mb)
        lo = QuarkRunner(pdg_lo, loops=self.loops)
        hi = QuarkRunner(pdg_hi, loops=self.loops)
        mu_MZ_l, md_MZ_l, ms_MZ_l, mc_MZ_l, mb_MZ_l = lo.compute_MZ_masses()
        mu_MZ_h, md_MZ_h, ms_MZ_h, mc_MZ_h, mb_MZ_h = hi.compute_MZ_masses()
        r_lo = {
            's/d': ms_MZ_l / md_MZ_l,
            'b/s': mb_MZ_l / ms_MZ_l,
            'c/u': mc_MZ_l / mu_MZ_l,
        }
        r_hi = {
            's/d': ms_MZ_h / md_MZ_h,
            'b/s': mb_MZ_h / ms_MZ_h,
            'c/u': mc_MZ_h / mu_MZ_h,
        }
        sig_as = {k: 0.5 * abs(r_hi[k] - r_lo[k]) for k in cen}
        # prior mass contributions: same fractional priors as at reference scales, largely cancel when evolved; we include them conservatively
        fu = self.pdg.sigma_m_u_2 / max(1e-16, self.pdg.m_u_2)
        fd = self.pdg.sigma_m_d_2 / max(1e-16, self.pdg.m_d_2)
        fs = self.pdg.sigma_m_s_2 / max(1e-16, self.pdg.m_s_2)
        fc = self.pdg.sigma_m_c_mc / max(1e-16, self.pdg.m_c_mc)
        fb = self.pdg.sigma_m_b_mb / max(1e-16, self.pdg.m_b_mb)
        sig_prior = {
            's/d': cen['s/d'] * ( (fs*fs + fd*fd) ** 0.5 ),
            'b/s': cen['b/s'] * ( (fb*fb + fs*fs) ** 0.5 ),
            'c/u': cen['c/u'] * ( (fc*fc + fu*fu) ** 0.5 ),
        }
        sig_tot = {k: (sig_as[k]**2 + sig_prior[k]**2) ** 0.5 for k in cen}
        return {k: (cen[k], sig_tot[k]) for k in cen}

    # --- Predictions at mu* (skeleton). Return empty until kernels are filled. ---
    def _predicted_ratios_mu_star_central(self) -> dict:
        """Central predicted ratios at μ* (no uncertainties)."""
        # Frozen sector constants and integers (parameter-free layer)
        PHI = (1.0 + math.sqrt(5.0)) / 2.0
        ecoh_eV = PHI ** (-5.0)
        ecoh_GeV = ecoh_eV * 1e-9
        A_U = (2.0 ** (-1)) * ecoh_GeV * (PHI ** 35)
        A_D = (2.0 ** (23)) * ecoh_GeV * (PHI ** (-5))
        r_up = {'u': 4, 'c': 15, 't': 21}
        r_dn = {'d': 4, 's': 15, 'b': 21}
        # Scales and charges
        mu_star = {'u': 2.0, 'd': 2.0, 's': 2.0, 'c': self.pdg.m_c_mc, 'b': self.pdg.m_b_mb}
        Q_up = {'u':  2.0/3.0, 'c': 2.0/3.0, 't': 2.0/3.0}
        Q_dn = {'d': -1.0/3.0, 's':-1.0/3.0, 'b':-1.0/3.0}
        try:
            from rs_predictor import RSPredictor
        except Exception:
            return {}
        pred_u = RSPredictor(A_U, r_up, self)
        pred_d = RSPredictor(A_D, r_dn, self)
        # Central fixed points
        m_u = pred_u.fixed_point('u', mu_star['u'], Q_up['u'])
        m_d = pred_d.fixed_point('d', mu_star['d'], Q_dn['d'])
        m_s = pred_d.fixed_point('s', mu_star['s'], Q_dn['s'])
        m_c = pred_u.fixed_point('c', mu_star['c'], Q_up['c'])
        m_b = pred_d.fixed_point('b', mu_star['b'], Q_dn['b'])
        return {
            's/d': (m_s/m_d, float('nan')),
            'b/s': (m_b/m_s, float('nan')),
            'b/d': (m_b/m_d, float('nan')),
            'c/u': (m_c/m_u, float('nan')),
        }

    def predicted_ratios_mu_star_with_unc(self) -> dict:
        """Return predicted ratios with 1σ at μ* using global variations (α_s, thresholds)."""
        cen = self._predicted_ratios_mu_star_central()
        if not cen:
            return {}
        # α_s variations
        pdg_lo = PDGInputs(alpha_s_MZ=self.pdg.alpha_s_MZ - self.pdg.sigma_alpha_s,
                           MZ=self.pdg.MZ, m_u_2=self.pdg.m_u_2, m_d_2=self.pdg.m_d_2, m_s_2=self.pdg.m_s_2,
                           m_c_mc=self.pdg.m_c_mc, m_b_mb=self.pdg.m_b_mb, mc_thr=self.pdg.mc_thr, mb_thr=self.pdg.mb_thr)
        pdg_hi = PDGInputs(alpha_s_MZ=self.pdg.alpha_s_MZ + self.pdg.sigma_alpha_s,
                           MZ=self.pdg.MZ, m_u_2=self.pdg.m_u_2, m_d_2=self.pdg.m_d_2, m_s_2=self.pdg.m_s_2,
                           m_c_mc=self.pdg.m_c_mc, m_b_mb=self.pdg.m_b_mb, mc_thr=self.pdg.mc_thr, mb_thr=self.pdg.mb_thr)
        lo = QuarkRunner(pdg_lo, loops=self.loops)
        hi = QuarkRunner(pdg_hi, loops=self.loops)
        r_lo = lo._predicted_ratios_mu_star_central()
        r_hi = hi._predicted_ratios_mu_star_central()
        # Threshold policy variations (±1%)
        pdg_thr_dn = PDGInputs(alpha_s_MZ=self.pdg.alpha_s_MZ, MZ=self.pdg.MZ,
                               m_u_2=self.pdg.m_u_2, m_d_2=self.pdg.m_d_2, m_s_2=self.pdg.m_s_2,
                               m_c_mc=self.pdg.m_c_mc, m_b_mb=self.pdg.m_b_mb,
                               mc_thr=self.pdg.mc_thr*0.99, mb_thr=self.pdg.mb_thr*0.99)
        pdg_thr_up = PDGInputs(alpha_s_MZ=self.pdg.alpha_s_MZ, MZ=self.pdg.MZ,
                               m_u_2=self.pdg.m_u_2, m_d_2=self.pdg.m_d_2, m_s_2=self.pdg.m_s_2,
                               m_c_mc=self.pdg.m_c_mc, m_b_mb=self.pdg.m_b_mb,
                               mc_thr=self.pdg.mc_thr*1.01, mb_thr=self.pdg.mb_thr*1.01)
        thr_dn = QuarkRunner(pdg_thr_dn, loops=self.loops)
        thr_up = QuarkRunner(pdg_thr_up, loops=self.loops)
        r_thr_dn = thr_dn._predicted_ratios_mu_star_central()
        r_thr_up = thr_up._predicted_ratios_mu_star_central()
        out = {}
        for key in ['s/d','b/s','b/d','c/u']:
            rc = cen[key][0]
            # α_s half-difference
            sig_as = 0.5 * abs(r_hi[key][0] - r_lo[key][0]) if key in r_hi and key in r_lo else 0.0
            # threshold half-difference
            sig_thr = 0.5 * abs(r_thr_up[key][0] - r_thr_dn[key][0]) if key in r_thr_up and key in r_thr_dn else 0.0
            sigma = (sig_as**2 + sig_thr**2) ** 0.5
            out[key] = (rc, sigma)
        return out


