"""Calibration curves.

The idea in plain terms
-----------------------
For each compound we have a handful of calibration runs with known
amounts. In each, the compound's peak is divided by its internal
standard's peak (the "response ratio"), which cancels run-to-run changes
in injection and detector sensitivity. A calibration curve is a line (or a
gentle curve) through response ratio against amount. To quantify a
sample, we read the curve backwards: from its response ratio to an amount.

Why weighting matters. Measurement scatter grows with signal: a peak ten
times bigger scatters roughly ten times more in absolute terms. An
ordinary least-squares fit is dominated by the biggest points, so the low
end of the curve (where most sample values sit) ends up badly placed.
Weighted least squares gives each point a weight of 1, 1/x or 1/x², and
we keep whichever makes the back-calculated amounts closest to the true
ones in relative terms (summed absolute percentage error; Almeida,
Castel-Branco and Falcão, 2002).

Straight or curved? A quadratic always fits at least as well, so we only
accept it when the improvement is bigger than chance would give: Mandel's
test compares the leftover scatter of the two fits with an F-test. The
choice can be overridden.

Through zero? The intercept is tested against zero with a t-test and
reported. The curve keeps its intercept unless you force it through zero.

What to trust. R² is reported but not used for any decision, because it
stays close to 1 even for clearly curved or badly weighted data. The useful
checks are:
- the back-calculated amounts, which should be within ±15% of the truth
  (±20% at the lowest level)
- the residual pattern
- studentised residuals and Cook's distance, which flag single points
  that pull the curve

Flagged points are not removed automatically; excluding one needs a
reason, which is recorded.

Limits. LOD and LOQ are 3.3 and 10 times the scatter at the lowest level,
divided by the slope there. Each quantified amount carries a standard
error from inverse prediction: the scatter of a single new measurement
plus the uncertainty of the curve itself, converted to amount through the
slope.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
from scipy import stats

WEIGHTS = ("1", "1/x", "1/x2")


@dataclass(frozen=True)
class CalibrationSettings:
    weights: str = "auto"          # auto, 1, 1/x, 1/x2
    model: str = "auto"            # auto, linear, quadratic
    force_zero: bool = False
    alpha: float = 0.05            # significance for Mandel's test and the intercept test
    accept_re: float = 15.0        # % back-calculation limit
    accept_re_lowest: float = 20.0
    min_levels: int = 3
    include_zero: bool = False     # use zero-amount points (blanks) in the fit

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "CalibrationSettings":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


def _raw_weight(x, scheme: str, xmin: float) -> np.ndarray:
    xs = np.maximum(np.asarray(x, dtype=float), xmin)
    return {"1": np.ones_like(xs), "1/x": 1.0 / xs, "1/x2": 1.0 / xs ** 2}[scheme]


def _weight(x: np.ndarray, scheme: str) -> np.ndarray:
    """Weights scaled to average 1 (so the residual SD is on the scale of a typical point)."""
    x = np.asarray(x, dtype=float)
    xmin = float(np.min(x[x > 0])) if np.any(x > 0) else 1.0
    w = _raw_weight(x, scheme, xmin)
    return w / w.mean()


def _design(x: np.ndarray, degree: int, zero: bool) -> np.ndarray:
    cols = [] if zero else [np.ones_like(x)]
    cols += [x ** k for k in range(1, degree + 1)]
    return np.column_stack(cols)


@dataclass
class _Fit:
    beta: np.ndarray
    cov: np.ndarray
    s2: float
    dof: int
    resid: np.ndarray
    hat: np.ndarray
    degree: int
    zero: bool


def _wls(x, y, w, degree, zero) -> _Fit:
    X = _design(x, degree, zero)
    sw = np.sqrt(w)
    Xw, yw = X * sw[:, None], y * sw
    XtX_inv = np.linalg.pinv(Xw.T @ Xw)
    beta = XtX_inv @ Xw.T @ yw
    resid = y - X @ beta
    dof = max(len(x) - X.shape[1], 1)
    s2 = float(np.sum(w * resid ** 2) / dof)
    hat = np.einsum("ij,jk,ik->i", Xw, XtX_inv, Xw)
    return _Fit(beta, XtX_inv * s2, s2, dof, resid, hat, degree, zero)


def _coefs(f: _Fit) -> tuple[float, float, float]:
    """(intercept, slope, curvature) whatever the design."""
    b = list(f.beta)
    a0 = 0.0 if f.zero else b.pop(0)
    a1 = b.pop(0) if b else 0.0
    a2 = b.pop(0) if b else 0.0
    return a0, a1, a2


def _invert(a0, a1, a2, y):
    """Amount for response y (the root on the rising branch nearest the data)."""
    y = np.asarray(y, dtype=float)
    if abs(a2) < 1e-15:
        return (y - a0) / a1 if a1 != 0 else np.full_like(y, np.nan)
    disc = a1 ** 2 - 4 * a2 * (a0 - y)
    root = np.sqrt(np.where(disc >= 0, disc, np.nan))
    x1 = (-a1 + root) / (2 * a2)
    x2 = (-a1 - root) / (2 * a2)
    # The branch where the curve rises with x (derivative > 0).
    return np.where(a1 + 2 * a2 * x1 > 0, x1, x2)


@dataclass
class CalibrationFit:
    compound: str
    x: np.ndarray                 # amounts used
    y: np.ndarray                 # responses used
    labels: list[str]             # run names for each point
    excluded: dict[str, str]      # run name -> reason (not used in the fit)
    x_excluded: np.ndarray
    y_excluded: np.ndarray
    weights: str
    model: str                    # linear or quadratic
    zero: bool
    intercept: float
    slope: float
    curvature: float
    cov: np.ndarray
    s: float                      # weighted residual SD (unit weight)
    dof: int
    r2: float
    mandel_p: float
    intercept_p: float
    back_calc: np.ndarray         # amounts recovered from each used point
    rel_error: np.ndarray         # % error of back_calc
    studentised: np.ndarray
    cooks: np.ndarray
    lod: float
    loq: float
    weight_scores: dict[str, float] = field(default_factory=dict)
    settings: CalibrationSettings = field(default_factory=CalibrationSettings)

    # ---------------------------------------------------------------- use
    def response(self, x) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        return self.intercept + self.slope * x + self.curvature * x ** 2

    def amount(self, y) -> np.ndarray:
        return _invert(self.intercept, self.slope, self.curvature, y)

    def _w_at(self, x) -> np.ndarray:
        """Weight a single new point at amount x would carry (same scaling as the fit)."""
        xmin = float(np.min(self.x[self.x > 0])) if np.any(self.x > 0) else 1.0
        mean = _raw_weight(self.x, self.weights, xmin).mean()
        return _raw_weight(x, self.weights, xmin) / mean

    def amount_with_se(self, y):
        """Amount and its standard error by inverse prediction (delta method)."""
        x0 = self.amount(y)
        x0c = np.clip(np.nan_to_num(x0, nan=0.0), 0, None)
        grad = np.column_stack(([] if self.zero else [np.ones_like(x0c)]) +
                               [x0c ** k for k in range(1, (2 if self.model == "quadratic" else 1) + 1)])
        var_fit = np.einsum("ij,jk,ik->i", grad, self.cov, grad)
        var_new = self.s ** 2 / self._w_at(np.where(x0c > 0, x0c, np.nan))
        deriv = self.slope + 2 * self.curvature * x0c
        se = np.sqrt(var_fit + np.nan_to_num(var_new, nan=self.s ** 2)) / np.abs(deriv)
        return x0, se

    # ------------------------------------------------------------ verdicts
    @property
    def n(self) -> int:
        return int(self.x.size)

    def passes(self) -> np.ndarray:
        lim = np.full(self.n, self.settings.accept_re)
        if self.n:
            lim[self.x == self.x.min()] = self.settings.accept_re_lowest
        return np.abs(self.rel_error) <= lim

    @property
    def ok(self) -> bool:
        return self.n >= self.settings.min_levels and bool(np.all(self.passes()))

    def influential(self) -> list[str]:
        """Points worth a second look: |studentised| > 3 or Cook's D > 4/n."""
        flag = (np.abs(self.studentised) > 3) | (self.cooks > 4.0 / max(self.n, 1))
        return [lab for lab, f in zip(self.labels, flag) if f]

    def summary(self) -> str:
        eq = f"y = {self.intercept:.4g} + {self.slope:.4g}·x"
        if self.model == "quadratic":
            eq += f" + {self.curvature:.4g}·x²"
        if self.zero:
            eq = eq.replace(f"{self.intercept:.4g} + ", "")
        w = {"1": "unweighted", "1/x": "weights 1/x", "1/x2": "weights 1/x²"}[self.weights]
        worst = float(np.max(np.abs(self.rel_error))) if self.n else float("nan")
        return (f"{self.model}, {w}; {eq}; R² {self.r2:.4f}; worst back-calculation {worst:.1f}%; "
                f"LOD {self.lod:.3g}, LOQ {self.loq:.3g}")


class CalibrationError(ValueError):
    pass


def fit_calibration(compound: str, x, y, labels=None, excluded: dict | None = None,
                    settings: CalibrationSettings = CalibrationSettings()) -> CalibrationFit:
    """Fit one calibration curve. `excluded` maps point labels to reasons;
    those points are kept for plotting but left out of the fit."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    labels = list(labels) if labels is not None else [str(i) for i in range(x.size)]
    excluded = dict(excluded or {})
    good = np.isfinite(x) & np.isfinite(y)
    if not settings.include_zero:
        good &= x > 0
    use = good & np.array([lab not in excluded for lab in labels], dtype=bool)
    xu, yu = x[use], y[use]
    lu = [lab for lab, u in zip(labels, use) if u]
    if len(np.unique(xu)) < max(settings.min_levels, 2):
        raise CalibrationError(f"{compound}: only {len(np.unique(xu))} usable levels")
    s = settings
    zero = s.force_zero

    def score(scheme, degree):
        f = _wls(xu, yu, _weight(xu, scheme), degree, zero)
        a0, a1, a2 = _coefs(f)
        back = _invert(a0, a1, a2, yu)
        re = 100.0 * (back - xu) / np.where(xu != 0, xu, np.nan)
        return f, back, re, float(np.nansum(np.abs(re[xu > 0])))

    schemes = WEIGHTS if s.weights == "auto" else (s.weights,)
    weight_scores = {sch: score(sch, 1)[3] for sch in schemes}
    scheme = min(weight_scores, key=weight_scores.get)
    w = _weight(xu, scheme)

    lin = _wls(xu, yu, w, 1, zero)
    mandel_p = float("nan")
    quad = None
    if len(np.unique(xu)) >= 4:
        quad = _wls(xu, yu, w, 2, zero)
        ss_lin, ss_quad = lin.s2 * lin.dof, quad.s2 * quad.dof
        if quad.s2 > 0 and quad.dof > 0:
            F = (ss_lin - ss_quad) / quad.s2
            mandel_p = float(stats.f.sf(F, 1, quad.dof))
    if s.model == "quadratic" and quad is not None:
        f = quad
    elif s.model == "auto" and quad is not None and mandel_p < s.alpha:
        f = quad
    else:
        f = lin
    a0, a1, a2 = _coefs(f)
    back = _invert(a0, a1, a2, yu)
    re = 100.0 * (back - xu) / np.where(xu != 0, xu, np.nan)

    if zero:
        intercept_p = float("nan")
    else:
        se0 = np.sqrt(f.cov[0, 0])
        intercept_p = float(2 * stats.t.sf(abs(a0 / se0), f.dof)) if se0 > 0 else float("nan")

    ybar = np.average(yu, weights=w)
    ss_tot = float(np.sum(w * (yu - ybar) ** 2))
    r2 = 1.0 - f.s2 * f.dof / ss_tot if ss_tot > 0 else float("nan")

    h = np.clip(f.hat, 0, 0.9999)
    p = len(f.beta)
    std_res = np.sqrt(w) * f.resid / np.sqrt(f.s2 * (1 - h)) if f.s2 > 0 else np.zeros_like(h)
    cooks = std_res ** 2 * h / (p * (1 - h))

    # Scatter at the lowest non-zero level, in response units.
    pos = xu[xu > 0]
    x_low = pos.min() if pos.size else 1.0
    w_low = float(_raw_weight(x_low, scheme, x_low) / _raw_weight(xu, scheme, x_low).mean())
    s_low = float(np.sqrt(f.s2 / w_low))
    slope_low = a1 + 2 * a2 * x_low
    lod = float(3.3 * s_low / slope_low) if slope_low > 0 else float("nan")
    loq = float(10.0 * s_low / slope_low) if slope_low > 0 else float("nan")

    ex_mask = good & ~use
    return CalibrationFit(
        compound=compound, x=xu, y=yu, labels=lu, excluded=excluded,
        x_excluded=x[ex_mask], y_excluded=y[ex_mask], weights=scheme,
        model="quadratic" if f.degree == 2 else "linear", zero=zero,
        intercept=a0, slope=a1, curvature=a2, cov=f.cov, s=float(np.sqrt(f.s2)), dof=f.dof,
        r2=float(r2), mandel_p=mandel_p, intercept_p=intercept_p, back_calc=back, rel_error=re,
        studentised=std_res, cooks=cooks, lod=lod, loq=loq, weight_scores=weight_scores,
        settings=settings)

