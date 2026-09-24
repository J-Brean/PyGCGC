"""Compound-class assignment for detected blobs.

The idea in plain terms
------------------------
Blob detection tells you where the peaks are. Nothing yet says whether a
peak is an alkane, an alkylbenzene, a PAH, a nitrogen heterocycle, a
siloxane... `core/identify.py` answers "which compound?" from a library; this
answers "which kind of compound?" from the spectrum alone, so it works on
compounds that are in no library.

Rules live in a TOML file (`presets/chemical_classes.toml`, read that file
first: it is the documentation) that says, for each class, what evidence to
look for:

    [[class]]
    name   = "alkylbenzene"
    need   = ["tropylium: 0.15 -> 0.30", "aryl: 0.02 -> 0.05"]
    prefer = ["dt2: 0.3 -> 1.0"]

`tropylium` is a named set of ions (91, 105, 119, 133); its "share" is the
fraction of the blob's ion current on them. `0.15 -> 0.30` reads "not
satisfied at 0.15, fully satisfied at 0.30, linear in between".

Scoring
-------
Every criterion becomes a number from 0 to 1. A class scores its weakest
`need` (a fuzzy AND, so one missing piece of key evidence sinks the class),
scaled from 0.5 to 1.0 by how well the `prefer` list is met (supporting
evidence: it can lift a class, never block it). Each blob gets a score for
every class; the best one wins if it clears a floor, otherwise the blob is
left "unclassified" rather than forced into a class.

Confidence tiers: confident (score >= 0.75), probable (>= 0.5), possible
(>= 0.3). A tier drops one notch when the runner-up class in a *different*
group is close (a call between alkene and cycloalkane stays confident as a
hydrocarbon; one between ketone and amine does not), and when the blob was
flagged as co-eluting, because its spectrum is then a mixture.

A criterion whose feature is unknown for a blob (no alkane ladder yet, no
molecular-ion candidate, an ion outside the run's m/z range) is skipped and
shown as n/a in `ClassResult.explain`. It is never counted as a failure.

Limits, stated once
--------------------
Unit-mass EI cannot tell every class apart. Alkenes and cycloalkanes share
their ion series. Aliphatic amines and ketones share 58, 72 and 86. m/z 29 to
31, the diagnostic ions of alcohols, ethers and primary amines, are below this
instrument's range. Where classes overlap the result says so (a small margin,
a lowered tier) instead of picking one quietly.

Tandem acquisitions (`soft_spectra_for_blobs`)
-----------------------------------------------
A 70 eV / 14 eV tandem file holds the same injection twice. At 14 eV the
molecular ion survives, so the same blobs can be given a trustworthy M, its
parity (the nitrogen rule) and its series (M mod 14). That settles amine
against ketone when M is found and corroborates hydrocarbon classes; it cannot
separate an alkene from a cycloalkane, which have the same formula. On the
supplied run 86% of confident alkylbenzene calls have M = 8 (mod 14), the
CnH2n-6 series, from ions that share nothing with the 70 eV call.
"""

from __future__ import annotations

import difflib
import re
import tomllib
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np

from .fold import Fold
from .hydrocarbons import AlkaneSeries
from .peaks import BlobSet, _spectra_by_label, make_frame
from .run import Run
from .siloxane import CYCLIC_REFERENCE, SiloxaneSettings, siloxane_score

DEFAULT_RULES = Path(__file__).resolve().parents[1] / "presets" / "chemical_classes.toml"
MOLECULAR = ("M_mass", "M_rel", "M_share", "M_odd", "M_plus2")
BUILTIN_FEATURES = MOLECULAR + ("siloxane", "dt2", "soft_M_mod14") + tuple("soft_" + m for m in MOLECULAR)
TIERS = ("possible", "probable", "confident")

_RAMP = re.compile(r"^\s*(\w+)\s*:\s*(-?[\d.]+)\s*->\s*(-?[\d.]+)\s*$")
_HARD = re.compile(r"^\s*(\w+)\s*(>=|<=|==|!=|>|<)\s*(-?[\d.]+)\s*$")


@dataclass(frozen=True)
class ClassSettings:
    m_min_rel: float = 0.08       # M candidate: highest ion at least this fraction of the base peak...
    soft_m_min_rel: float = 0.03  # the same, for the low-energy spectrum: an n-alkane M+ is only 4 to 7% of its base peak there
    m_min_mass: int = 60          # ...and above this m/z (keeps CO2 at 44 and other background from being "M")
    isotope_ratio: float = 3.0    # step M down a mass if the ion below is this many times taller (13C isotope)
    min_ions: int = 5             # a spectrum with fewer non-zero ions is left unclassified
    prefer_floor: float = 0.5     # score of a class whose needs are all met and whose preferences all fail
    confident: float = 0.75
    probable: float = 0.5
    possible: float = 0.3
    min_margin: float = 0.15      # best minus best-of-another-group below this drops the tier a notch
    ladder_slack_carbons: float = 1.0   # dt2 is unknown further than this outside the alkane ladder

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "ClassSettings":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


# ------------------------------------------------------------------ rules
@dataclass(frozen=True)
class Criterion:
    text: str
    feature: str
    a: float
    b: float
    op: str = ""          # "" for a ramp; ">=", "<=", "==", "!=", ">" or "<" for a hard cut

    def membership(self, x: np.ndarray) -> np.ndarray:
        """0 to 1 for each value of the feature; NaN stays NaN (unknown)."""
        x = np.asarray(x, dtype=float)
        if self.op:
            ok = {">=": x >= self.a, "<=": x <= self.a, "==": x == self.a, "!=": x != self.a, ">": x > self.a, "<": x < self.a}[self.op]
            return np.where(np.isnan(x), np.nan, ok.astype(float))
        return np.clip((x - self.a) / (self.b - self.a), 0.0, 1.0)


@dataclass
class ChemClass:
    name: str
    group: str
    note: str
    need: list
    prefer: list
    enabled: bool = True     # a disabled class is read and checked but never scores


@dataclass
class RuleSet:
    ions: dict            # name -> tuple of m/z
    ratios: dict          # name -> (numerator set, denominator set)
    losses: dict          # name -> mass lost from M
    classes: list

    def feature_names(self) -> list:
        return (list(BUILTIN_FEATURES) + list(self.ions) + list(self.ratios)
                + [f"loss_{k}" for k in self.losses] + [f"soft_loss_{k}" for k in self.losses])


def _parse_criterion(text: str, where: str) -> Criterion:
    m = _RAMP.match(text)
    if m:
        a, b = float(m.group(2)), float(m.group(3))
        if a == b:
            raise ValueError(f"{where}: '{text}' has the same start and end; use '>= {a:g}' for a hard cut")
        return Criterion(text=text.strip(), feature=m.group(1), a=a, b=b)
    m = _HARD.match(text)
    if m:
        return Criterion(text=text.strip(), feature=m.group(1), a=float(m.group(3)), b=float(m.group(3)), op=m.group(2))
    raise ValueError(f"{where}: can't read '{text}'. Write 'name: a -> b' or 'name >= value'.")


def load_rules(path=None) -> RuleSet:
    """Read a rule file (default: the shipped `presets/chemical_classes.toml`).
    Every mistake is reported with the class and the criterion it is in."""
    path = Path(path) if path else DEFAULT_RULES
    raw = tomllib.loads(path.read_text(encoding="utf-8-sig"))     # -sig: Windows editors add a BOM that tomllib rejects
    ions = {k: tuple(int(m) for m in v) for k, v in raw.get("ions", {}).items()}
    ratios = {}
    for k, v in raw.get("ratios", {}).items():
        if len(v) != 2 or any(s not in ions for s in v):
            raise ValueError(f"[ratios] {k}: needs two names from [ions], got {v}")
        ratios[k] = (v[0], v[1])
    losses = {k: int(v) for k, v in raw.get("losses", {}).items()}
    rules = RuleSet(ions=ions, ratios=ratios, losses=losses, classes=[])
    valid = rules.feature_names()
    for entry in raw.get("class", []):
        name = entry.get("name", "").strip()
        if not name:
            raise ValueError(f"{path.name}: a [[class]] has no name")
        if any(c.name == name for c in rules.classes):
            raise ValueError(f"{path.name}: class '{name}' appears twice")
        crit = {}
        for key in ("need", "prefer"):
            crit[key] = [_parse_criterion(t, f"class '{name}' {key}") for t in entry.get(key, [])]
            for c in crit[key]:
                if c.feature not in valid:
                    close = difflib.get_close_matches(c.feature, valid, n=3)
                    hint = f" Did you mean {', '.join(close)}?" if close else ""
                    raise ValueError(f"class '{name}' {key}: unknown feature '{c.feature}' in '{c.text}'.{hint}")
        if not crit["need"]:
            raise ValueError(f"class '{name}' has no 'need' criteria, so it would match everything")
        rules.classes.append(ChemClass(name=name, group=entry.get("group", "").strip() or name,
                                       note=entry.get("note", ""), need=crit["need"], prefer=crit["prefer"],
                                       enabled=bool(entry.get("enabled", True))))
    return rules


# ------------------------------------------------------------------ features
def _molecular_ion(spec: np.ndarray, mz: np.ndarray, s: ClassSettings) -> int:
    """Index of the molecular-ion candidate, or -1."""
    base = spec.max()
    if base <= 0:
        return -1
    idx = np.flatnonzero((spec >= s.m_min_rel * base) & (mz >= s.m_min_mass))
    if idx.size == 0:
        return -1
    j = int(idx[-1])
    for _ in range(2):      # an ion much smaller than the one below it is probably its 13C isotope
        if j > 0 and spec[j - 1] >= s.isotope_ratio * spec[j]:
            j -= 1
    return j


def _molecular_features(spectra: np.ndarray, mz: np.ndarray, losses: dict, settings: ClassSettings) -> dict:
    """M_mass, M_rel, M_share, M_odd, M_plus2 and loss_<name> for every spectrum (NaN where unknown)."""
    n = spectra.shape[0]
    total = spectra.sum(axis=1)
    live = total > 0
    m_idx = np.array([_molecular_ion(spectra[i], mz, settings) for i in range(n)])
    have_m = m_idx >= 0
    safe = np.where(have_m, m_idx, 0)
    m_int = np.where(have_m, spectra[np.arange(n), safe], np.nan)
    out = {"M_mass": np.where(have_m, mz[safe], np.nan).astype(float),
           "M_rel": np.where(have_m, m_int / np.where(live, spectra.max(axis=1), np.nan), np.nan),
           "M_share": np.where(have_m & live, m_int / np.where(live, total, 1.0), np.nan),
           "M_odd": np.where(have_m, (mz[safe] % 2).astype(float), np.nan)}

    def at_offset(offset: int):
        """I(M + offset) / I(M), NaN where M is unknown or the ion is outside the range."""
        pos = safe + offset
        ok = have_m & (pos >= 0) & (pos < mz.size) & (mz[np.clip(pos, 0, mz.size - 1)] == mz[safe] + offset)
        return np.where(ok, spectra[np.arange(n), np.clip(pos, 0, mz.size - 1)] / np.where(m_int > 0, m_int, np.nan), np.nan)

    out["M_plus2"] = at_offset(2)
    for name, lost in losses.items():
        out[f"loss_{name}"] = at_offset(-lost)
    return out


def compute_features(spectra: np.ndarray, mz: np.ndarray, rules: RuleSet, t1=None, t2=None,
                     series: AlkaneSeries | None = None, settings: ClassSettings = ClassSettings(),
                     soft_spectra: np.ndarray | None = None, soft_mz: np.ndarray | None = None) -> dict:
    """Every named feature, one array per feature with one value per spectrum
    (NaN where it can't be known). `soft_spectra` (same rows, on `soft_mz`) is
    the same blobs at a lower ionisation energy, where the molecular ion
    survives; it feeds the `soft_*` features and nothing else."""
    spectra = np.clip(np.asarray(spectra, dtype=float), 0, None)
    mz = np.asarray(mz, dtype=np.int64)
    n = spectra.shape[0]
    where = {int(m): j for j, m in enumerate(mz)}
    total = spectra.sum(axis=1)
    live = total > 0
    feats: dict = {}

    def share(ions):
        cols = [where[m] for m in ions if m in where]
        if not cols:
            return np.full(n, np.nan)
        return np.where(live, spectra[:, cols].sum(axis=1) / np.where(live, total, 1.0), np.nan)

    for name, ions in rules.ions.items():
        feats[name] = share(ions)
    for name, (num, den) in rules.ratios.items():
        # Where both sides are negligible the ratio is noise, not a measurement.
        den_share = np.maximum(feats[den], 0.005)
        feats[name] = np.where((feats[num] < 0.005) & (feats[den] < 0.005), np.nan, feats[num] / den_share)

    feats.update(_molecular_features(spectra, mz, rules.losses, settings))

    soft_names = ["soft_" + m for m in MOLECULAR] + [f"soft_loss_{k}" for k in rules.losses] + ["soft_M_mod14"]
    if soft_spectra is not None:
        soft = np.clip(np.asarray(soft_spectra, dtype=float), 0, None)
        soft_axis = np.asarray(soft_mz, dtype=np.int64)
        soft_settings = replace(settings, m_min_rel=settings.soft_m_min_rel)
        for name, v in _molecular_features(soft, soft_axis, rules.losses, soft_settings).items():
            feats["soft_" + name] = v
        feats["soft_M_mod14"] = feats["soft_M_mass"] % 14
        soft_sparse = (soft > 0).sum(axis=1) < settings.min_ions
        for name in soft_names:
            feats[name] = np.where(soft_sparse, np.nan, feats[name])
    else:
        for name in soft_names:
            feats[name] = np.full(n, np.nan)
    bleed, cyclic = SiloxaneSettings(), SiloxaneSettings.from_reference(CYCLIC_REFERENCE)
    mz0 = int(mz[0])
    feats["siloxane"] = np.array([max(siloxane_score(spectra[i], mz0, bleed), siloxane_score(spectra[i], mz0, cyclic))
                                  if live[i] else np.nan for i in range(n)])

    feats["dt2"] = np.full(n, np.nan)
    if series is not None and len(series) >= 2 and t1 is not None and t2 is not None:
        t1, t2 = np.asarray(t1, dtype=float), np.asarray(t2, dtype=float)
        carbons = series.carbons()
        for i in range(n):
            c = series.carbon_of(float(t1[i]))
            if carbons.min() - settings.ladder_slack_carbons <= c <= carbons.max() + settings.ladder_slack_carbons:
                feats["dt2"][i] = t2[i] - series.t2_of(c)

    # Too little signal to say anything about: everything unknown.
    sparse = (spectra > 0).sum(axis=1) < settings.min_ions
    for name in feats:
        if not name.startswith("soft_"):
            feats[name] = np.where(sparse, np.nan, feats[name])
    return feats


# ------------------------------------------------------------------ scoring
def _class_score(cls: ChemClass, feats: dict, s: ClassSettings):
    need = np.array([c.membership(feats[c.feature]) for c in cls.need])          # (k, n)
    known = ~np.isnan(need)
    weakest = np.where(known, need, np.inf).min(axis=0)
    weakest = np.where(known.any(axis=0), weakest, 0.0)                          # no needs known: no claim
    if cls.prefer:
        prefer = np.array([c.membership(feats[c.feature]) for c in cls.prefer])
        pk = ~np.isnan(prefer)
        mean = np.where(pk.any(axis=0), np.nansum(np.where(pk, prefer, 0.0), axis=0) / np.maximum(pk.sum(axis=0), 1), 1.0)
        factor = s.prefer_floor + (1.0 - s.prefer_floor) * mean
    else:
        prefer, factor = np.empty((0, need.shape[1])), np.ones(need.shape[1])
    return weakest * factor, need, prefer


@dataclass
class ClassCall:
    """One blob's result."""
    index: int
    cls: str | None            # None: unclassified
    group: str | None
    score: float
    tier: str                  # confident, probable, possible or unclassified
    runner_up: str | None
    runner_up_score: float
    margin: float              # best minus the best class in another group
    limiting: str | None       # the need that held the score down, if it is below 1
    co_elution: bool = False


@dataclass
class ClassResult:
    rules: RuleSet
    classes: list                           # the enabled classes, in score-column order
    settings: ClassSettings
    features: dict
    scores: np.ndarray                      # (n_blobs, len(classes))
    calls: list
    _detail: list = field(default_factory=list, repr=False)   # per class: (need memberships, prefer memberships)

    def __len__(self) -> int:
        return len(self.calls)

    def explain(self, i: int) -> str:
        """Why blob i got the class it did: every criterion, its value and its score."""
        call = self.calls[i]
        head = (f"blob {i}: {call.cls} ({call.score:.2f}, {call.tier})" if call.cls
                else f"blob {i}: unclassified (best score {call.score:.2f})")
        lines = [head]
        if call.runner_up:
            lines.append(f"  runner-up {call.runner_up} ({call.runner_up_score:.2f}), margin to another group {call.margin:.2f}")
        if call.co_elution:
            lines.append("  flagged as co-eluting: the spectrum is a mixture, so the tier is lowered")
        k = self.scores[i].argmax()
        cls = self.classes[k]
        need_m, prefer_m = self._detail[k]
        for label, crits, memb in (("need", cls.need, need_m), ("prefer", cls.prefer, prefer_m)):
            for c, m in zip(crits, memb):
                v = self.features[c.feature][i]
                lines.append(f"  {label:6s} {c.text:32s} value {'n/a' if np.isnan(v) else f'{v:.3g}':>7s}  -> "
                             f"{'n/a' if np.isnan(m[i]) else f'{m[i]:.2f}'}")
        return "\n".join(lines)

    def summary(self, volume: np.ndarray, by: str = "class") -> list:
        """(name, n blobs, volume, percent of all volume), largest first.
        `by` is "class" or "group"; unclassified blobs are listed too."""
        volume = np.asarray(volume, dtype=float)
        key = [(getattr(c, "cls" if by == "class" else "group") or "unclassified") for c in self.calls]
        total = float(volume.sum()) or 1.0
        rows = []
        for name in sorted(set(key)):
            sel = np.array([k == name for k in key])
            rows.append((name, int(sel.sum()), float(volume[sel].sum()), 100.0 * float(volume[sel].sum()) / total))
        return sorted(rows, key=lambda r: -r[2])


def classify_spectra(spectra: np.ndarray, mz: np.ndarray, rules: RuleSet | None = None, t1=None, t2=None,
                     series: AlkaneSeries | None = None, co_elution=None,
                     settings: ClassSettings = ClassSettings(), soft_spectra=None, soft_mz=None) -> ClassResult:
    """Class of every row of `spectra` (n, n_mz), on the integer m/z axis `mz`.
    `t1`, `t2` (s) and an alkane `series` are optional; without them `dt2` is
    unknown and only the spectrum decides. `soft_spectra` and `soft_mz` are the
    same rows measured at a low ionisation energy (see `soft_spectra_for_blobs`)."""
    rules = rules or load_rules()
    feats = compute_features(spectra, mz, rules, t1, t2, series, settings, soft_spectra, soft_mz)
    n = np.asarray(spectra).shape[0]
    co = np.zeros(n, dtype=bool) if co_elution is None or np.size(co_elution) != n else np.asarray(co_elution, dtype=bool)

    classes = [c for c in rules.classes if c.enabled]
    scores, detail = np.zeros((n, len(classes))), []
    for k, cls in enumerate(classes):
        scores[:, k], need, prefer = _class_score(cls, feats, settings)
        detail.append((need, prefer))
    groups = np.array([c.group for c in classes])

    calls = []
    for i in range(n):
        order = np.argsort(scores[i])[::-1]
        b = int(order[0])
        best = float(scores[i, b])
        other = scores[i][groups != groups[b]]
        rival = float(other.max()) if other.size else 0.0
        second = next((int(j) for j in order[1:] if scores[i, j] > 0), None)
        if best < settings.possible:
            calls.append(ClassCall(i, None, None, best, "unclassified", classes[b].name if best > 0 else None,
                                   best, best - rival, None, bool(co[i])))
            continue
        tier = "confident" if best >= settings.confident else "probable" if best >= settings.probable else "possible"
        drops = int(best - rival < settings.min_margin) + int(co[i])
        tier = TIERS[max(TIERS.index(tier) - drops, 0)]
        need_m = detail[b][0][:, i]
        weak = int(np.nanargmin(np.where(np.isnan(need_m), np.inf, need_m))) if (~np.isnan(need_m)).any() else None
        limiting = classes[b].need[weak].text if weak is not None and need_m[weak] < 0.999 else None
        calls.append(ClassCall(i, classes[b].name, classes[b].group, best, tier,
                               classes[second].name if second is not None else None,
                               float(scores[i, second]) if second is not None else 0.0,
                               best - rival, limiting, bool(co[i])))
    return ClassResult(rules=rules, classes=classes, settings=settings, features=feats, scores=scores, calls=calls,
                       _detail=detail)


def classify_blobs(blobs: BlobSet, rules: RuleSet | None = None, series: AlkaneSeries | None = None,
                   settings: ClassSettings = ClassSettings(), soft=None) -> ClassResult:
    """Class of every blob in a `BlobSet`, from its apex spectrum (the cleanest
    one available; a footprint sum picks up a co-eluting neighbour). A blob
    whose apex spectrum is empty falls back to its footprint spectrum: on the
    supplied run 319 of 4068 blobs (11.7% of the volume, all high S/N) have no
    apex spectrum at all. Pass the run's alkane `series` to switch on the
    t2-above-the-alkane-line evidence. `soft` is `(spectra, mz)` from
    `soft_spectra_for_blobs` for the same run at a lower ionisation energy."""
    use_apex = blobs.apex_spectra.shape[0] == blobs.t1.size and blobs.apex_spectra.size
    spectra = blobs.spectra
    if use_apex:
        apex_empty = (blobs.apex_spectra > 0).sum(axis=1) < settings.min_ions
        spectra = np.where(apex_empty[:, None], blobs.spectra, blobs.apex_spectra)
    soft_spectra, soft_mz = soft if soft is not None else (None, None)
    return classify_spectra(spectra, blobs.mz, rules, blobs.t1, blobs.t2, series, blobs.co_elution, settings,
                            soft_spectra, soft_mz)


def soft_spectra_for_blobs(soft_run: Run, fold: Fold, blobs: BlobSet, settings: ClassSettings = ClassSettings()):
    """Spectra of `blobs` (found on a 70 eV run) measured at the lower energy
    of a tandem acquisition: `(spectra, mz)`, one row per blob.

    The two channels are interleaved scans of one injection, about 10 ms
    apart, so a blob's pixels on the 70 eV grid are the pixels of the same
    material on the 14 eV grid. `fold` is the 70 eV run's; the 14 eV run is
    folded with the same resolved spec (same absolute t0, period and phase).
    Each row is the nearest-scan apex spectrum, or the footprint sum where
    that has fewer than `settings.min_ions` ions (the same fallback as
    `classify_blobs`). `soft_run` should have had its background removed.
    """
    soft_fold = Fold(soft_run.times, fold.spec)
    frame = make_frame(soft_fold, blobs.settings)
    if frame.rows != blobs.labels.shape[0]:
        raise ValueError(f"the {soft_run.name} run folds to {frame.rows} rows, the blobs to {blobs.labels.shape[0]}: "
                         "they are not the same acquisition")
    n = int(blobs.row.size)
    footprint, _ = _spectra_by_label(soft_run, blobs.labels, frame, n)
    apex_labels = np.zeros_like(blobs.labels)
    apex_labels[blobs.row, blobs.col] = np.arange(1, n + 1)
    apex, _ = _spectra_by_label(soft_run, apex_labels, frame, n)
    thin = (apex > 0).sum(axis=1) < settings.min_ions
    return np.where(thin[:, None], footprint, apex), soft_run.mz

def write_classes_csv(result: ClassResult, blobs, path) -> Path:
    """One row per blob: where it is, how big, what class, how sure and why not more sure."""
    import csv
    path = Path(path)
    f = result.features
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["blob", "t1_min", "t2_s", "volume", "class", "group", "tier", "score", "runner_up", "runner_up_score",
                    "margin", "held_back_by", "co_eluting", "dt2_s", "M_70eV", "M_14eV"])
        for i, c in enumerate(result.calls):
            w.writerow([i, f"{blobs.t1[i] / 60.0:.4f}", f"{blobs.t2[i]:.3f}", f"{blobs.volume[i]:.6g}", c.cls or "", c.group or "",
                        c.tier, f"{c.score:.3f}", c.runner_up or "", f"{c.runner_up_score:.3f}", f"{c.margin:.3f}",
                        c.limiting or "", int(c.co_elution), *(("" if np.isnan(v) else f"{v:.6g}")
                                                                for v in (f["dt2"][i], f["M_mass"][i], f["soft_M_mass"][i]))])
    return path
