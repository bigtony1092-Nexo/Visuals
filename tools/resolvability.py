#!/usr/bin/env python3
"""Who can be told apart from whom, in units of measurement sigma.

In all-ion fragmentation nothing is isolated: every precursor fragments at once.
A fragment carries its own m/z plus its PARENT's retention time and drift — it
does not carry the parent's m/z. So assigning a fragment to a parent needs the
parents resolved in (rt, drift) alone. Precursor m/z helps only indirectly, by
inducing a drift difference: at fixed charge drift tracks m/z closely, so two
precursors apart in m/z are apart in drift, and that difference is inherited by
their fragments. The correlation between m/z and drift is the mechanism here,
not a redundancy.

Everything below is a RELATIVE distance. An absolute tolerance says nothing
without the measurement scale beside it: 7 drift bins is eight sigma on one
instrument and noise on another. Separations are reported in sigma, and sigma is
measured from the data where the data can supply it.

    python3 tools/resolvability.py ions.json --sigma-drift 0.9 --sigma-rt 0.15
    python3 tools/resolvability.py ions.json --estimate-sigma
    python3 tools/resolvability.py ions.json --estimate-sigma --k 2.0 --worst 10

Input is JSON lines or a JSON list, one object per precursor observation:

    {"run": "...", "uid": "SEQ|2", "mz": 614.34, "z": 2,
     "drift": 69.4, "rt": 27.8}

`uid` identifies the same precursor across runs, which is what --estimate-sigma
uses: the scatter of one precursor's own repeated measurements is the resolution
scale. That estimate includes run-to-run calibration drift, so it OVERSTATES the
within-run peak width and the occupancy it reports is pessimistic.

Exit status is 1 if any precursor is unresolved from another, because those
fragments cannot be assigned by coordinates alone and need a different handle.
"""

import argparse
import json
import math
import statistics as st
import sys
from collections import defaultdict

K_DEFAULT = 2.0      # separations below this many sigma are called unresolved
# A sigma that is a large share of the axis's whole range is not a resolution,
# it is a confound. Retention time across runs of different gradient length is
# the standard case: the same peptide elutes at minute 28 in a 60-minute run and
# minute 56 in a 120-minute one, and the "scatter" is the gradient, not the peak.
SIGMA_RANGE_WARN = 0.10


def load(path):
    text = open(path).read().strip()
    if text.startswith("["):
        rows = json.loads(text)
    else:
        rows = [json.loads(l) for l in text.splitlines() if l.strip()]
    out = []
    for r in rows:
        try:
            out.append({"run": str(r.get("run", "")), "uid": str(r["uid"]),
                        "mz": float(r["mz"]), "z": int(r.get("z", 0)),
                        "drift": float(r["drift"]), "rt": float(r["rt"])})
        except (KeyError, TypeError, ValueError):
            continue
    return out


def estimate_sigma(rows, axis):
    """Resolution scale from one precursor's own repeated measurements.

    Returns (sigma, n_uids). Pooled within-uid SD, which is the right quantity:
    the spread of DIFFERENT precursors on an axis says how predictable that axis
    is, not how finely two of them can be told apart.
    """
    by_uid = defaultdict(list)
    for r in rows:
        by_uid[r["uid"]].append(r[axis])
    num = den = 0.0
    used = 0
    for v in by_uid.values():
        if len(v) < 3:
            continue
        s = st.stdev(v)
        num += (len(v) - 1) * s * s
        den += len(v) - 1
        used += 1
    if den <= 0:
        return None, 0
    sigma = math.sqrt(num / den)
    # A pooled SD of zero is not a fine instrument, it is an axis that did not
    # vary across the repeats — a run-level constant, or a value rounded flat.
    # Dividing by it would make every pair infinitely separated.
    if sigma <= 0:
        return 0.0, used
    return sigma, used


def rank_rt(rows):
    """Replace rt with its rank within each run.

    MARCH's own rule: minutes change between runs, the ORDER does not. Across
    runs of different gradient length the raw minute is incomparable, so any
    cross-run scale built on it measures the gradient rather than the peak.
    """
    by_run = defaultdict(list)
    for r in rows:
        by_run[r["run"]].append(r)
    for group in by_run.values():
        group.sort(key=lambda r: r["rt"])
        n = len(group)
        for i, r in enumerate(group):
            # Scaled to 0-100 so the number is comparable between runs of
            # different length, and so --sigma-rt reads as "percent of elution".
            r["rt"] = 100.0 * i / (n - 1) if n > 1 else 50.0
    return rows


def occupancy(rows, sd, sr, k):
    """For each precursor, which others sit within k sigma on BOTH axes.

    Both, not either: a pair separated on one axis is resolved, so the collision
    set is the intersection. Charge is not used — fragments do not carry it any
    more than they carry the parent's m/z.
    """
    # One coordinate per precursor, as the MEDIAN of its observations. Keeping
    # the first row seen compared coordinates from different runs and made the
    # verdict depend on input line order.
    grouped = defaultdict(list)
    for r in rows:
        grouped[r["uid"]].append(r)
    ions = []
    for uid, group in grouped.items():
        ions.append({"uid": uid, "z": group[0]["z"],
                     "mz": st.median([g["mz"] for g in group]),
                     "drift": st.median([g["drift"] for g in group]),
                     "rt": st.median([g["rt"] for g in group]),
                     "n_obs": len(group)})
    ions.sort(key=lambda r: r["drift"])
    out = {}
    for i, a in enumerate(ions):
        hits = []
        for b in ions:
            if b is a:
                continue
            dd = abs(a["drift"] - b["drift"]) / sd
            dr = abs(a["rt"] - b["rt"]) / sr
            if dd < k and dr < k:
                hits.append((math.hypot(dd, dr), b["uid"], dd, dr))
        hits.sort()
        out[a["uid"]] = hits
    return ions, out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("ions", help="JSON or JSON-lines precursor observations")
    p.add_argument("--sigma-drift", type=float, help="drift measurement SD, same units as `drift`")
    p.add_argument("--sigma-rt", type=float, help="retention measurement SD, same units as `rt`")
    p.add_argument("--estimate-sigma", action="store_true",
                   help="measure both from repeated observations of the same uid")
    p.add_argument("--k", type=float, default=K_DEFAULT,
                   help="separations below this many sigma count as unresolved (default %.1f)" % K_DEFAULT)
    p.add_argument("--worst", type=int, default=8, help="list this many most-crowded precursors")
    p.add_argument("--rank-rt", action="store_true",
                   help="replace rt with within-run elution rank, scaled 0-100, before"
                        " anything else; use whenever runs differ in gradient length")
    args = p.parse_args(argv)

    rows = load(args.ions)
    if not rows:
        print("No usable precursor rows.", file=sys.stderr)
        return 1

    if args.rank_rt:
        rows = rank_rt(rows)
    sd, sr = args.sigma_drift, args.sigma_rt
    notes = []
    spans = {axis: (max(r[axis] for r in rows) - min(r[axis] for r in rows))
             for axis in ("drift", "rt")}
    if args.estimate_sigma:
        for axis, name in (("drift", "sigma_drift"), ("rt", "sigma_rt")):
            est, used = estimate_sigma(rows, axis)
            if est is None:
                notes.append("could not estimate %s: no uid has 3+ observations" % name)
            elif est == 0.0:
                notes.append("%s estimated at ZERO from %d precursors: that axis did not"
                             % (name, used))
                notes.append("  vary across the repeats at all, so it carries no measurement")
                notes.append("  scale that can be measured.")
                if axis == "rt" and args.rank_rt:
                    # Not a defect: --rank-rt reduces rt to elution order, and an
                    # order that never changes between runs is the best case. There
                    # is simply nothing to estimate from, so the tolerance has to be
                    # chosen — how much rank slip counts as the same position —
                    # rather than measured.
                    notes.append("  With --rank-rt that means the elution ORDER was identical in")
                    notes.append("  every run, which is the best case and not a problem. Rank has")
                    notes.append("  no scatter to measure, so choose the tolerance instead: pass")
                    notes.append("  --sigma-rt as the rank slip you count as the same position")
                    notes.append("  (1.0 is one position out of the 0-100 scale).")
                else:
                    notes.append("  A run-level constant or a coarsely rounded value does this.")
                    notes.append("  Supply --%s explicitly or fix the source."
                                 % name.replace("_", "-"))
            else:
                notes.append("%s estimated at %.4f from %d precursors seen 3+ times"
                             % (name, est, used))
                span = spans[axis]
                if span > 0 and est / span > SIGMA_RANGE_WARN:
                    notes.append("  REFUSED: that is %.0f%% of the whole %s range (%.3f)."
                                 % (100 * est / span, axis, span))
                    notes.append("  A scale that large is not a peak width, it is a")
                    notes.append("  confound inside the repeats. For rt across runs of")
                    notes.append("  different gradient length this is the gradient itself:")
                    notes.append("  pass --rank-rt to use elution order instead, as MARCH")
                    notes.append("  does, or estimate within a single run.")
                    continue
                if axis == "drift" and sd is None:
                    sd = est
                if axis == "rt" and sr is None:
                    sr = est
    if not sd or not sr:
        print("Need a measurement scale on both axes: pass --sigma-drift and "
              "--sigma-rt, or --estimate-sigma over data with repeated uids.",
              file=sys.stderr)
        for n in notes:
            print("  " + n, file=sys.stderr)
        return 1

    ions, hits = occupancy(rows, sd, sr, args.k)
    counts = sorted(len(v) for v in hits.values())
    n = len(ions)
    clean = sum(1 for c in counts if c == 0)

    out = []
    out.append("RESOLVABILITY IN (rt, drift), THE SPACE A FRAGMENT CAN BE ASSIGNED IN")
    out.append("  %d distinct precursors from %d observations" % (n, len(rows)))
    out.append("  sigma_drift %.4f   sigma_rt %.4f   threshold %.1f sigma on both axes"
               % (sd, sr, args.k))
    for note in notes:
        out.append("  " + note)
    if args.estimate_sigma:
        out.append("  Estimated sigma includes run-to-run calibration drift, so it is an")
        out.append("  upper bound on the within-run peak width and this occupancy is")
        out.append("  pessimistic: real resolution is finer and collisions fewer.")
    out.append("")
    out.append("  uniquely resolved        %5d  (%.1f%%)" % (clean, 100.0 * clean / n))
    out.append("  sharing with 1+ others   %5d  (%.1f%%)" % (n - clean, 100.0 * (n - clean) / n))
    if counts:
        out.append("  co-residents per precursor: median %d, p90 %d, max %d"
                   % (counts[len(counts) // 2], counts[min(int(0.9 * len(counts)), len(counts) - 1)],
                      counts[-1]))
    out.append("")
    out.append("  Those %.1f%% are the fragments that coordinates alone cannot assign."
               % (100.0 * (n - clean) / n))
    out.append("  They are what a fragment-physics prior has to resolve, and the count")
    out.append("  is the honest measure of how hard that job is in this data.")

    crowded = sorted(hits.items(), key=lambda kv: -len(kv[1]))[:args.worst]
    if crowded and crowded[0][1]:
        out.append("")
        out.append("  Most crowded precursors (nearest competitor, in sigma):")
        for uid, h in crowded:
            if not h:
                continue
            d, other, dd, dr = h[0]
            out.append("    %-34s %2d competitors   nearest %-22s drift %.2f, rt %.2f sigma"
                       % (uid[:34], len(h), other[:22], dd, dr))

    print("\n".join(out))
    return 1 if clean < n else 0


if __name__ == "__main__":
    sys.exit(main())
