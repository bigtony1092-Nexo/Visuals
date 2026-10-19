#!/usr/bin/env python3
"""Draw the panels an export can actually support, and leave the rest visibly blank.

The plan tells you which panels could render. This draws them. One page, inline
SVG, no scripts and no network: a viewer that needs a CDN is not a deliverable,
and a chart that silently drops the points it could not plot is worse than no
chart.

So the rules this follows, each of which exists because the opposite is easy:

  A blank panel is DRAWN, with the reason on it. A panel that renders nowhere
  and a panel that renders empty must not look alike.

  The drift axis is labelled from drift_kind. On SONAR it reads "quadrupole
  window position - NOT mobility", because an axis labelled "mobility" over a
  quadrupole position is a false figure, not a mislabelled one.

  Every point count is stated: drawn, and dropped with why. A log axis cannot
  show a zero, and a reader cannot see the absence.

  The tier is on every panel. Mixing a decoy-controlled tier with an untested
  one in unmarked points is how a confident wrong answer gets drawn.

  Quant is relative abundance unless an absolute standard was declared. No axis
  says concentration.

    python3 tools/render_panels.py export/ --method method.json -o panels.html

Renderability is decided by tools/plan_panels.py, not here: two sources of truth
about whether a panel may be drawn is how a blanked panel gets drawn anyway.
"""

import argparse
import html
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import plan_panels  # noqa: E402

W, H = 560, 360
PAD = {"l": 74, "r": 18, "t": 16, "b": 52}
TIER_MARK = {"atlas_consensus": ("#1b6ca8", "circle"),
             "pass2_confirmed": ("#2f8f4e", "square"),
             "bare_proximity": ("#b06000", "triangle")}
TIER_ORDER = ["atlas_consensus", "pass2_confirmed", "bare_proximity"]


def esc(s):
    return html.escape(str(s), quote=True)


def read_jsonl(path):
    rows = []
    if not Path(path).exists():
        return rows
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line:
            rows.append(json.loads(line))
    return rows


def dig(rec, dotted):
    cur = rec
    for key in dotted.split("."):
        if not isinstance(cur, dict):
            return None
        cur = cur.get(key)
    return cur


def nice_ticks(lo, hi, n=5):
    """Ticks on round numbers. A tick at 0.8333 tells the reader nothing."""
    if hi <= lo:
        return [lo], lo - 0.5, hi + 0.5
    span = hi - lo
    step = 10 ** math.floor(math.log10(span / max(n, 2)))
    for mult in (1, 2, 2.5, 5, 10):
        if span / (step * mult) <= n:
            step *= mult
            break
    first = math.floor(lo / step) * step
    ticks = []
    t = first
    while t <= hi + step * 1e-9:
        if t >= lo - step * 1e-9:
            ticks.append(t)
        t += step
    return ticks or [lo], min(lo, ticks[0] if ticks else lo), max(hi, ticks[-1] if ticks else hi)


def fmt(v):
    if v == 0:
        return "0"
    a = abs(v)
    if a >= 1e5 or a < 1e-3:
        return "%.0e" % v
    if a >= 100:
        return "%.0f" % v
    if a >= 1:
        return "%.1f" % v
    return "%.3f" % v


def marker(x, y, kind, colour, r=3.2):
    if kind == "square":
        return ('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" '
                'fill-opacity="0.78"/>' % (x - r, y - r, 2 * r, 2 * r, colour))
    if kind == "triangle":
        return ('<polygon points="%.1f,%.1f %.1f,%.1f %.1f,%.1f" fill="%s" '
                'fill-opacity="0.78"/>' % (x, y - r * 1.2, x - r * 1.1, y + r,
                                           x + r * 1.1, y + r, colour))
    return ('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" fill-opacity="0.78"/>'
            % (x, y, r, colour))


def scatter(points, xlab, ylab, ylog=False, xlog=False):
    """points: [(x, y, tier, label)]. Returns (svg, note_lines)."""
    notes = []
    usable, dropped_x, dropped_y = [], 0, 0
    for x, y, tier, label in points:
        if x is None or (isinstance(x, float) and math.isnan(x)):
            dropped_x += 1
            continue
        if y is None or (isinstance(y, float) and math.isnan(y)):
            dropped_y += 1
            continue
        if ylog and y <= 0:
            dropped_y += 1
            continue
        if xlog and x <= 0:
            dropped_x += 1
            continue
        usable.append((x, y, tier, label))
    # A log axis cannot show a zero and the reader cannot see what is missing,
    # so the count is stated rather than left to the gap in the cloud.
    if dropped_x:
        notes.append("%d point(s) had no %s%s and are not drawn"
                     % (dropped_x, xlab.split(" (")[0],
                        " or it was not positive, which a log axis cannot show"
                        if xlog else ""))
    if dropped_y:
        notes.append("%d point(s) had no %s%s and are not drawn"
                     % (dropped_y, ylab.split(" (")[0],
                        " or it was not positive, which a log axis cannot show"
                        if ylog else ""))
    if not usable:
        return None, notes + ["nothing left to draw"]

    xs = [p[0] for p in usable]
    ys = [p[1] for p in usable]
    tx = [math.log10(v) for v in xs] if xlog else xs
    ty = [math.log10(v) for v in ys] if ylog else ys
    xt, x0, x1 = nice_ticks(min(tx), max(tx))
    yt, y0, y1 = nice_ticks(min(ty), max(ty))
    if x1 == x0:
        x0, x1 = x0 - 1, x1 + 1
    if y1 == y0:
        y0, y1 = y0 - 1, y1 + 1
    pw = W - PAD["l"] - PAD["r"]
    ph = H - PAD["t"] - PAD["b"]

    def sx(v):
        return PAD["l"] + (v - x0) / (x1 - x0) * pw

    def sy(v):
        return PAD["t"] + ph - (v - y0) / (y1 - y0) * ph

    out = ['<svg viewBox="0 0 %d %d" width="100%%" role="img">' % (W, H)]
    out.append('<rect x="%d" y="%d" width="%.1f" height="%.1f" fill="#fcfcfa" '
               'stroke="#d8d8d2"/>' % (PAD["l"], PAD["t"], pw, ph))
    for t in yt:
        y = sy(t)
        out.append('<line x1="%d" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#ececE6"/>'
                   % (PAD["l"], y, PAD["l"] + pw, y))
        lab = ("1e%d" % round(t)) if ylog else fmt(t)
        out.append('<text x="%d" y="%.1f" font-size="10" text-anchor="end" '
                   'fill="#555">%s</text>' % (PAD["l"] - 6, y + 3, esc(lab)))
    for t in xt:
        x = sx(t)
        out.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%.1f" stroke="#ececE6"/>'
                   % (x, PAD["t"], x, PAD["t"] + ph))
        lab = ("1e%d" % round(t)) if xlog else fmt(t)
        out.append('<text x="%.1f" y="%.1f" font-size="10" text-anchor="middle" '
                   'fill="#555">%s</text>' % (x, PAD["t"] + ph + 15, esc(lab)))
    for x, y, tier, label in usable:
        colour, kind = TIER_MARK.get(tier, ("#777", "circle"))
        px = sx(math.log10(x) if xlog else x)
        py = sy(math.log10(y) if ylog else y)
        out.append('<g><title>%s</title>%s</g>'
                   % (esc(label), marker(px, py, kind, colour)))
    out.append('<text x="%.1f" y="%d" font-size="11" text-anchor="middle" '
               'fill="#222">%s</text>'
               % (PAD["l"] + pw / 2, H - 10, esc(xlab)))
    out.append('<text x="14" y="%.1f" font-size="11" text-anchor="middle" '
               'fill="#222" transform="rotate(-90 14 %.1f)">%s</text>'
               % (PAD["t"] + ph / 2, PAD["t"] + ph / 2, esc(ylab)))
    out.append('</svg>')
    notes.insert(0, "%d of %d point(s) drawn" % (len(usable), len(points)))
    return "\n".join(out), notes


# --------------------------------------------------------------------- panels

def axis_kind(drift_kinds):
    """Which single drift kind a drift axis may carry.

    Real mobility wins when it is present, because that is the kind the drift
    panels are about; everything else is then excluded rather than relabelled.
    """
    kinds = set(drift_kinds)
    for kind in ("ion_mobility", "quad_window"):
        if kind in kinds:
            return kind
    return "absent"


def drift_label(kind, mixed):
    """The axis label IS the claim. Getting it wrong makes the figure false."""
    base = {"ion_mobility": "drift (ion mobility)",
            "quad_window": "quadrupole window position - NOT mobility",
            "absent": "drift (no drift axis in this acquisition)"}[kind]
    return base + (" - other drift kinds excluded" if mixed else "")


def charge_str(coord):
    """Null charge is "z?" on a label, never 1+. The two are not the same claim."""
    z = coord.get("charge")
    return "z?" if z is None else "%d+" % z


def label_of(ion):
    ann = ion.get("annotation") or {}
    pep = ann.get("peptide_sequence") or "unassigned"
    return "%s  %s %s  rt %.2f  %s" % (
        ion["ion_uid"], pep, charge_str(ion["coord"]), ion["coord"]["rt_min"],
        ion["confidence"]["tier"])


def panel_intensity_vs_time(ions, ctx):
    pts = [(i["coord"]["rt_min"], dig(i, "intensity.matched_intensity"),
            i["confidence"]["tier"], label_of(i)) for i in ions]
    return scatter(pts, "retention time (min)",
                   "matched intensity - %s" % ctx["quant_words"], ylog=True)


def on_one_drift_axis(ions, ctx):
    """Keep only the ions whose drift means the same thing as the axis says.

    A mixed export is PARTIAL in the plan, which means the panel EXCLUDES the
    samples whose drift is something else. Excluding them is the behaviour; a
    "MIXED" caption over points from both kinds is not an exclusion, it is a
    false point with a disclaimer next to it.
    """
    if not ctx["mixed_drift"]:
        return ions, []
    keep = [i for i in ions if i["coord"]["drift_kind"] == ctx["axis_kind"]]
    dropped = len(ions) - len(keep)
    notes = []
    if dropped:
        notes.append("%d ion(s) excluded: their drift_kind is not %s, so their"
                     " drift is a different physical quantity and cannot share"
                     " this axis" % (dropped, ctx["axis_kind"]))
    return keep, notes


def panel_drift_vs_time(ions, ctx):
    ions, extra = on_one_drift_axis(ions, ctx)
    pts = [(i["coord"]["rt_min"], i["coord"].get("drift"),
            i["confidence"]["tier"], label_of(i)) for i in ions]
    svg, notes = scatter(pts, "retention time (min)", ctx["drift_label"])
    return svg, (extra + notes if svg is None else notes + extra)


def panel_mz_vs_drift(ions, ctx):
    ions, extra = on_one_drift_axis(ions, ctx)
    pts = [(i["coord"]["mz"], i["coord"].get("drift"),
            i["confidence"]["tier"], label_of(i)) for i in ions]
    svg, notes = scatter(pts, "m/z", ctx["drift_label"])
    return svg, (extra + notes if svg is None else notes + extra)


def panel_seqstart_vs_intensity(ions, ctx):
    # The panel's own rule: primary cluster only. Pooling clusters was found to
    # corrupt the pattern, so a missing cluster_rank is excluded, not assumed 1.
    keep, no_rank = [], 0
    for i in ions:
        rank = dig(i, "intensity.cluster_rank")
        if rank is None:
            no_rank += 1
            continue
        if rank != 1:
            continue
        keep.append(i)
    pts = [(dig(i, "annotation.seq_start"),
            dig(i, "intensity.summed_precursor_intensity")
            or dig(i, "intensity.matched_intensity"),
            i["confidence"]["tier"], label_of(i)) for i in keep]
    svg, notes = scatter(pts, "seq_start (residue position)",
                         "summed precursor intensity - %s" % ctx["quant_words"],
                         ylog=True)
    if no_rank:
        why = ("%d ion(s) carry no cluster_rank and are excluded: this panel is"
               " primary-cluster only, and pooling clusters corrupts the pattern"
               % no_rank)
        notes.insert(0, why) if svg is None else notes.append(why)
    return svg, notes


def panel_forensic_track(ions, ctx):
    """The event track as a table. Replaying it in order is the flight."""
    rows = []
    for i in ions[:40]:
        for n, e in enumerate(i["events"]):
            snap = e.get("coord_snapshot") or {}
            rows.append("<tr><td class=m>%s</td><td>%s</td><td>%s</td>"
                        "<td class=m>%s</td><td class=m>%s</td><td>%s</td></tr>"
                        % (esc(i["ion_uid"]) if n == 0 else "",
                           esc(e["stage"]), esc(e["at"]),
                           fmt(snap["mz"]) if snap.get("mz") is not None else "",
                           fmt(snap["rt_min"]) if snap.get("rt_min") is not None else "",
                           esc(e.get("note") or "")))
    if not rows:
        return None, ["no events on any ion"]
    note = ["%d event(s) over %d ion(s)"
            % (sum(len(i["events"]) for i in ions[:40]), min(len(ions), 40))]
    if len(ions) > 40:
        note.append("showing the first 40 of %d ions" % len(ions))
    return ("<table class=track><thead><tr><th>ion_uid</th><th>stage</th>"
            "<th>at</th><th>m/z</th><th>rt</th><th>note</th></tr></thead>"
            "<tbody>%s</tbody></table>" % "".join(rows)), note


def panel_ion_flight(ions, ctx):
    """The flight is the events replayed. With one snapshot there is no path."""
    paths, moved = [], 0
    for i in ions:
        snaps = [e["coord_snapshot"] for e in i["events"]
                 if e.get("coord_snapshot")
                 and e["coord_snapshot"].get("mz") is not None
                 and e["coord_snapshot"].get("rt_min") is not None]
        if len(snaps) >= 2:
            moved += 1
        paths.append((i, snaps))
    if moved == 0:
        return None, ["every ion carries one coordinate snapshot, so there is no"
                      " path to fly. A second snapshot appears when a stage MOVES"
                      " an ion -- recalibration is the usual one -- and this panel"
                      " is what makes that move visible instead of silent."]
    pts = []
    for i, snaps in paths:
        for s in snaps:
            pts.append((s["rt_min"], s["mz"], i["confidence"]["tier"], label_of(i)))
    svg, notes = scatter(pts, "retention time (min)", "m/z")
    notes.append("%d of %d ion(s) moved between stages" % (moved, len(ions)))
    return svg, notes


C13_SPACING = 1.003355

CSS = """
:root{--fg:#1a1a18;--mut:#666;--line:#ddd;--card:#fff;--bg:#f4f3ef;
--warn:#8a5a00;--warnbg:#fff8e8;--blank:#f7f7f5}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
.wrap{max-width:1180px;margin:0 auto;padding:24px 16px 56px}
h1{font-size:23px;margin:0 0 4px;letter-spacing:-.01em}
.sub{color:var(--mut);font-size:13px;margin-bottom:18px}
.meta{display:flex;flex-wrap:wrap;gap:8px;margin:0 0 18px}
.chip{background:var(--card);border:1px solid var(--line);border-radius:999px;
padding:3px 11px;font-size:12px}
.banner{background:var(--warnbg);border:1px solid #e8d49a;border-left:4px solid var(--warn);
border-radius:6px;padding:11px 14px;margin:0 0 18px;font-size:13px;color:#4a3500}
.banner b{color:var(--warn)}
.banner ul{margin:6px 0 0;padding-left:20px}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,460px),1fr));gap:16px}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;
padding:14px 16px 12px;min-width:0}
.card.wide{grid-column:1/-1}
.card h2{font-size:15px;margin:0 0 2px;font-weight:600}
.fam{font-size:11px;color:var(--mut);text-transform:uppercase;letter-spacing:.07em}
.blankbox{background:var(--blank);border:1px dashed #c4c4bc;border-radius:6px;
padding:26px 18px;text-align:center;color:var(--mut);font-size:13px;margin-top:10px}
.blankbox .why{display:block;margin-top:7px;color:#8a3a00;font-size:12.5px;text-align:left}
.notes{margin:8px 0 0;padding-left:18px;color:var(--mut);font-size:12px}
.notes li{margin:1px 0}
.legend{display:flex;flex-wrap:wrap;gap:14px;font-size:12px;color:var(--mut);
margin:14px 0 0;padding-top:10px;border-top:1px solid var(--line)}
.legend span{display:inline-flex;align-items:center;gap:5px}
.sw{width:10px;height:10px;display:inline-block}
table.track{width:100%;border-collapse:collapse;font-size:11.5px;margin-top:8px}
table.track th{text-align:left;border-bottom:1px solid var(--line);padding:4px 6px;
color:var(--mut);font-weight:600}
table.track td{padding:3px 6px;border-bottom:1px solid #f0f0ec;vertical-align:top}
td.m,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
footer{margin-top:30px;color:var(--mut);font-size:12px;border-top:1px solid var(--line);
padding-top:12px}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
--fg:#e8e8e4;--mut:#9b9b95;--line:#333;--card:#1c1c1a;--bg:#131312;
--warnbg:#2a2110;--blank:#1a1a18}
:root:not([data-theme="light"]) .banner{color:#e8d8b0;border-color:#5a4a20}
:root:not([data-theme="light"]) .blankbox{border-color:#3a3a36}
:root:not([data-theme="light"]) .blankbox .why{color:#e0a060}
:root:not([data-theme="light"]) svg rect[fill="#fcfcfa"]{fill:#242422;stroke:#3a3a36}
:root:not([data-theme="light"]) svg line{stroke:#2e2e2b}
:root:not([data-theme="light"]) svg text{fill:#bbb}}
:root[data-theme="dark"]{--fg:#e8e8e4;--mut:#9b9b95;--line:#333;--card:#1c1c1a;
--bg:#131312;--warnbg:#2a2110;--blank:#1a1a18}
"""


ION_COLOUR = {"b": "#1b6ca8", "y": "#b3312c", None: "#9a9a92"}


def series_of(peak):
    """b, y, or unassigned. The colour convention is the one practitioners read."""
    ion = peak.get("ion")
    if not ion:
        return None
    first = str(ion)[0].lower()
    return first if first in ("b", "y") else None


def bar_plot(series, xlab, ylab, mirrored=False, annotate=True):
    """Bars on a shared m/z axis. mirrored draws series[1] downward.

    series: [(name, [(mz, rel, label, colour_key, tip)], ...)]
    """
    allpts = [p for _, pts in series for p in pts]
    if not allpts:
        return None, ["no peaks to draw"]
    xs = [p[0] for p in allpts]
    hi = max(abs(p[1]) for p in allpts) or 1.0
    xt, x0, x1 = nice_ticks(min(xs), max(xs))
    span = (x1 - x0) or 1.0
    x0 -= span * 0.04
    x1 += span * 0.04
    pw = W - PAD["l"] - PAD["r"]
    ph = H - PAD["t"] - PAD["b"]
    mid = PAD["t"] + (ph / 2 if mirrored else ph)

    def sx(v):
        return PAD["l"] + (v - x0) / (x1 - x0) * pw

    def sy(rel, down=False):
        reach = (ph / 2 if mirrored else ph) * 0.88
        return mid + (reach * rel / hi if down else -reach * rel / hi)

    out = ['<svg viewBox="0 0 %d %d" width="100%%" role="img">' % (W, H)]
    out.append('<rect x="%d" y="%d" width="%.1f" height="%.1f" fill="#fcfcfa" '
               'stroke="#d8d8d2"/>' % (PAD["l"], PAD["t"], pw, ph))
    for frac in (0.25, 0.5, 0.75, 1.0):
        for down in ((False, True) if mirrored else (False,)):
            y = sy(hi * frac, down)
            out.append('<line x1="%d" y1="%.1f" x2="%.1f" y2="%.1f" '
                       'stroke="#ececE6"/>' % (PAD["l"], y, PAD["l"] + pw, y))
            out.append('<text x="%d" y="%.1f" font-size="9" text-anchor="end" '
                       'fill="#777">%d</text>'
                       % (PAD["l"] - 5, y + 3, round(hi * frac)))
    out.append('<line x1="%d" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#999"/>'
               % (PAD["l"], mid, PAD["l"] + pw, mid))
    for t in xt:
        x = sx(t)
        out.append('<text x="%.1f" y="%.1f" font-size="10" text-anchor="middle" '
                   'fill="#555">%s</text>' % (x, PAD["t"] + ph + 15, fmt(t)))
    for n, (name, pts) in enumerate(series):
        down = mirrored and n > 0
        for mz, rel, label, key, tip in pts:
            x = sx(mz)
            y = sy(rel, down)
            colour = ION_COLOUR.get(key, "#9a9a92")
            out.append('<g><title>%s</title><line x1="%.1f" y1="%.1f" x2="%.1f" '
                       'y2="%.1f" stroke="%s" stroke-width="1.7"/></g>'
                       % (esc(tip), x, mid, x, y, colour))
            if annotate and label and abs(rel) >= hi * 0.18:
                out.append('<text x="%.1f" y="%.1f" font-size="9" '
                           'text-anchor="middle" fill="%s">%s</text>'
                           % (x, y + (10 if down else -4), colour, esc(label)))
        if mirrored:
            out.append('<text x="%d" y="%.1f" font-size="10" fill="#666">%s</text>'
                       % (PAD["l"] + 6, mid + (16 if down else -8), esc(name)))
    out.append('<text x="%.1f" y="%d" font-size="11" text-anchor="middle" '
               'fill="#222">%s</text>' % (PAD["l"] + pw / 2, H - 10, esc(xlab)))
    out.append('<text x="14" y="%.1f" font-size="11" text-anchor="middle" '
               'fill="#222" transform="rotate(-90 14 %.1f)">%s</text>'
               % (PAD["t"] + ph / 2, PAD["t"] + ph / 2, esc(ylab)))
    out.append('</svg>')
    return "\n".join(out), []


def line_plot(curves, xlab, ylab, marks=()):
    """curves: [(name, [(x, y)], colour)]. marks: [(x, label)] vertical rules."""
    allpts = [p for _, pts, _ in curves for p in pts]
    if not allpts:
        return None, ["no points to draw"]
    xt, x0, x1 = nice_ticks(min(p[0] for p in allpts), max(p[0] for p in allpts))
    yt, y0, y1 = nice_ticks(0.0, max(p[1] for p in allpts))
    pw = W - PAD["l"] - PAD["r"]
    ph = H - PAD["t"] - PAD["b"]
    if x1 == x0:
        x1 = x0 + 1
    if y1 == y0:
        y1 = y0 + 1

    def sx(v):
        return PAD["l"] + (v - x0) / (x1 - x0) * pw

    def sy(v):
        return PAD["t"] + ph - (v - y0) / (y1 - y0) * ph

    out = ['<svg viewBox="0 0 %d %d" width="100%%" role="img">' % (W, H)]
    out.append('<rect x="%d" y="%d" width="%.1f" height="%.1f" fill="#fcfcfa" '
               'stroke="#d8d8d2"/>' % (PAD["l"], PAD["t"], pw, ph))
    for t in yt:
        y = sy(t)
        out.append('<line x1="%d" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#ececE6"/>'
                   % (PAD["l"], y, PAD["l"] + pw, y))
        out.append('<text x="%d" y="%.1f" font-size="10" text-anchor="end" '
                   'fill="#555">%s</text>' % (PAD["l"] - 6, y + 3, fmt(t)))
    for t in xt:
        x = sx(t)
        out.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%.1f" stroke="#ececE6"/>'
                   % (x, PAD["t"], x, PAD["t"] + ph))
        out.append('<text x="%.1f" y="%.1f" font-size="10" text-anchor="middle" '
                   'fill="#555">%s</text>' % (x, PAD["t"] + ph + 15, fmt(t)))
    for x, label in marks:
        px = sx(x)
        out.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%.1f" stroke="#b3312c" '
                   'stroke-dasharray="3 3"/>' % (px, PAD["t"], px, PAD["t"] + ph))
        out.append('<text x="%.1f" y="%d" font-size="9" text-anchor="middle" '
                   'fill="#b3312c">%s</text>' % (px, PAD["t"] - 4, esc(label)))
    for name, pts, colour in curves:
        d = " ".join("%s%.1f,%.1f" % ("M" if n == 0 else "L", sx(x), sy(y))
                     for n, (x, y) in enumerate(sorted(pts)))
        out.append('<g><title>%s</title><path d="%s" fill="none" stroke="%s" '
                   'stroke-width="1.6"/></g>' % (esc(name), d, colour))
        for x, y in sorted(pts):
            out.append('<circle cx="%.1f" cy="%.1f" r="2" fill="%s"/>'
                       % (sx(x), sy(y), colour))
    out.append('<text x="%.1f" y="%d" font-size="11" text-anchor="middle" '
               'fill="#222">%s</text>' % (PAD["l"] + pw / 2, H - 10, esc(xlab)))
    out.append('<text x="14" y="%.1f" font-size="11" text-anchor="middle" '
               'fill="#222" transform="rotate(-90 14 %.1f)">%s</text>'
               % (PAD["t"] + ph / 2, PAD["t"] + ph / 2, esc(ylab)))
    out.append('</svg>')
    return "\n".join(out), []


def by_uid(rows, kind=None):
    out = {}
    for r in rows:
        if kind and r.get("kind") != kind:
            continue
        out.setdefault(r["ion_uid"], []).append(r)
    return out


def pick_ion(ions, ctx, need):
    """Which ion to show on a single-ion panel, and why that one.

    The most intense ion at the best tier WITH the required table, so the panel
    is never blank because an arbitrary first row happened to lack a spectrum.
    """
    have = [i for i in ions if i["ion_uid"] in need]
    if not have:
        return None
    have.sort(key=lambda i: (-TIER_ORDER[::-1].index(i["confidence"]["tier"]),
                             -(dig(i, "intensity.matched_intensity") or 0)))
    return have[0]


def panel_msms_mirror(ions, ctx):
    """The load-bearing panel: experimental against the Atlas library.

    If a geometry-selected identification is right, this is where it shows in a
    form nobody has to be taught. So it is also where a wrong one shows.
    """
    exp = by_uid(ctx["spectra"], "ms2_experimental")
    lib = by_uid(ctx["spectra"], "ms2_library")
    both = set(exp) & set(lib)
    if not both:
        why = []
        if exp and not lib:
            why.append("%d experimental spectrum(a) but no ms2_library: there is"
                       " nothing to mirror against, and a one-sided plot is not a"
                       " comparison" % len(exp))
        elif lib and not exp:
            why.append("library spectra only, no ms2_experimental")
        else:
            why.append("no ion carries both an experimental and a library spectrum")
        return None, why
    ion = pick_ion(ions, ctx, both)
    if ion is None:
        return None, ["spectra are present but none of their ion_uids appears in"
                      " ions.jsonl, so no spectrum can be attributed to an ion"]
    uid = ion["ion_uid"]
    e, l = exp[uid][0], lib[uid][0]

    def pts(spec, side):
        out = []
        for pk in spec["peaks"]:
            key = series_of(pk)
            ppm = pk.get("ppm_error")
            tip = "%s %.4f  rel %.1f%s" % (pk.get("ion") or "unassigned", pk["mz"],
                                           pk["intensity"],
                                           "  %+.1f ppm" % ppm if ppm is not None else "")
            out.append((pk["mz"], pk["intensity"], pk.get("ion"), key, side + ": " + tip))
        return out

    svg, _ = bar_plot([("experimental", pts(e, "experimental")),
                       ("Atlas library", pts(l, "library"))],
                      "m/z", "relative intensity (base peak = 100)", mirrored=True)
    m = e.get("match") or {}
    ann = ion["annotation"] or {}
    notes = ["%s %s  %s  %s" % (ann.get("peptide_sequence") or "unassigned",
                                charge_str(ion["coord"]),
                                 ann.get("prot_accession") or "",
                                 ion["confidence"]["tier"])]
    if m.get("dot_product") is not None:
        notes.append("dot product %.3f on %s, %s of %s library peaks matched"
                     % (m["dot_product"], m.get("basis") or "unstated basis",
                        m.get("n_matched"), m.get("n_library")))
        if m.get("basis") != "sqrt_intensity":
            notes.append("the dot product is not on sqrt intensities, so it is"
                          " dominated by the base peak and is not the number the"
                          " panel's convention describes")
    else:
        notes.append("no dot product on this spectrum: the mirror is drawn but"
                     " nothing scores the agreement")
    ppms = [pk["ppm_error"] for pk in e["peaks"] if pk.get("ppm_error") is not None]
    if ppms:
        notes.append("matched-peak error %+.1f to %+.1f ppm" % (min(ppms), max(ppms)))
    unassigned = [pk for pk in e["peaks"] if not pk.get("ion")]
    if unassigned:
        tot = sum(pk["intensity"] for pk in e["peaks"]) or 1.0
        notes.append("%d experimental peak(s), %.0f%% of the intensity, match no"
                     " ladder ion and are drawn grey"
                     % (len(unassigned),
                        100.0 * sum(pk["intensity"] for pk in unassigned) / tot))
    missed = [pk.get("ion") for pk in l["peaks"]
              if pk.get("ion") not in {p.get("ion") for p in e["peaks"]}]
    if missed:
        notes.append("library ion(s) absent from the experimental spectrum: %s"
                     % ", ".join(str(x) for x in missed if x))
    if m.get("library_source"):
        notes.append("library: %s" % m["library_source"])
    if len(both) > 1:
        notes.append("%d ion(s) have both spectra; this is the most intense at the"
                     " best tier" % len(both))
    return svg, notes


def panel_isotope_envelope(ions, ctx):
    """Observed envelope with the expected pattern over it.

    A modelled (averagine) envelope must never be drawn as if it were measured,
    so a modelled peak is hatched and counted in the notes.
    """
    iso = by_uid(ctx["spectra"], "ms1_isotope")
    with_env = {i["ion_uid"] for i in ions if i.get("isotopes")}
    usable = set(iso) | with_env
    if not usable:
        return None, ["no ms1_isotope spectrum and no ion carries a resolved"
                      " isotope envelope. A consumer would have to model one, and"
                      " a modelled envelope is not evidence of a charge state."]
    ion = pick_ion(ions, ctx, usable)
    if ion is None:
        return None, ["isotope data present but attributable to no ion in ions.jsonl"]
    uid = ion["ion_uid"]
    z = ion["coord"]["charge"]
    if z is None:
        return None, ["this ion's charge was not determined, and the expected"
                      " spacing is 1.00335/z. Without z there is no spacing to"
                      " check the envelope against -- and if the envelope were"
                      " resolved, its spacing would have given the charge."]
    env = ion.get("isotopes") or []
    if not env and uid in iso:
        env = [{"k": pk.get("idx", n), "mz": pk["mz"], "rel": pk["intensity"],
                "observed": True} for n, pk in enumerate(iso[uid][0]["peaks"])]
    modelled = [e for e in env if e.get("observed") is False]
    obs = [(e["mz"], e["rel"], "k=%d" % e["k"], None,
            "k=%d  %.4f  rel %.1f%s" % (e["k"], e["mz"], e["rel"],
                                        "  MODELLED, not measured"
                                        if e.get("observed") is False else ""))
           for e in env]
    series = [("observed", obs)]
    fit = ion.get("isotope_fit") or {}
    expected = []
    for e in env:
        delta = next((r["delta"] for r in (fit.get("residual") or [])
                      if r["k"] == e["k"]), None)
        if delta is not None:
            expected.append((e["mz"], e["rel"] - delta, None, "y",
                             "expected k=%d  rel %.1f" % (e["k"], e["rel"] - delta)))
    if expected:
        series.append(("expected (%s)" % fit.get("expected_basis", "?"), expected))
    svg, _ = bar_plot(series, "m/z", "relative intensity (most intense = 100)",
                      mirrored=bool(expected))
    spacing = C13_SPACING / z
    notes = ["charge %d+: isotope spacing should be 1.00335/%d = %.4f Th"
             % (z, z, spacing)]
    if len(env) >= 2:
        gaps = [env[n + 1]["mz"] - env[n]["mz"] for n in range(len(env) - 1)]
        notes.append("observed spacing %s Th"
                     % ", ".join("%.4f" % g for g in gaps))
        worst = max(abs(g - spacing) for g in gaps)
        if worst > 0.01:
            notes.append("the widest gap is %.4f Th off the expected spacing, which"
                         " is more than a charge state apart at this m/z -- the"
                         " envelope is not one species" % worst)
    if fit.get("rms") is not None:
        notes.append("residual RMS %.2f percentage points, basis %s"
                     % (fit["rms"], fit.get("expected_basis")))
    if not expected:
        notes.append("no expected envelope on this ion, so the observed one is"
                     " drawn alone and nothing here tests it")
    for r in (fit.get("ratio_cv") or []):
        if r["cv"] > 0.25:
            notes.append("isotope ratio k=%d varies %.0f%% across the elution peak:"
                         " something enters or leaves during it"
                         % (r["k"], 100 * r["cv"]))
    if modelled:
        notes.append("%d of %d isotope(s) are MODELLED (averagine), not measured,"
                     " and must not be read as observation" % (len(modelled), len(env)))
    return svg, notes


def panel_xic(ions, ctx):
    traces = by_uid(ctx["traces"], "xic")
    if not traces:
        return None, ["no xic traces in the export"]
    ion = pick_ion(ions, ctx, set(traces))
    if ion is None:
        return None, ["traces present but attributable to no ion in ions.jsonl"]
    uid = ion["ion_uid"]
    rows = traces[uid]
    curves, marks, notes = [], [], []
    palette = ["#1b6ca8", "#2f8f4e", "#b06000", "#8a3a8a"]
    for n, t in enumerate(rows[:4]):
        pts = [(p[0], p[1]) for p in t["points"] if len(p) >= 2]
        if len(pts) < 2:
            notes.append("a trace carries %d point(s): a peak needs at least two"
                         " consecutive scans to be a peak" % len(pts))
            continue
        curves.append(("%s window %.1f ppm" % (t["sample_id"],
                                               t.get("extraction_window_ppm") or 0),
                       pts, palette[n % len(palette)]))
        if t.get("apex_rt_min") is not None and n == 0:
            marks.append((t["apex_rt_min"], "apex %.3f" % t["apex_rt_min"]))
    if not curves:
        return None, notes or ["no trace had enough points to draw"]
    svg, _ = line_plot(curves, "retention time (min)",
                       "intensity - %s" % ctx["quant_words"], marks)
    ann = ion["annotation"] or {}
    notes.insert(0, "%s %s at m/z %.4f, %s"
                 % (ann.get("peptide_sequence") or "unassigned",
                    charge_str(ion["coord"]), ion["coord"]["mz"],
                    ion["confidence"]["tier"]))
    windows = {t.get("extraction_window_ppm") for t in rows[:4]}
    stated = sorted(w for w in windows if isinstance(w, (int, float)) and w > 0)
    if stated:
        notes.append("extraction window%s %s ppm -- a wider window admits more"
                     " interference, so the window is part of the peak's shape"
                     % ("s" if len(stated) > 1 else "",
                        ", ".join("%g" % w for w in stated)))
    if len(stated) < len(windows):
        notes.append("%d trace(s) state no extraction window (or state zero, which"
                     " is not a window). Without it there is no telling how much"
                     " of this peak is the ion and how much is whatever else fell"
                     " inside the extraction" % (len(windows) - len(stated)))
    shapes = {}
    for t in rows[:4]:
        if t.get("fwhm_min") is not None:
            key = (t["sample_id"], round(t["fwhm_min"], 4), len(t["points"]))
            shapes[key] = shapes.get(key, 0) + 1
    for (sid, fwhm, n), count in sorted(shapes.items()):
        notes.append("%s FWHM %.3f min over %d point(s)%s"
                     % (sid, fwhm, n,
                        " (%d traces identical on these three)" % count
                        if count > 1 else ""))
    if len(set(t["sample_id"] for t in rows)) == 1 and len(rows) > 1:
        notes.append("all %d trace(s) carry the same sample_id and ion_uid, so they"
                     " are repeats of one extraction, not a comparison across"
                     " samples" % len(rows))
    scans = dig(ion, "intensity.n_consecutive_scans")
    if scans is not None:
        notes.append("the ion claims %d consecutive scans; this trace shows %d"
                     " point(s)" % (scans, len(rows[0]["points"])))
    if len(rows) > 4:
        notes.append("%d traces on this ion, showing 4" % len(rows))
    return svg, notes


def cross_sample(rows, ctx, need, what):
    """A cross-sample panel needs several samples IN ITS OWN TABLE.

    plan_panels counts the manifest's sample list, which is what the export
    describes. A table can still hold one sample's rows, and a cross-sample
    figure drawn from one sample is false however many samples the manifest
    declares.
    """
    seen = sorted({r["sample_id"] for r in rows})
    if len(seen) < need:
        return None, ["%s covers %d sample(s) (%s) but this panel compares across"
                      " %d. The manifest declares %d samples; the table does not"
                      " hold them, and a cross-sample figure drawn from %d is"
                      " false." % (what, len(seen), ", ".join(seen) or "none",
                                   need, len(ctx["samples"]), len(seen))]
    return seen, []


def panel_top3_percentile(ions, ctx):
    rows = ctx["proteins"]
    if not rows:
        return None, ["no proteins.jsonl in the export"]
    seen, why = cross_sample(rows, ctx, 2, "proteins.jsonl")
    if seen is None:
        return None, why
    pts = [(r.get("top3avg"), r.get("percentile"), r.get("tier", "bare_proximity"),
            "%s  %s  top3 %.3f  pct %.1f" % (r["sample_id"], r["prot_accession"],
                                             r.get("top3avg") or 0,
                                             r.get("percentile") or 0))
           for r in rows]
    svg, notes = scatter(pts, "Top3avg - %s" % ctx["quant_words"],
                         "percentile within sample", xlog=True)
    notes.append("%d sample(s): %s" % (len(seen), ", ".join(seen)))
    uncontrolled = [r for r in rows if not r.get("decoy_controlled")]
    if uncontrolled:
        notes.append("%d of %d protein row(s) are not decoy controlled, so their"
                     " rank and percentile are positions in a list that includes"
                     " untested identifications" % (len(uncontrolled), len(rows)))
    return svg, notes


def panel_stoichiometry(ions, ctx):
    """Ratios within ONE run only, carrying n_peptides_used on both sides.

    The panel's own rule: a 1-peptide against a 3-peptide ratio is not the same
    measurement, so the count rides with the number instead of being dropped.
    """
    rows = ctx["proteins"]
    if not rows:
        return None, ["no proteins.jsonl in the export"]
    by_sample = {}
    for r in rows:
        by_sample.setdefault(r["sample_id"], []).append(r)
    runs = [(sid, rs) for sid, rs in sorted(by_sample.items()) if len(rs) >= 2]
    if not runs:
        return None, ["no run holds two proteins, and a ratio needs two. Rows per"
                      " run: %s. Ratios across runs are not drawn here: that mixes"
                      " loading into the ratio."
                      % ", ".join("%s=%d" % (s, len(r))
                                  for s, r in sorted(by_sample.items()))]
    pts, notes, thin = [], [], 0
    for sid, rs in runs:
        rs = sorted(rs, key=lambda r: -(r.get("top3avg") or 0))
        base = rs[0]
        for other in rs[1:]:
            if not base.get("top3avg") or not other.get("top3avg"):
                continue
            ratio = other["top3avg"] / base["top3avg"]
            n_lo = min(base.get("n_peptides_used") or 0,
                       other.get("n_peptides_used") or 0)
            if n_lo < 2:
                thin += 1
            pts.append((n_lo, ratio, other.get("tier", "bare_proximity"),
                        "%s  %s / %s = %.3f  (n_peptides %s vs %s)"
                        % (sid, other["prot_accession"], base["prot_accession"],
                           ratio, other.get("n_peptides_used"),
                           base.get("n_peptides_used"))))
    svg, snotes = scatter(pts, "peptides behind the weaker side of the ratio",
                          "ratio to the strongest protein in the run", ylog=True)
    notes = snotes + ["%d ratio(s) from %d run(s): %s"
                      % (len(pts), len(runs), ", ".join(s for s, _ in runs))]
    if thin:
        notes.append("%d ratio(s) rest on a single peptide on one side. A 1-peptide"
                     " against a 3-peptide ratio is not the same measurement, which"
                     " is why the count is the x axis rather than a footnote" % thin)
    return svg, notes


def panel_flux(ions, ctx):
    """Change across runs, on a STATED normalization.

    Raw Top3avg across runs mixes biology with how much was loaded, so a run
    without a declared normalization is not drawn.
    """
    rows = ctx["proteins"]
    if not rows:
        return None, ["no proteins.jsonl in the export"]
    norms = (ctx["manifest_quant"].get("normalizations") or [])
    if not norms:
        return None, ["the manifest declares no normalization. Raw Top3avg across"
                      " runs mixes biology with loading, so there is no quantity"
                      " here that a flux panel could honestly plot."]
    field = next((n for n in norms
                  if any((r.get("normalized") or {}).get(n) is not None
                         for r in rows)), None)
    if field is None:
        return None, ["the manifest declares %s but no protein row carries any of"
                      " them under `normalized`, so the stated normalization was"
                      " not actually applied" % ", ".join(norms)]
    seen, why = cross_sample(rows, ctx, 2, "proteins.jsonl")
    if seen is None:
        return None, why + ["normalization on offer: %s" % field]
    order = {sid: n for n, sid in enumerate(seen)}
    curves, palette = [], ["#1b6ca8", "#2f8f4e", "#b06000", "#8a3a8a"]
    by_prot = {}
    for r in rows:
        v = (r.get("normalized") or {}).get(field)
        if v is not None:
            by_prot.setdefault(r["prot_accession"], []).append(
                (order[r["sample_id"]], v))
    for n, (prot, pts) in enumerate(sorted(by_prot.items())):
        if len(pts) >= 2:
            curves.append((prot, pts, palette[n % len(palette)]))
    if not curves:
        return None, ["no protein appears in two runs under %s, so nothing has a"
                      " flux to draw" % field]
    svg, notes = line_plot(curves, "run (%s)" % " - ".join(seen), field)
    notes.append("normalization %s, as declared in the manifest" % field)
    notes.append("%d protein(s) seen in 2+ runs of %d in the table"
                 % (len(curves), len(by_prot)))
    return svg, notes


def panel_fragment_hits(ions, ctx):
    """Fragment-first candidates: what was queried, matched, and controlled."""
    rows = ctx["fragment_hits"]
    if not rows:
        return None, ["no fragment_hits.jsonl in the export"]
    body = ["<table class=track><thead><tr><th>hit</th><th>peptide</th>"
            "<th>k of n</th><th>coelution</th><th>delta mass</th>"
            "<th>controls</th><th>evidence</th></tr></thead><tbody>"]
    notes, uncontrolled = [], 0
    for r in rows[:30]:
        k = r.get("k_of_n") or {}
        co = r.get("coelution") or {}
        dm = r.get("delta_mass") or {}
        ctl = r.get("control") or {}
        if ctl.get("reversed_matched") is None and ctl.get("entrapment_matched") is None:
            uncontrolled += 1
        cand = ", ".join("%s %.4f" % (c.get("name"), c.get("monoisotopic") or 0)
                         for c in (dm.get("candidates") or [])[:2])
        body.append("<tr><td class=m>%s</td><td class=m>%s %s</td>"
                    "<td>%s/%s (need %s)</td><td>%s</td><td>%s</td>"
                    "<td>rev %s, entrap %s</td><td>%s</td></tr>"
                    % (esc(r.get("hit_id")),
                       esc((r.get("peptide") or {}).get("sequence")),
                       esc((r.get("peptide") or {}).get("prot_accession") or ""),
                       k.get("matched"), k.get("queried"), k.get("required"),
                       ("%.3f min, %d scans" % (co.get("rt_span_min") or 0,
                                                co.get("consecutive_scans") or 0)),
                       ("%+.5f (%+.1f ppm)%s" % (dm.get("value") or 0,
                                                 dm.get("ppm") or 0,
                                                 "  " + cand if cand else "")
                        if dm.get("value") is not None else ""),
                       ctl.get("reversed_matched"), ctl.get("entrapment_matched"),
                       esc(r.get("evidence") or "")))
    body.append("</tbody></table>")
    notes.append("%d hit(s)%s" % (len(rows),
                                  ", showing 30" if len(rows) > 30 else ""))
    if uncontrolled:
        notes.append("%d hit(s) carry no reversed or entrapment count. k-of-n on an"
                     " uncontrolled query is a count, not a false-discovery rate."
                     % uncontrolled)
    offsets = {(r.get("coelution") or {}).get("drift_offset_bins_applied")
               for r in rows}
    offsets = {o for o in offsets if o}
    if offsets:
        notes.append("a drift offset was applied to the coelution test (%s bins):"
                     " the fragments were moved to meet the precursor, so the"
                     " agreement is partly constructed"
                     % ", ".join(str(o) for o in sorted(offsets)))
    return "".join(body), notes


DRAW = {"intensity_vs_time": panel_intensity_vs_time,
        "msms_mirror": panel_msms_mirror,
        "isotope_envelope": panel_isotope_envelope,
        "xic": panel_xic,
        "top3_percentile": panel_top3_percentile,
        "fragment_hits": panel_fragment_hits,
        "stoichiometry": panel_stoichiometry,
        "flux": panel_flux,
        "drift_vs_time": panel_drift_vs_time,
        "mz_vs_drift": panel_mz_vs_drift,
        "seqstart_vs_intensity": panel_seqstart_vs_intensity,
        "forensic_track": panel_forensic_track,
        "ion_flight": panel_ion_flight}


def build(manifest, method, ions, catalog, selected, tables=None):
    tables = tables or {}
    samples = manifest.get("samples", [])
    kinds = [s["drift_kind"] for s in samples]
    quant = manifest.get("quant") or {}
    absolute = quant.get("absolute_standard_declared") or \
        (method.get("quant", {}) or {}).get("absolute_standard")
    one_kind = axis_kind(kinds)
    ctx = {"manifest_quant": quant,
           "spectra": tables.get("spectra", []),
           "traces": tables.get("traces", []),
           "proteins": tables.get("proteins", []),
           "fragment_hits": tables.get("fragment_hits", []),
           "samples": samples,
           "axis_kind": one_kind,
           "mixed_drift": len(set(kinds)) > 1,
           "drift_label": drift_label(one_kind, len(set(kinds)) > 1),
           "quant_words": "absolute (standard declared)" if absolute
                          else "RELATIVE abundance, not concentration"}

    tiers_present = [t for t in TIER_ORDER
                     if any(i["confidence"]["tier"] == t for i in ions)]
    head = []
    head.append("<h1>%s</h1>" % esc(method.get("name") or manifest["export_id"]))
    head.append('<div class=sub>%s &middot; %s &middot; drawn by tools/render_panels.py'
                ' from %d ion record(s)</div>'
                % (esc(manifest["export_id"]), esc(manifest.get("created_at", "")),
                   len(ions)))
    chips = ["%d sample(s)" % len(samples),
             "%s ions exported of %s" % (f"{manifest['counts']['ions_exported']:,}",
                                         f"{manifest['counts']['ions_total']:,}")]
    for s in samples[:4]:
        chips.append("%s: %s / %s / %s" % (s["sample_id"], s["vendor"],
                                           s["acquisition_method"], s["drift_kind"]))
    head.append('<div class=meta>%s</div>'
                % "".join('<span class=chip>%s</span>' % esc(c) for c in chips))

    alerts = []
    if not absolute:
        alerts.append("No absolute standard was declared, so every quant axis is"
                      " RELATIVE abundance. No panel here may be read as a"
                      " concentration.")
    if "quad_window" in kinds:
        alerts.append("A sample is SONAR: its drift column is a scanning-quadrupole"
                      " position, physically unrelated to ion mobility. The drift"
                      " axis is labelled accordingly and must not be compared with"
                      " a mobility axis.")
    if len(set(kinds)) > 1:
        alerts.append("Samples disagree about what drift means (%s). Points from"
                      " different drift kinds share no axis, so the drift panels"
                      " carry %s only and say how many ions they excluded."
                      % (", ".join(sorted(set(kinds))), one_kind))
    if tiers_present == ["bare_proximity"]:
        # The rate is read from the export, never asserted. A caption that
        # quoted a figure from elsewhere would be describing a different run,
        # and for most of this project's life it described retired code.
        pc = manifest.get("proximity_control") or {}
        ent = ((pc.get("arms") or {}).get("entrapment") or {}) if pc else {}
        if ent.get("queried"):
            axes = ", ".join(pc.get("axes_matched") or []) or "unstated axes"
            tail = ("That window's false-positive rate, measured on this run by"
                    " querying peptides that cannot be in the sample, is %.1f%%"
                    " (%d of %d probes, 95%% CI %.1f-%.1f%%), on %s."
                    % (100 * ent["rate"], ent["matched"], ent["queried"],
                       100 * ent["ci"][0], 100 * ent["ci"][1], axes))
        else:
            tail = ("That window's false-positive rate was NOT measured on this"
                    " run, so nothing here says how often it is wrong.")
        alerts.append("Every ion drawn is bare_proximity: an m/z window plus"
                      " consecutive scans, with no decoy control on the"
                      " individual ion. " + tail)
    elif len(tiers_present) > 1:
        alerts.append("Several confidence tiers are drawn together (%s). They are"
                      " marked differently; read the legend before comparing"
                      " points." % ", ".join(tiers_present))
    dec = manifest.get("decimation") or {}
    if dec.get("strategy") and dec["strategy"] != "none" and dec.get("of"):
        alerts.append("Decimated: %s kept %s of %s (%.2f%%). Every panel below shows"
                      " that fraction, not the data."
                      % (dec["strategy"], f"{dec['kept']:,}", f"{dec['of']:,}",
                         100.0 * dec["kept"] / dec["of"]))
    for w in manifest.get("warnings", []):
        # Skip a manifest warning the page has already made in its own words.
        # Two sources for one fact reads as two problems.
        if any(mark in w and any(mark in a for a in alerts)
               for mark in ("false-positive rate", "RELATIVE abundance", "quadrupole")):
            continue
        if w not in alerts:
            alerts.append(esc(w))
    if alerts:
        head.append('<div class=banner><b>Read before the charts</b><ul>%s</ul></div>'
                    % "".join("<li>%s</li>" % a for a in alerts))

    cards, drawn, blanked = [], 0, 0
    for sel in selected:
        pid = sel["id"]
        panel = catalog.get(pid)
        if panel is None:
            continue
        status, reasons = plan_panels.assess(panel, method, manifest)
        title = panel.get("title", pid)
        fam = panel.get("familiarity", "-")
        wide = pid in ("forensic_track", "ion_flight")
        body, notes = None, []
        if status == plan_panels.BLANK:
            why = reasons
        elif pid not in DRAW:
            why = ["no renderer for this panel yet: %s needs a table this tool does"
                   " not draw from" % pid]
            status = plan_panels.BLANK
        else:
            body, notes = DRAW[pid](ions, ctx)
            if body is None:
                why = notes or ["nothing to draw"]
                status = plan_panels.BLANK
                notes = []
            else:
                why = []
        card = ['<div class="card%s">' % (" wide" if wide else "")]
        card.append('<div class=fam>%s &middot; %s</div>' % (esc(fam), esc(status)))
        card.append("<h2>%s</h2>" % esc(title))
        if body is None:
            blanked += 1
            card.append('<div class=blankbox>not drawn%s</div>'
                        % "".join('<span class=why>%s</span>' % esc(w) for w in why))
        else:
            drawn += 1
            card.append(body)
            if reasons:
                notes = list(notes) + ["partial: " + r for r in reasons]
            if notes:
                card.append('<ul class=notes>%s</ul>'
                            % "".join("<li>%s</li>" % esc(n) for n in notes))
        card.append("</div>")
        cards.append("\n".join(card))

    legend = []
    for t in tiers_present:
        colour, _ = TIER_MARK[t]
        legend.append('<span><i class=sw style="background:%s"></i>%s</span>'
                      % (colour, esc(t)))
    legend.append("<span>hover a point for its ion_uid, peptide and tier</span>")

    return """<!doctype html>
<html lang=en><head><meta charset=utf-8>
<meta name=viewport content="width=device-width,initial-scale=1">
<title>%s</title><style>%s</style></head>
<body><div class=wrap>
%s
<div class=grid>
%s
</div>
<div class=legend>%s</div>
<footer>%d panel(s) drawn, %d left blank with the reason on them. Renderability is
decided by tools/plan_panels.py against schema/panels.json, so a panel blanked in
the plan is blanked here. No external scripts, no network: this file is the whole
viewer.</footer>
</div></body></html>
""" % (esc(method.get("name") or manifest["export_id"]), CSS,
       "\n".join(head), "\n".join(cards), "".join(legend), drawn, blanked)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("export", help="export directory or its manifest.json")
    ap.add_argument("--method", help="method JSON naming which panels to draw")
    ap.add_argument("-o", "--out", help="HTML to write (default <export>/panels.html)")
    args = ap.parse_args(argv)

    path = Path(args.export)
    man_path = path if path.is_file() else path / "manifest.json"
    if not man_path.exists():
        print("no manifest at %s" % man_path, file=sys.stderr)
        return 2
    export_dir = man_path.parent
    manifest = json.loads(man_path.read_text())
    ions = read_jsonl(export_dir / "ions.jsonl")
    if not ions:
        # "Nothing to draw" is wrong when a proteins or spectra table is sitting
        # right there. The real reason is the contract: tables.ions is const
        # true, so quant without the ion list it came from is a figure, not an
        # export -- and that is the rule, not a gap in this tool.
        others = [t for t in ("proteins", "spectra", "traces", "fragment_hits")
                  if (export_dir / ("%s.jsonl" % t)).exists()]
        print("no ions.jsonl in %s." % export_dir, file=sys.stderr)
        if others:
            print("The directory does hold %s. Panels could be drawn from %s"
                  " alone, but the contract makes tables.ions const true: quant"
                  " and spectra are derived from the ion list, and without it"
                  " there is no evidence trail under them. This is a set of"
                  " figures, not an export."
                  % (", ".join("%s.jsonl" % t for t in others), others[0]),
                  file=sys.stderr)
        else:
            print("Nothing here to draw, and an empty page would not say so.",
                  file=sys.stderr)
        return 2

    declared = manifest.get("tables") or {}
    tables, missing = {}, []
    for name in ("spectra", "traces", "proteins", "fragment_hits"):
        path = export_dir / ("%s.jsonl" % name)
        rows = read_jsonl(path)
        tables[name] = rows
        if declared.get(name) and not rows:
            missing.append(name)
    want_kinds = set(declared.get("spectrum_kinds") or [])
    have_kinds = {r.get("kind") for r in tables["spectra"]}
    absent_kinds = sorted(want_kinds - have_kinds)
    if absent_kinds and tables["spectra"]:
        manifest.setdefault("warnings", []).append(
            "the manifest declares spectrum kind(s) %s but no row in spectra.jsonl"
            " carries them. The panels that need them are blank below."
            % ", ".join(absent_kinds))
        declared["spectrum_kinds"] = sorted(have_kinds & want_kinds)
    if missing:
        # validate_export is the tool that fails on this; the page says it too,
        # because a reader looking at charts is not going to re-run the validator.
        manifest.setdefault("warnings", []).append(
            "the manifest declares %s but %s holds no such rows. The panels that"
            " need them are blank below, and the export does not describe itself."
            % (", ".join(missing),
               ", ".join("%s.jsonl" % m for m in missing)))
        for name in missing:
            declared[name] = False

    catalog = {p["id"]: p for p in
               json.loads((ROOT / "schema" / "panels.json").read_text())["panels"]}
    if args.method:
        method = json.loads(Path(args.method).read_text())
        selected = method.get("panels", [])
    else:
        # No method: offer every panel this tool can draw, so the absence of a
        # method never silently narrows the page.
        method = {"name": manifest.get("export_id"), "tier_policy": {}, "panels": []}
        selected = [{"id": pid} for pid in catalog if pid in DRAW]
        print("no --method: offering the %d panel(s) this tool can draw"
              % len(selected))

    out = Path(args.out) if args.out else export_dir / "panels.html"
    out.write_text(build(manifest, method, ions, catalog, selected, tables))
    print("wrote %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
