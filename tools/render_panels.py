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


def label_of(ion):
    ann = ion.get("annotation") or {}
    pep = ann.get("peptide_sequence") or "unassigned"
    return "%s  %s %d+  rt %.2f  %s" % (
        ion["ion_uid"], pep, ion["coord"]["charge"], ion["coord"]["rt_min"],
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


DRAW = {"intensity_vs_time": panel_intensity_vs_time,
        "drift_vs_time": panel_drift_vs_time,
        "mz_vs_drift": panel_mz_vs_drift,
        "seqstart_vs_intensity": panel_seqstart_vs_intensity,
        "forensic_track": panel_forensic_track,
        "ion_flight": panel_ion_flight}

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


def build(manifest, method, ions, catalog, selected):
    samples = manifest.get("samples", [])
    kinds = [s["drift_kind"] for s in samples]
    quant = manifest.get("quant") or {}
    absolute = quant.get("absolute_standard_declared") or \
        (method.get("quant", {}) or {}).get("absolute_standard")
    one_kind = axis_kind(kinds)
    ctx = {"axis_kind": one_kind,
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
        alerts.append("Every ion drawn is bare_proximity: an m/z window plus"
                      " consecutive scans, with no decoy control. That bare"
                      " mechanism passed reversed-tryptic decoys at 33-70%"
                      " elsewhere in this project.")
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
        print("no ions.jsonl in %s: there is nothing to draw, and an empty page"
              " would not say so" % export_dir, file=sys.stderr)
        return 2

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
    out.write_text(build(manifest, method, ions, catalog, selected))
    print("wrote %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
