#!/usr/bin/env python3
"""Name the vendor and acquisition method from what is actually in a raw folder.

The pipeline cannot pick tolerances until it knows these two, and getting the
second one wrong is worse than not knowing it: `drift` is ion mobility on
HDMSE and diaPASEF, a scanning-quadrupole position on Waters SONAR, and absent
on DDA. A tolerance keyed on vendor alone treats a quadrupole position as a
mobility measurement.

So this reports evidence, not a guess. Every answer carries the files it was
read from, and anything it cannot establish comes back as `unknown` rather than
as a default.

    python3 tools/detect_vendor.py /path/to/raw_folder
    python3 tools/detect_vendor.py /path/to/raw_folder --json

What each vendor looks like on disk:

    Waters    a .raw DIRECTORY holding _FUNC00n.DAT and _extern.inf
    Bruker    a .d directory holding analysis.tdf (+ .tdf_bin) or analysis.baf
    Thermo    a .raw FILE, not a directory
    Sciex     .wiff with a .wiff.scan beside it
    Agilent   a .d directory holding AcqData/
    mzML      open format; the vendor is whatever wrote it, named inside

Acquisition method, where the folder says so:

    Waters    _extern.inf text: a SONAR block, an HDMSE/MSE function table, or
              neither, which means DDA
    Bruker    the tdf's own metadata names diaPASEF, PASEF or DDA
    mzML      the instrument configuration and scan descriptions

Exit status is 0 when both vendor and acquisition method are established, 1 when
either is unknown, 2 when the path is not a readable directory.
"""

import argparse
import json
import os
import re
import sys

# drift_kind follows from the acquisition method, never from the vendor alone.
DRIFT_KIND = {
    "HDMSE": "ion_mobility", "UDMSE": "ion_mobility", "MSE": "absent",
    "SONAR": "quad_window", "diaPASEF": "ion_mobility", "PASEF": "ion_mobility",
    "DDA": "absent",
}


def listdir(path):
    try:
        return sorted(os.listdir(path))
    except OSError:
        return []


def read_text(path, limit=200_000):
    """Vendor .inf files are latin-1 or UTF-16; read defensively."""
    try:
        with open(path, "rb") as fh:
            raw = fh.read(limit)
    except OSError:
        return ""
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        encodings = ("utf-16", "utf-8", "latin-1")
    elif raw.count(b"\x00") > len(raw) // 4:
        # BOM-less UTF-16 decodes cleanly as UTF-8 and comes back
        # NUL-interleaved, so every keyword test silently misses. The NUL
        # density, not the decoder's verdict, is what gives it away.
        encodings = ("utf-16-le", "utf-16-be", "utf-8", "latin-1")
    else:
        encodings = ("utf-8", "utf-16", "latin-1")
    for enc in encodings:
        try:
            text = raw.decode(enc)
        except (UnicodeDecodeError, UnicodeError):
            continue
        if text.count("\x00") > len(text) // 8:
            continue            # decoded, but not into readable text
        return text
    return ""


def waters_method(folder, evidence):
    """Read _extern.inf, then fall back to counting functions."""
    inf = next((f for f in listdir(folder) if f.lower() == "_extern.inf"), None)
    if inf:
        text = read_text(os.path.join(folder, inf))
        if text:
            evidence.append("read %s (%d chars)" % (inf, len(text)))
            low = text.lower()
            if "sonar" in low:
                return "SONAR", "the _extern.inf names SONAR"
            if re.search(r"\budmse\b", low):
                return "UDMSE", "the _extern.inf names UDMSE"
            if re.search(r"\bhdmse\b|\bhd-mse\b", low):
                return "HDMSE", "the _extern.inf names HDMSE"
            if re.search(r"\bmse\b", low):
                return "MSE", "the _extern.inf names MSE, and names no mobility mode"
            if re.search(r"\bdda\b|data[- ]dependent", low):
                return "DDA", "the _extern.inf names DDA"
            # A bare mention of mobility is not a mode. "Mobility Trap: off"
            # mentions it too, and reading that as HDMSE puts a drift tolerance
            # on a run with no drift axis.
            if re.search(r"mobility|\bims\b", low):
                on = re.search(r"mobility[^\n]{0,40}?\b(on|enabled|true|yes)\b", low)
                off = re.search(r"mobility[^\n]{0,40}?\b(off|disabled|false|no|none)\b", low)
                if on and not off:
                    return "HDMSE", "the _extern.inf says mobility is on"
                if off and not on:
                    return "MSE", "the _extern.inf says mobility is off"
                return "unknown", ("the _extern.inf mentions mobility without saying"
                                   " whether it was on, and names no acquisition mode")
    funcs = [f for f in listdir(folder) if re.match(r"_FUNC\d+\.DAT$", f, re.I)]
    if funcs:
        evidence.append("found %d _FUNCnn.DAT files" % len(funcs))
        if len(funcs) >= 2:
            # Two functions is the MSE/HDMSE shape (low and high energy). It does
            # not say whether mobility was on, so this is not enough.
            return "unknown", ("%d acquisition functions, which is the MSE family,"
                               " but nothing here says whether mobility was on"
                               % len(funcs))
        return "DDA", "a single acquisition function"
    return "unknown", "no _extern.inf and no _FUNCnn.DAT files to read"


def bruker_method(folder, evidence):
    names = {f.lower(): f for f in listdir(folder)}
    if "analysis.tdf" in names:
        tdf = names["analysis.tdf"]      # open the real name, not the folded one
        evidence.append("found %s" % tdf)
        # The tdf is SQLite; its GlobalMetadata table names the scan mode.
        try:
            import sqlite3
            con = sqlite3.connect("file:%s?mode=ro"
                                  % os.path.join(folder, tdf), uri=True)
            rows = dict(con.execute(
                "SELECT Key, Value FROM GlobalMetadata").fetchall())
            con.close()
            blob = " ".join(str(v) for v in rows.values()).lower()
            evidence.append("read GlobalMetadata (%d keys)" % len(rows))
            if "diapasef" in blob or "dia-pasef" in blob:
                return "diaPASEF", "GlobalMetadata names diaPASEF"
            if "pasef" in blob:
                return "PASEF", "GlobalMetadata names PASEF"
            if "dda" in blob:
                return "DDA", "GlobalMetadata names DDA"
            return "unknown", "analysis.tdf read but no scan mode named in GlobalMetadata"
        except Exception as exc:                      # a tdf we cannot open
            return "unknown", "analysis.tdf present but unreadable (%s)" % type(exc).__name__
    if "analysis.baf" in names:
        evidence.append("found %s" % names["analysis.baf"])
        return "unknown", "analysis.baf is the older format and names no scan mode here"
    return "unknown", "no analysis.tdf or analysis.baf"


def detect(path):
    """Return a record naming vendor, acquisition method, drift kind and evidence."""
    out = {"path": path, "vendor": "unknown", "acquisition_method": "unknown",
           "drift_kind": "unknown", "evidence": [], "why_unknown": []}
    if not os.path.isdir(path):
        out["why_unknown"].append("%s is not a readable directory" % path)
        return out
    names = listdir(path)
    lower = {n.lower(): n for n in names}
    ev = out["evidence"]

    # A Waters .raw is a directory; a Thermo .raw is a file. The extension alone
    # does not distinguish them, which is why this looks at what the entry IS.
    if path.rstrip("/\\").lower().endswith(".raw"):
        out["vendor"] = "Waters"
        ev.append("the folder itself ends in .raw and is a directory")
        out["acquisition_method"], why = waters_method(path, ev)
    elif path.rstrip("/\\").lower().endswith(".d"):
        if "analysis.tdf" in lower or "analysis.baf" in lower:
            out["vendor"] = "Bruker"
            ev.append("a .d directory holding %s"
                      % ("analysis.tdf" if "analysis.tdf" in lower else "analysis.baf"))
            out["acquisition_method"], why = bruker_method(path, ev)
        elif "acqdata" in lower:
            out["vendor"] = "Agilent"
            ev.append("a .d directory holding AcqData/")
            out["acquisition_method"], why = "unknown", "Agilent method reading is not implemented"
        else:
            out["why_unknown"].append("a .d directory with neither analysis.tdf,"
                                      " analysis.baf nor AcqData/")
            return out
    else:
        # A folder of runs rather than one run: report what is inside.
        kinds = {}
        for n in names:
            full = os.path.join(path, n)
            low = n.lower()
            if low.endswith(".raw"):
                kinds.setdefault("Waters" if os.path.isdir(full) else "Thermo", []).append(n)
            elif low.endswith(".d") and os.path.isdir(full):
                kinds.setdefault("Bruker/Agilent", []).append(n)
            elif low.endswith(".wiff"):
                kinds.setdefault("Sciex", []).append(n)
            elif low.endswith((".mzml", ".mzml.gz")):
                kinds.setdefault("mzML", []).append(n)
        if not kinds:
            out["why_unknown"].append("no vendor run found directly inside %s" % path)
            return out
        if len(kinds) > 1:
            out["why_unknown"].append(
                "more than one vendor inside one folder: %s. Each needs its own"
                " tolerances, so they cannot share a run."
                % ", ".join("%s (%d)" % (k, len(v)) for k, v in sorted(kinds.items())))
            out["contents"] = {k: v for k, v in kinds.items()}
            return out
        vendor, runs = next(iter(kinds.items()))
        out["vendor"] = vendor
        out["runs"] = runs
        ev.append("%d %s run(s) inside the folder" % (len(runs), vendor))
        out["why_unknown"].append(
            "this is a folder of runs, not one run. Detect each run separately:"
            " the acquisition method is per run and must not be assumed shared.")
        return out

    out["drift_kind"] = DRIFT_KIND.get(out["acquisition_method"], "unknown")
    if out["acquisition_method"] == "unknown":
        out["why_unknown"].append(why)
    else:
        ev.append(why)
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("path")
    p.add_argument("--json", action="store_true", help="machine-readable record")
    args = p.parse_args(argv)

    rec = detect(args.path)
    if args.json:
        print(json.dumps(rec, indent=1))
    else:
        print("path                %s" % rec["path"])
        print("vendor              %s" % rec["vendor"])
        print("acquisition_method  %s" % rec["acquisition_method"])
        print("drift_kind          %s" % rec["drift_kind"])
        if rec["drift_kind"] == "quad_window":
            print("                    NOT mobility: a scanning-quadrupole position.")
            print("                    A mobility tolerance must not be applied to it.")
        elif rec["drift_kind"] == "absent":
            print("                    no drift dimension; drift panels stay blank")
        for e in rec["evidence"]:
            print("  evidence          %s" % e)
        for w in rec["why_unknown"]:
            print("  UNKNOWN           %s" % w)
        for k in ("runs", "contents"):
            if rec.get(k):
                print("  %-17s %s" % (k, rec[k]))
    if not os.path.isdir(args.path):
        return 2
    return 0 if "unknown" not in (rec["vendor"], rec["acquisition_method"]) else 1


if __name__ == "__main__":
    sys.exit(main())
