#!/usr/bin/env python3
"""Prove the detector reads evidence and refuses to guess.

The case that matters most: Waters SONAR puts a scanning-quadrupole position in
the drift column, so a run detected as "Waters" without its acquisition method
would get a mobility tolerance applied to a quadrupole position. The detector
must therefore refuse "Waters, method unknown" rather than default.

    python3 tools/test_detect_vendor.py
"""

import json
import os
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "tools" / "detect_vendor.py"
sys.path.insert(0, str(ROOT / "tools"))
import detect_vendor as dv  # noqa: E402


def run(path, *extra):
    p = subprocess.run([sys.executable, str(TOOL), str(path)] + list(extra),
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def waters(tmp, name, inf_text=None, n_funcs=2):
    d = tmp / name
    d.mkdir(parents=True, exist_ok=True)
    for i in range(1, n_funcs + 1):
        (d / ("_FUNC%03d.DAT" % i)).write_bytes(b"\x00" * 16)
    if inf_text is not None:
        (d / "_extern.inf").write_text(inf_text)
    return d


def bruker(tmp, name, mode=None):
    d = tmp / name
    d.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(d / "analysis.tdf")
    con.execute("CREATE TABLE GlobalMetadata (Key TEXT, Value TEXT)")
    con.execute("INSERT INTO GlobalMetadata VALUES ('InstrumentName','timsTOF Pro')")
    if mode:
        con.execute("INSERT INTO GlobalMetadata VALUES ('AcquisitionMode', ?)", (mode,))
    con.commit()
    con.close()
    (d / "analysis.tdf_bin").write_bytes(b"\x00" * 8)
    return d


def sonar_is_not_mobility(tmp):
    """The single most consequential call this tool makes."""
    d = waters(tmp, "sonar_run.raw", "Acquisition: SONAR\nQuadrupole scanning enabled\n")
    rc, out = run(d)
    bad = []
    if rc != 0:
        bad.append("exit %d on a run it could fully identify" % rc)
    if "vendor              Waters" not in out:
        bad.append("did not identify Waters")
    if "acquisition_method  SONAR" not in out:
        bad.append("did not read SONAR from the _extern.inf")
    if "drift_kind          quad_window" not in out:
        bad.append("gave SONAR a drift kind other than quad_window")
    if "NOT mobility" not in out:
        bad.append("did not warn that the drift column is not mobility")
    return bad


def hdmse_is_mobility(tmp):
    d = waters(tmp, "hdmse_run.raw", "Acquisition: HDMSE\nMobility: on\n")
    rc, out = run(d)
    bad = []
    if "acquisition_method  HDMSE" not in out:
        bad.append("did not read HDMSE")
    if "drift_kind          ion_mobility" not in out:
        bad.append("HDMSE did not map to ion_mobility")
    if "NOT mobility" in out:
        bad.append("warned about quad_window on a mobility run")
    return bad


def two_functions_alone_is_not_enough(tmp):
    """Two functions is the MSE family but says nothing about mobility.

    Defaulting here is how a quadrupole position acquires a mobility tolerance.
    """
    d = waters(tmp, "nomethod.raw", inf_text=None, n_funcs=2)
    rc, out = run(d)
    bad = []
    if rc != 1:
        bad.append("exit %d: guessed a method from the function count" % rc)
    if "acquisition_method  unknown" not in out:
        bad.append("did not report the method as unknown")
    if "whether mobility was on" not in out:
        bad.append("did not say what the evidence failed to establish")
    if "vendor              Waters" not in out:
        bad.append("lost the vendor, which the folder shape does establish")
    return bad


def bruker_reads_its_own_metadata(tmp):
    d = bruker(tmp, "pasef_run.d", "diaPASEF")
    rc, out = run(d)
    bad = []
    if rc != 0:
        bad.append("exit %d on a readable tdf" % rc)
    if "acquisition_method  diaPASEF" not in out:
        bad.append("did not read diaPASEF from GlobalMetadata")
    if "drift_kind          ion_mobility" not in out:
        bad.append("diaPASEF did not map to ion_mobility")
    return bad


def a_tdf_naming_no_mode_is_unknown(tmp):
    d = bruker(tmp, "silent_run.d", None)
    rc, out = run(d)
    bad = []
    if rc != 1:
        bad.append("exit %d: invented a mode the tdf does not name" % rc)
    if "no scan mode named" not in out:
        bad.append("did not say the metadata named no scan mode")
    if "vendor              Bruker" not in out:
        bad.append("lost the vendor")
    return bad


def thermo_raw_is_a_file_not_a_directory(tmp):
    """.raw means Waters when it is a directory and Thermo when it is a file."""
    d = tmp / "folder_of_runs"
    d.mkdir(exist_ok=True)
    (d / "sample_one.raw").write_bytes(b"\x01" * 32)
    (d / "sample_two.raw").write_bytes(b"\x01" * 32)
    rc, out = run(d)
    bad = []
    if "Thermo" not in out:
        bad.append("a .raw FILE was not read as Thermo")
    if "Waters" in out:
        bad.append("called .raw files Waters, which is the directory form")
    if "folder of runs" not in out:
        bad.append("did not say this is a folder of runs rather than one run")
    return bad


def mixed_vendors_in_one_folder_are_refused(tmp):
    d = tmp / "mixed"
    d.mkdir(exist_ok=True)
    (d / "a.raw").write_bytes(b"\x01" * 8)
    (d / "b.d").mkdir(exist_ok=True)
    (d / "b.d" / "analysis.tdf").write_bytes(b"\x00" * 8)
    rc, out = run(d)
    bad = []
    if rc != 1:
        bad.append("exit %d: accepted two vendors in one folder" % rc)
    if "more than one vendor" not in out:
        bad.append("did not refuse a mixed folder")
    if "own tolerances" not in out:
        bad.append("did not say why they cannot share a run")
    return bad


def a_folder_of_runs_does_not_claim_one_method(tmp):
    """The method is per run; assuming a folder shares one is how SONAR slips in."""
    d = tmp / "many.raw_folder"
    d.mkdir(exist_ok=True)
    for n in ("r1.raw", "r2.raw"):
        sub = d / n
        sub.mkdir(exist_ok=True)
        (sub / "_FUNC001.DAT").write_bytes(b"\x00" * 8)
    rc, out = run(d)
    bad = []
    if rc != 1:
        bad.append("exit %d: claimed a method for a folder of runs" % rc)
    if "must not be assumed shared" not in out:
        bad.append("did not say the method is per run")
    return bad


def a_missing_path_exits_two(tmp):
    rc, out = run(tmp / "nope")
    return [] if rc == 2 else ["exit %d for a path that is not a directory" % rc]


def json_output_is_machine_readable(tmp):
    d = waters(tmp, "json_run.raw", "SONAR\n")
    rc, out = run(d, "--json")
    try:
        rec = json.loads(out)
    except ValueError:
        return ["--json did not emit parseable JSON"]
    bad = []
    for k in ("vendor", "acquisition_method", "drift_kind", "evidence"):
        if k not in rec:
            bad.append("--json record has no %r" % k)
    if rec.get("drift_kind") != "quad_window":
        bad.append("--json drift_kind is %r" % rec.get("drift_kind"))
    return bad


CASES = [sonar_is_not_mobility, hdmse_is_mobility, two_functions_alone_is_not_enough,
         bruker_reads_its_own_metadata, a_tdf_naming_no_mode_is_unknown,
         thermo_raw_is_a_file_not_a_directory, mixed_vendors_in_one_folder_are_refused,
         a_folder_of_runs_does_not_claim_one_method, a_missing_path_exits_two,
         json_output_is_machine_readable]


def main():
    failed = 0
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for case in CASES:
            bad = case(tmp)
            failed += bool(bad)
            print("%s %s" % ("ok  " if not bad else "FAIL", case.__name__))
            for line in bad:
                print("     %s" % line)
    print("\n%d/%d passed" % (len(CASES) - failed, len(CASES)))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
