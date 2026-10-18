#!/usr/bin/env python3
"""Run the checks in order, and stop where a gate says the next step is unreadable.

The tools in this repository each answer one question, and the order matters
because an earlier answer can make a later one meaningless. Running them
separately leaves it to a person to remember that — which is the failure this
repository has hit five times: a number computed without asking what else could
have produced it.

The chain, and why each step precedes the next:

  1 EXPORT     validate_export.py — does the export describe itself honestly?
                 A contract violation upstream makes every later number suspect.
  2 PANELS     plan_panels.py — which panels can render from what is there?
                 No point measuring what no panel will draw.
  3 CHEMISTRY  fragment_chemistry.py — do the fragment sources carry chemistry?
                 GATE: a flat source is not a reference. Nothing downstream that
                 compares patterns against it can be read, in either direction.
  4 RESOLVE    resolvability.py — which precursors can be told apart at all?
                 GATE: without a measurement scale there is no answer, and the
                 collision count bounds what any assignment can achieve.
  5 STATS      ladder_stats.py — intervals on the real-vs-decoy rates.
                 GATE: asymmetric ladders make every rate a floor, so the
                 intervals describe a handicapped comparison.
  6 OFFSET     offset_check.py — is a shared offset hiding behind the controls?
                 GATE: a shared offset is invisible to every spread-based check
                 and gets more confidently wrong as sample count rises, so a
                 tightening interval above means nothing until this has run.

    python3 tools/aria_check.py --export export/ --panels method.json manifest.json
    python3 tools/aria_check.py --fragments a.jsonl b.jsonl --ions ions.json
    python3 tools/aria_check.py --all-of examples/

A step with no inputs is skipped and said to be skipped. A gate that fails stops
the chain there: later steps are reported as NOT REACHED rather than silently
omitted, because "we did not run it" and "it passed" must never look alike.

Exit status is 0 only if every step that had inputs ran and passed. 1 if a gate
failed or a step reported findings, 2 if a tool could not be run at all.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
RULE = "=" * 78


class Chain:
    def __init__(self, verbose):
        self.verbose = verbose
        self.results = []          # (name, state, detail)
        self.stopped_by = None

    def step(self, name, argv, gate=False, why_gate=""):
        """Run one tool. A gate that fails stops everything after it."""
        if self.stopped_by:
            self.results.append((name, "NOT REACHED",
                                 "stopped by %s" % self.stopped_by.strip()))
            return None
        print("%s\n%s\n%s" % (RULE, name, RULE))
        proc = subprocess.run([sys.executable] + argv, capture_output=True, text=True)
        body = (proc.stdout or "") + (proc.stderr or "")
        print(body.rstrip() or "(no output)")
        print()
        if proc.returncode == 2:
            self.results.append((name, "COULD NOT RUN", "exit 2"))
            self.stopped_by = name
            return proc
        state = "passed" if proc.returncode == 0 else "findings"
        if proc.returncode != 0 and gate:
            state = "GATE FAILED"
            self.stopped_by = name
            self.results.append((name, state, why_gate))
            return proc
        self.results.append((name, state, ""))
        return proc

    def skip(self, name, reason):
        self.results.append((name, "skipped", reason))

    def report(self):
        print("%s\nSUMMARY\n%s" % (RULE, RULE))
        width = max(len(n) for n, _, _ in self.results)
        for name, state, detail in self.results:
            print("  %-*s  %-14s %s" % (width, name, state, detail))
        bad = [n for n, s, _ in self.results if s in ("GATE FAILED", "findings", "COULD NOT RUN")]
        unreached = [n for n, s, _ in self.results if s == "NOT REACHED"]
        print()
        if self.stopped_by:
            print("  The chain stopped at %s. %d step(s) were NOT REACHED: %s."
                  % (self.stopped_by.strip(), len(unreached),
                     ", ".join(u.strip() for u in unreached) or "none"))
            print("  Not reached is not passed. Fix the gate and run again.")
        elif bad:
            print("  Every step ran. %d reported findings: %s." % (len(bad), ", ".join(bad)))
        else:
            ran = [n for n, s, _ in self.results if s == "passed"]
            print("  %d step(s) passed, %d skipped for want of inputs."
                  % (len(ran), len(self.results) - len(ran)))
            if len(ran) < len(self.results):
                print("  A skipped step has not passed. Supply its inputs to close it.")
        return 2 if any(s == "COULD NOT RUN" for _, s, _ in self.results) else (1 if bad else 0)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--export", help="export directory or manifest for validate_export.py")
    p.add_argument("--excerpt", action="store_true", help="pass --excerpt to validate_export.py")
    p.add_argument("--limit", type=int, help="pass --limit to validate_export.py")
    p.add_argument("--panels", nargs=2, metavar=("METHOD", "MANIFEST"),
                   help="method and manifest for plan_panels.py")
    p.add_argument("--fragments", nargs="+", metavar="JSONL",
                   help="fragment tables for the chemistry gate")
    p.add_argument("--ions", help="precursor coordinates for resolvability.py")
    p.add_argument("--sigma-drift", type=float)
    p.add_argument("--sigma-rt", type=float)
    p.add_argument("--rank-rt", action="store_true")
    p.add_argument("--counts", help="gates/ions JSON for ladder_stats.py")
    p.add_argument("--ppm", help="per-PSM ppm errors for offset_check.py")
    p.add_argument("--anchor", help="an independent reference's ppm, for offset_check.py")
    p.add_argument("--ppm-target", dest="ppm_target", type=float,
                   help="the ppm the Atlas aims for, for offset_check.py")
    p.add_argument("--covariate", action="append", default=[])
    p.add_argument("--all-of", metavar="DIR",
                   help="look in DIR for the conventional filenames and use what is there")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    if args.all_of:
        d = Path(args.all_of)
        args.export = args.export or (str(d / "manifest.example.json")
                                      if (d / "manifest.example.json").exists()
                                      else str(d) if (d / "manifest.json").exists() else None)
        if args.export and "example" in args.export:
            args.excerpt = True
        if not args.panels:
            m = d / "method.example.json"
            mm = d / "manifest.example.json"
            if m.exists() and mm.exists():
                args.panels = [str(m), str(mm)]
        if not args.fragments:
            found = sorted(str(f) for f in d.glob("*fragment*chem*.jsonl"))
            args.fragments = found or None
        if not args.ions and (d / "ions.json").exists():
            args.ions = str(d / "ions.json")
        if not args.counts and (d / "counts.json").exists():
            args.counts = str(d / "counts.json")
        if not args.ppm:
            for name in ("ppm.jsonl", "ppm_errors.jsonl"):
                if (d / name).exists():
                    args.ppm = str(d / name)
                    break
        if not args.anchor and (d / "anchor.jsonl").exists():
            args.anchor = str(d / "anchor.jsonl")

    c = Chain(args.verbose)

    if args.export:
        a = [str(TOOLS / "validate_export.py"), args.export]
        if args.excerpt:
            a.append("--excerpt")
        if args.limit:
            a += ["--limit", str(args.limit)]
        c.step("1 EXPORT    contract", a)
    else:
        c.skip("1 EXPORT    contract", "no --export given")

    if args.panels:
        c.step("2 PANELS    renderability",
               [str(TOOLS / "plan_panels.py")] + list(args.panels))
    else:
        c.skip("2 PANELS    renderability", "no --panels given")

    if args.fragments:
        c.step("3 CHEMISTRY gate", [str(TOOLS / "fragment_chemistry.py")] + list(args.fragments),
               gate=True,
               why_gate="a source carries no fragmentation chemistry, so a pattern "
                        "comparison against it cannot be read in either direction")
    else:
        c.skip("3 CHEMISTRY gate", "no --fragments given")

    if args.ions:
        a = [str(TOOLS / "resolvability.py"), args.ions]
        if args.sigma_drift is not None:
            a += ["--sigma-drift", str(args.sigma_drift)]
        if args.sigma_rt is not None:
            a += ["--sigma-rt", str(args.sigma_rt)]
        if args.sigma_drift is None or args.sigma_rt is None:
            a.append("--estimate-sigma")
        if args.rank_rt:
            a.append("--rank-rt")
        # Not a gate: collisions are a fact to carry forward, not a reason to
        # stop. Only a MISSING measurement scale is fatal, and the tool itself
        # exits 2 for that, which the chain treats as COULD NOT RUN.
        c.step("4 RESOLVE   collisions", a)
    else:
        c.skip("4 RESOLVE   collisions", "no --ions given")

    if args.counts:
        a = [str(TOOLS / "ladder_stats.py"), args.counts]
        for name in args.covariate:
            a += ["--covariate", name]
        c.step("5 STATS     intervals", a, gate=True,
               why_gate="the ladders are asymmetric, so every rate above is a floor "
                        "and the intervals describe a handicapped comparison")
    else:
        c.skip("5 STATS     intervals", "no --counts given")

    if args.ppm:
        a = [str(TOOLS / "offset_check.py"), args.ppm]
        if args.anchor:
            a += ["--anchor", args.anchor]
        if args.ppm_target is not None:
            a += ["--target", str(args.ppm_target)]
        # A gate, because the controls above cannot see what this looks for: a
        # tightening interval on the mass dimension means nothing until a shared
        # offset has been ruled out, and it gets MORE convincing as one persists.
        c.step("6 OFFSET    shared shift", a, gate=True,
               why_gate="a shared offset is present, so every spread-based control "
                        "above reports it as health and more samples make the "
                        "coordinate more confidently wrong")
    else:
        c.skip("6 OFFSET    shared shift", "no --ppm given")

    return c.report()


if __name__ == "__main__":
    sys.exit(main())
