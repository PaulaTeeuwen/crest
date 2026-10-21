#!/usr/bin/env python3
"""run_seeded_search.py -- single entry point for an nci_metal search,
with isomer/ligand seeding as a toggle instead of a separate manual
pipeline.

Without seeding, this is just "run crest once on this structure." With
seeding (the default), it additionally: enumerates every symmetry-distinct
way the input's own ligands could be arranged (generate_isomer_structures
.generate, itself built on isomer_enumerate.distinct_isomers), runs an
independent nci_metal search from the original input AND from every
generated isomer, then merges every resulting ensemble into one
deduplicated set via merge_isomer_ensembles (irmsd-based -compare, so the
same equal-mass-axis-aware alignment used throughout this project decides
what counts as a duplicate).

This replaces the manual three-tool pipeline used throughout this
project's own test suite (generate_isomer_structures.py by hand -> one
crest invocation per isomer directory, run by hand -> merge_isomer_
ensembles.py by hand) with one command and one flag.

Usage:
    # seeded (default): explore every distinct starting arrangement
    python run_seeded_search.py input.xyz --metal 1 --toml-template base.toml \\
        --outdir search_out --crest-bin /path/to/crest

    # unseeded: a single plain run from the given structure only
    python run_seeded_search.py input.xyz --metal 1 --toml-template base.toml \\
        --outdir search_out --crest-bin /path/to/crest --no-seed

The toml template must NOT contain an `input = ...` line (one is written
per run-directory automatically); everything else (runtype, calculation
.level, an optional calculation.nci_metal restraint block, etc.) is
copied through unchanged.
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_isomer_structures as gis
import merge_isomer_ensembles as mie


def read_template(path: str) -> str:
    with open(path) as f:
        text = f.read()
    if re.search(r'^\s*input\s*=', text, re.MULTILINE):
        raise ValueError(
            f"{path}: toml template must not set 'input' itself -- "
            f"run_seeded_search.py writes a per-run-directory input= line"
        )
    return text


def write_run_dir(run_dir: str, template_text: str, xyz_src: str, xyz_name: str,
                   restrain_metal: str | None = None, uff4mof_type: str | None = None) -> None:
    """restrain_metal: when given (an element symbol), appends a
    [[calculation.nci_metal]] restrain=true block for that element. Needed
    for generated isomer structures specifically -- generate_isomer_
    structures.py deliberately EXPANDS donor distances past realistic
    bond lengths to avoid clashes while repositioning a ligand (see its
    own docstring/expand_by), relying on a restraint to pull them back in
    during the run. Without it, GFN-FF's own topology perception can fail
    to recognize the expanded donor as bonded at all -- confirmed
    directly: an unrestrained run on an expanded isomer here fragmented
    into N disconnected pieces and produced a meaningless, unbounded
    "conformer" count from freely tumbling fragments rather than a real
    ensemble. The original (non-generated) input structure has no such
    expansion and does not need this.
    """
    os.makedirs(run_dir, exist_ok=True)
    shutil.copy(xyz_src, os.path.join(run_dir, xyz_name))
    with open(os.path.join(run_dir, "nci_metal.toml"), "w") as f:
        f.write(f'input   = "{xyz_name}"\n')
        f.write(template_text)
        if restrain_metal:
            f.write(f'\n[[calculation.nci_metal]]\natom     = "{restrain_metal}"\nrestrain = true\n')
            if uff4mof_type:
                f.write(f'uff4mof  = "{uff4mof_type}"\nforce    = true\n')


def detect_uff4mof_type(crest_bin: str, outdir: str, template_text: str, input_xyz: str, metal_symbol: str) -> str | None:
    """Runs a `-dry` probe on the ORIGINAL (non-expanded, real) input
    structure with a plain restrain=true block, and parses which UFF4MOF
    type crest's own geometry-based auto-detection suggests for the
    metal. That type is then reused, explicitly (uff4mof=... force=true),
    for every GENERATED isomer -- generate_isomer_structures.py's
    deliberately-expanded donor distances make CN=0 there, so the SAME
    auto-detection run directly on an isomer sees no candidates at all
    (confirmed directly: "CN=0: no unambiguous auto-detection -> none")
    even with restrain=true set; running the probe on the real,
    non-expanded structure instead sidesteps that chicken-and-egg problem.
    Returns None if even the real structure's CN is ambiguous/undetected
    (the caller then falls back to plain restrain=true, same as before).
    """
    # NOTE: -dry mode skips the UFF4MOF coordination-sphere assignment
    # printout entirely (confirmed directly: only present in a real run's
    # output) -- so this has to be a real, if fast, run rather than a
    # dry-run preview.
    probe_dir = os.path.join(outdir, "_uff4mof_probe")
    write_run_dir(probe_dir, template_text + f'\n[[calculation.nci_metal]]\natom     = "{metal_symbol}"\nrestrain = true\n',
                  input_xyz, os.path.basename(input_xyz))
    result = subprocess.run([crest_bin, "nci_metal.toml"], cwd=probe_dir,
                             capture_output=True, text=True)
    m = re.search(r'geometry (?:suggests|best matches) (\S+?)(?:\s|$|->)', result.stdout)
    return m.group(1) if m else None


def run_crest(crest_bin: str, run_dir: str) -> bool:
    """Returns True if crest produced a conformer ensemble."""
    log_path = os.path.join(run_dir, "run.log")
    with open(log_path, "w") as log:
        subprocess.run([crest_bin, "nci_metal.toml"], cwd=run_dir, stdout=log, stderr=subprocess.STDOUT)
    return os.path.exists(os.path.join(run_dir, "crest_conformers.xyz"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_xyz")
    ap.add_argument("--metal", type=int, required=True, help="1-based index of the metal atom")
    ap.add_argument("--toml-template", required=True,
                     help="toml file with everything EXCEPT the input= line "
                          "(runtype, calculation.level, optional restraint block, ...)")
    ap.add_argument("--outdir", default="seeded_search")
    ap.add_argument("--crest-bin", default="crest")
    ap.add_argument("--max-bite-angle", type=float, default=100.0)
    ap.add_argument("--seed", dest="seed", action="store_true", default=True,
                     help="enumerate and run every distinct starting isomer (default)")
    ap.add_argument("--no-seed", dest="seed", action="store_false",
                     help="run a single plain search from input_xyz only, no isomer enumeration")
    ap.add_argument("--rthr", type=float, default=0.125, help="merge RMSD threshold, Angstrom")
    ap.add_argument("--enantiomers", choices=["distinguish", "merge", "auto"], default="distinguish")
    args = ap.parse_args(argv)

    template_text = read_template(args.toml_template)
    os.makedirs(args.outdir, exist_ok=True)

    if not args.seed:
        run_dir = os.path.join(args.outdir, "unseeded")
        xyz_name = os.path.basename(args.input_xyz)
        write_run_dir(run_dir, template_text, args.input_xyz, xyz_name)
        print(f"--no-seed: running a single plain search in {run_dir}")
        ok = run_crest(args.crest_bin, run_dir)
        if not ok:
            print(f"crest did not produce crest_conformers.xyz in {run_dir} -- check {run_dir}/run.log", file=sys.stderr)
            return 1
        final = os.path.join(args.outdir, "final_ensemble.xyz")
        shutil.copy(os.path.join(run_dir, "crest_conformers.xyz"), final)
        print(f"wrote {final}")
        return 0

    # seeded: enumerate isomers, one run-directory per isomer (+ the
    # original input itself), run each, then merge every ensemble found.
    isomer_dir = os.path.join(args.outdir, "isomers")
    print(f"enumerating distinct starting isomers into {isomer_dir} ...")
    gis.generate(args.input_xyz, args.metal - 1, isomer_dir, args.max_bite_angle)

    run_dirs = []

    base_dir = os.path.join(args.outdir, "run_input")
    write_run_dir(base_dir, template_text, args.input_xyz, os.path.basename(args.input_xyz))
    run_dirs.append(base_dir)

    metal_symbol, _ = gis.read_xyz(args.input_xyz)
    metal_symbol = metal_symbol[args.metal - 1]

    uff4mof_type = detect_uff4mof_type(args.crest_bin, args.outdir, template_text, args.input_xyz, metal_symbol)
    if uff4mof_type:
        print(f"detected UFF4MOF type {uff4mof_type} from the real input geometry -- "
              f"forcing it explicitly for every generated isomer (their expanded donor "
              f"distances make auto-detection see CN=0, so it can't be trusted there)")
    else:
        print("WARNING: could not detect a UFF4MOF type from the input structure -- "
              "isomer runs will use plain restrain=true, which may apply no restraint "
              "at all if their own CN also comes out ambiguous/undetected", file=sys.stderr)

    isomer_files = sorted(f for f in os.listdir(isomer_dir) if f.endswith(".xyz"))
    for fname in isomer_files:
        tag = os.path.splitext(fname)[0]
        run_dir = os.path.join(args.outdir, f"run_{tag}")
        write_run_dir(run_dir, template_text, os.path.join(isomer_dir, fname), fname,
                      restrain_metal=metal_symbol, uff4mof_type=uff4mof_type)
        run_dirs.append(run_dir)

    print(f"running {len(run_dirs)} crest searches ({len(isomer_files)} seeded isomer(s) + the original input) ...")
    succeeded = []
    for run_dir in run_dirs:
        print(f"  running {run_dir} ...")
        if run_crest(args.crest_bin, run_dir):
            succeeded.append(os.path.join(run_dir, "crest_conformers.xyz"))
        else:
            print(f"  WARNING: {run_dir} produced no crest_conformers.xyz -- skipped in merge "
                  f"(check {run_dir}/run.log)", file=sys.stderr)

    if not succeeded:
        print("no run produced an ensemble -- nothing to merge", file=sys.stderr)
        return 1

    final = os.path.join(args.outdir, "final_ensemble.xyz")
    if len(succeeded) == 1:
        shutil.copy(succeeded[0], final)
        print(f"only one successful run -- wrote it directly to {final}")
    else:
        print(f"merging {len(succeeded)} ensembles (enantiomers: {args.enantiomers}) ...")
        mie.main([final, *succeeded, "--crest-bin", args.crest_bin,
                  "--rthr", str(args.rthr), "--enantiomers", args.enantiomers])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
