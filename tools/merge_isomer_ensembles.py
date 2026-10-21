#!/usr/bin/env python3
"""merge_isomer_ensembles.py -- combine the separate crest_conformers.xyz
outputs from multiple isomer-seeded nci_metal runs (one per distinct
isomer from generate_isomer_structures.py) into ONE deduplicated
conformer database, using crest's own `-compare` tool as the comparison
engine rather than a hand-rolled RMSD.

Why not just concatenate the files: different starting isomers can
relax into the SAME actual conformer (we saw this happen for Berry
pseudorotation and Bailar-twist test cases this session), so a naive
concatenation double-counts real duplicates. But a NAIVE deduplication
is just as wrong in the other direction: comparing two mirror-image
structures (e.g. a Bailar-twist Delta/Lambda pair) with the wrong
settings can report RMSD~0 and silently merge two genuinely distinct
enantiomers into one -- confirmed directly earlier this session (crest's
own `-irmsd` on a Delta/Lambda Co(en)3 pair: ~8e-8 A under auto-detected
stereocheck, because has_stereo only ever recognizes 4-coordinate
centers, never an octahedral metal's). This script always passes
`-inversion off` to `-compare` for exactly that reason: never allow a
false-enantiomer match, so chirality the search actually found never
gets discarded as a "duplicate" of its own mirror image.

`-compare` itself is a pure analysis tool (writes rmsdmatch.dat, doesn't
emit a merged ensemble) -- this script drives it pairwise and builds the
merged ensemble from its match report: keep everything from the running
merged set, append only the structures from the next ensemble that have
NO match (within --rthr) against anything already kept.

Requirement (inherited from crest's own compare_ensembles): all input
ensembles must have IDENTICAL atom count AND atom ORDER -- satisfied
automatically for any set of isomers built via generate_isomer_structures
.py, since it only ever moves existing atoms by their fixed index, never
reorders them.

Usage:
    python merge_isomer_ensembles.py out_merged.xyz \\
        run1/crest_conformers.xyz run2/crest_conformers.xyz ...
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

ENERGY_RE = re.compile(r"energy=([-+0-9.eEdD]+)")


def split_frames(path: str) -> list[tuple[int, str, list[str]]]:
    """Returns [(natoms, comment_line, atom_lines), ...] for an xyz ensemble."""
    with open(path) as f:
        lines = f.readlines()
    frames = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].strip())
        comment = lines[i + 1]
        atoms = lines[i + 2:i + 2 + n]
        frames.append((n, comment, atoms))
        i += 2 + n
    return frames


def frame_energy(comment: str) -> float:
    m = ENERGY_RE.search(comment)
    if m:
        return float(m.group(1).replace("D", "E").replace("d", "e"))
    # fall back: bare float somewhere on the comment line (plain crest convention)
    m = re.search(r"[-+]?\d+\.\d+(?:[eEdD][-+]?\d+)?", comment)
    if m:
        return float(m.group(0).replace("D", "E").replace("d", "e"))
    raise ValueError(f"could not find an energy value in comment line: {comment!r}")


def write_frames(path: str, frames: list[tuple[int, str, list[str]]]) -> None:
    with open(path, "w") as f:
        for n, comment, atoms in frames:
            f.write(f"{n}\n{comment}" if comment.endswith("\n") else f"{n}\n{comment}\n")
            f.writelines(atoms)


def run_compare(crest_bin: str, ens_a: str, ens_b: str, rthr: float, maxcompare: int,
                 workdir: str, enantiomers: str) -> list[tuple[int, int, float]]:
    """Runs `crest -compare ens_a ens_b`, returns the rmsdmatch.dat
    content as [(index_in_a, index_in_b, rmsd), ...] (1-based indices,
    matching crest's own convention).

    `enantiomers` controls whether a mirror-image match is allowed to
    count as a duplicate:
      "distinguish" (default, recommended) -- passes `-inversion off`:
        NEVER treat a structure as a duplicate of its own mirror image,
        so a genuine enantiomer pair (e.g. a Bailar-twist Delta/Lambda)
        the search actually found survives the merge as two entries.
      "merge" -- passes `-inversion on`: ALWAYS allow a mirror-image
        match, i.e. deliberately treat enantiomers as the same entry.
        Use this only if chirality genuinely doesn't matter for your
        purpose and you want the smaller, chirality-collapsed database.
      "auto" -- no override; defers to crest's own has_stereo detection.
        NOT recommended for metal complexes: has_stereo only ever
        recognizes 4-coordinate (tetrahedral) stereocenters, never a 5-
        or 6-coordinate metal one, so "auto" silently behaves like
        "merge" for exactly the chiral octahedral/TBP cases this tool
        exists for (confirmed directly: Delta/Lambda Co(en)3 gives
        RMSD~8e-8 A under "auto", vs the correct ~0.93 A under
        "distinguish").

    A real structure must be passed to crest as a bare positional
    argument alongside `-compare` -- crest's generic startup always
    looks for one, even though compare_ensembles() itself never uses
    it (it reads ens_a/ens_b directly). The first frame of ens_a is
    reused for this; its content is otherwise irrelevant here.
    """
    dummy_input = os.path.join(workdir, "_dummy_main_input.xyz")
    write_frames(dummy_input, split_frames(ens_a)[:1])
    inversion_arg = {"distinguish": "off", "merge": "on", "auto": None}[enantiomers]
    # the CLI flag is "-maxcomp", not "-maxcompare" (an earlier version
    # of this script used the wrong name, which crest silently rejected
    # as an unrecognized argument -- not a crest bug, just this script's
    # own mistake. -maxcomp confirmed working directly: "N lowest" in
    # crest's own printout matches the full ensemble size, not the
    # default-100 cap).
    cmd = [crest_bin, dummy_input, "-compare", os.path.abspath(ens_a), os.path.abspath(ens_b),
           "-rthr", str(rthr), "-maxcomp", str(maxcompare)]
    if inversion_arg is not None:
        cmd += ["-inversion", inversion_arg]
    result = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    matchfile = os.path.join(workdir, "rmsdmatch.dat")
    if not os.path.exists(matchfile):
        raise RuntimeError(f"crest -compare did not produce rmsdmatch.dat; stderr:\n{result.stderr[-2000:]}")
    matches = []
    with open(matchfile) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                matches.append((int(parts[0]), int(parts[1]), float(parts[2])))
    return matches


def merge_two(crest_bin: str, merged_frames: list, next_path: str, rthr: float,
               maxcompare: int, scratch: str, enantiomers: str) -> list:
    """merged_frames (already-accumulated, energy-sorted) + next_path's
    unmatched structures -> new merged_frames list."""
    merged_path = os.path.join(scratch, "merged_current.xyz")
    write_frames(merged_path, merged_frames)

    next_frames = split_frames(next_path)
    matches = run_compare(crest_bin, merged_path, next_path, rthr, maxcompare, scratch, enantiomers)
    # matches are (index in merged, index in next) over the ncomp1 x ncomp2
    # lowest-maxcompare window -- anything beyond that window was never
    # checked and is conservatively treated as unmatched (kept).
    matched_in_next = {j for (_, j, _) in matches}
    kept_from_next = [fr for k, fr in enumerate(next_frames, start=1) if k not in matched_in_next]
    print(f"  {next_path}: {len(next_frames)} structures, {len(matched_in_next)} matched "
          f"an existing one (discarded), {len(kept_from_next)} new")
    combined = merged_frames + kept_from_next
    combined.sort(key=lambda fr: frame_energy(fr[1]))
    return combined


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output", help="path for the merged, deduplicated ensemble xyz")
    ap.add_argument("ensembles", nargs="+", help="crest_conformers.xyz files to merge, one per isomer-seeded run")
    ap.add_argument("--crest-bin", default="crest")
    ap.add_argument("--rthr", type=float, default=0.125, help="RMSD threshold in Angstrom (crest's own default)")
    ap.add_argument("--maxcompare", type=int, default=5000,
                     help="max structures compared per ensemble (crest's own -maxcomp default is only 100)")
    ap.add_argument("--enantiomers", choices=["distinguish", "merge", "auto"], default="distinguish",
                     help='whether to treat mirror-image structures as duplicates: "distinguish" (default) '
                          'never merges a structure with its own enantiomer, so real chirality the search '
                          'found (e.g. a Bailar-twist Delta/Lambda pair) survives as two entries; "merge" '
                          'deliberately collapses enantiomer pairs into one if chirality doesn\'t matter for '
                          'your purpose; "auto" defers to crest\'s own stereocenter detection, which is NOT '
                          'recommended here since it never recognizes a 5-/6-coordinate metal center as '
                          'chiral and so behaves like "merge" for exactly the cases this tool is for.')
    args = ap.parse_args(argv)

    if len(args.ensembles) < 2:
        print("need at least 2 ensembles to merge", file=sys.stderr)
        return 1

    merged = split_frames(args.ensembles[0])
    merged.sort(key=lambda fr: frame_energy(fr[1]))
    print(f"starting from {args.ensembles[0]}: {len(merged)} structures")
    print(f"enantiomer handling: {args.enantiomers}")

    with tempfile.TemporaryDirectory(prefix="merge_isomer_") as scratch:
        for path in args.ensembles[1:]:
            merged = merge_two(args.crest_bin, merged, path, args.rthr, args.maxcompare, scratch, args.enantiomers)

    write_frames(args.output, merged)
    print(f"\nwrote {len(merged)} deduplicated structures to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
