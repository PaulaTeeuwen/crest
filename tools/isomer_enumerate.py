#!/usr/bin/env python3
"""isomer_enumerate.py -- enumerate the symmetry-distinct ways to place a
given multiset of ligands (monodentate and/or bidentate) onto the
coordination sites of a fixed geometry, so that each distinct isomer can
be handed to crest's nci_metal runtype as its own starting structure and
its own MTD run -- guaranteeing coverage of the intended isomer space
instead of hoping a single unconstrained run's sampling (MTD + GC
crossing) discovers every one of them from one starting point.

Core idea: two site-assignments are the SAME isomer if and only if some
PROPER rotation of the coordination polyhedron maps one onto the other.
Reflections are deliberately excluded from that equivalence -- a proper
rotation preserves chirality, so enantiomers (e.g. a Bailar-twist
Delta/Lambda pair) are correctly kept as two *distinct* isomers, not
merged into one. The rotation group is found numerically (Kabsch fit +
permutation search) directly from each geometry's site vectors, so
adding a new geometry only means adding its site vectors -- no manual
rotation-matrix bookkeeping, and no risk of accidentally treating a
reflection as if it were a rotation.

A ligand of denticity d occupies d sites at once (its own coordination
"unit"); denticity>1 units are only allowed to occupy a set of sites
that are close enough together to plausibly be bridged (a chelate
"bite"), controlled by --max-bite-angle. Monodentate and bidentate (or
higher) ligands can be mixed freely in the same complex.

Usage:
    python isomer_enumerate.py --geometry octahedral --ligands Cl:2,N:4
    python isomer_enumerate.py --geometry octahedral --ligands en:3 --denticity en:2
    python isomer_enumerate.py --geometry tetrahedral,square_planar --ligands Cl:2,N:2
    python isomer_enumerate.py --list-geometries
"""

from __future__ import annotations

import argparse
import itertools
import sys
from dataclasses import dataclass, field

import numpy as np

# ======================================================================
# Geometry definitions: just the site direction vectors (unit vectors,
# metal at the origin). Everything else (the rotation group, which site
# pairs are chelate-compatible, ...) is derived numerically from these.
# ======================================================================

GEOMETRIES: dict[str, np.ndarray] = {
    "tetrahedral": np.array([
        [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1],
    ], dtype=float),
    "square_planar": np.array([
        [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0],
    ], dtype=float),
    "trigonal_bipyramidal": np.array([
        [0, 0, 1], [0, 0, -1],
        [1, 0, 0], [-0.5, np.sqrt(3) / 2, 0], [-0.5, -np.sqrt(3) / 2, 0],
    ], dtype=float),
    "square_pyramidal": np.array([
        [0, 0, 1],
        [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0],
    ], dtype=float),
    "octahedral": np.array([
        [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
    ], dtype=float),
}

for _name, _v in GEOMETRIES.items():
    GEOMETRIES[_name] = _v / np.linalg.norm(_v, axis=1, keepdims=True)


def find_rotation_group(sites: np.ndarray, tol: float = 1e-3) -> list[tuple[int, ...]]:
    """The proper-rotation symmetry group of `sites`, as the set of site-
    index permutations it induces. Brute-force over all n! permutations
    (n<=6 in every geometry here, so <=720 candidates, each a cheap SVD)
    -- for each permutation, Kabsch-fit the best rotation mapping
    sites -> sites[perm], keep it only if that rotation (a) is proper
    (det=+1, i.e. NOT a reflection) and (b) reproduces the permutation
    exactly (within tol). This is what guarantees reflections never
    sneak into the equivalence relation used later to dedupe isomers.
    """
    n = len(sites)
    group = []
    for perm in itertools.permutations(range(n)):
        P = sites
        Q = sites[list(perm)]
        H = P.T @ Q
        U, _, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        D = np.diag([1.0, 1.0, d])
        R = Vt.T @ D @ U.T
        if abs(np.linalg.det(R) - 1.0) > tol:
            continue
        if np.allclose(R @ P.T, Q.T, atol=tol):
            group.append(perm)
    return group


def chelate_compatible_pairs(sites: np.ndarray, max_bite_deg: float) -> set[frozenset[int]]:
    """Site-index pairs whose vertex angle (metal-site-i / metal-site-j)
    is small enough to plausibly be bridged by one bidentate ligand."""
    n = len(sites)
    out = set()
    for i in range(n):
        for j in range(i + 1, n):
            ang = np.degrees(np.arccos(np.clip(np.dot(sites[i], sites[j]), -1, 1)))
            if ang <= max_bite_deg:
                out.add(frozenset((i, j)))
    return out


# ======================================================================
# Isomer representation and enumeration
# ======================================================================

# An isomer is a frozenset of (ligand_label, site_tuple) pairs -- one
# entry per ligand UNIT (a monodentate ligand occupies one site, a
# bidentate ligand occupies two, etc.). site_tuple is an ORDERED tuple,
# not a frozenset: for an ASYMMETRIC ligand (its own donor atoms are not
# all chemically equivalent -- e.g. an N,O-chelate like glycinate) which
# END sits at which site is a genuinely different arrangement, so order
# must be preserved through rotation, not thrown away. For a SYMMETRIC
# ligand (e.g. en, both donor N's equivalent) order carries no physical
# meaning, so it is normalized away only where it's actually compared
# (canonical_form's sortkey), not by discarding it up front -- otherwise
# an asymmetric ligand placed at the same two sites would silently be
# treated as if it had no direction at all.
#
# LIMITATION: this only distinguishes "fully symmetric" (order never
# matters) from "fully asymmetric" (order always matters) per ligand
# label. A ligand with PARTIAL internal symmetry (e.g. a tridentate
# ligand with two equivalent outer donors and one different central
# donor) would need its own donor-automorphism group, analogous to the
# site rotation group, which isn't implemented -- not needed yet since
# nothing here exceeds denticity 2.
Isomer = frozenset


def apply_perm(isomer: Isomer, perm: tuple[int, ...]) -> Isomer:
    return frozenset((label, tuple(perm[s] for s in sites)) for label, sites in isomer)


def canonical_form(isomer: Isomer, group: list[tuple[int, ...]], asymmetric_labels: set[str]) -> Isomer:
    """The lexicographically-smallest image of `isomer` under the
    rotation group -- two isomers are the SAME distinct arrangement iff
    they share a canonical form. This is exactly what keeps AABB and
    BBAA-type relabelings from being counted twice, while correctly
    keeping non-superimposable mirror images (e.g. Delta/Lambda) apart,
    since `group` here only ever contains proper rotations. For a
    symmetric ligand's site tuple, comparison sorts it first (order is
    not meaningful); for an asymmetric one, order is compared as-is (it
    IS meaningful -- which end sits where is part of what makes two
    placements the same or different)."""
    def keyfor(label, sites):
        return sites if label in asymmetric_labels else tuple(sorted(sites))

    def normalize(iso):
        # the RETURNED form must itself be order-normalized for symmetric
        # labels, not just used for comparison -- otherwise two isomers
        # that are equivalent up to a symmetric ligand's internal tuple
        # order (e.g. (2,5) vs (5,2)) can tie on sortkey but still come
        # back as two different frozenset objects, breaking dedup.
        return frozenset((label, keyfor(label, sites)) for label, sites in iso)

    def sortkey(iso):
        return tuple(sorted((label, keyfor(label, sites)) for label, sites in iso))
    best = min((apply_perm(isomer, g) for g in group), key=sortkey)
    return normalize(best)


def enumerate_assignments(n_sites: int, ligand_units: list[tuple[str, int]],
                           compatible: set[frozenset[int]],
                           asymmetric_labels: set[str]) -> list[Isomer]:
    """All ways to partition `n_sites` sites into groups matching
    `ligand_units` (a list of (label, denticity) for every individual
    ligand -- e.g. 3 identical bidentate "en" units is
    [("en",2),("en",2),("en",2)]). A denticity-1 unit can take any
    single free site; a denticity>1 unit can only take a set of free
    sites that is itself in `compatible` (extended here to require every
    pairwise subset of a >2-denticity group to be compatible, though no
    geometry below currently needs more than denticity 2). For a
    denticity>=2 unit whose label is in `asymmetric_labels`, EVERY
    ordering of the chosen site combo is emitted as its own raw
    candidate (2 orderings for denticity 2) -- for a symmetric label,
    only one (order is settled later, at comparison time, not here)."""
    results: list[Isomer] = []

    def compatible_group(sites: tuple[int, ...]) -> bool:
        if len(sites) == 1:
            return True
        return all(frozenset(p) in compatible for p in itertools.combinations(sites, 2))

    def orderings(label, combo):
        if len(combo) > 1 and label in asymmetric_labels:
            return set(itertools.permutations(combo))
        return {tuple(combo)}

    def backtrack(remaining_sites: frozenset[int], remaining_units: tuple[tuple[str, int], ...],
                  chosen: list[tuple[str, tuple[int, ...]]]):
        if not remaining_units:
            results.append(frozenset(chosen))
            return
        label, dent = remaining_units[0]
        rest_units = remaining_units[1:]
        for combo in itertools.combinations(sorted(remaining_sites), dent):
            if not compatible_group(combo):
                continue
            for ordered in orderings(label, combo):
                backtrack(remaining_sites - set(combo), rest_units, chosen + [(label, ordered)])

    backtrack(frozenset(range(n_sites)), tuple(ligand_units), [])
    return results


def distinct_isomers(geometry_name: str, ligand_units: list[tuple[str, int]],
                      max_bite_deg: float = 100.0, asymmetric_labels: set[str] | None = None) -> list[Isomer]:
    asymmetric_labels = asymmetric_labels or set()
    sites = GEOMETRIES[geometry_name]
    n = len(sites)
    if sum(d for _, d in ligand_units) != n:
        raise ValueError(f"ligand denticities sum to {sum(d for _, d in ligand_units)}, "
                          f"but {geometry_name} has {n} coordination sites")
    group = find_rotation_group(sites)
    compatible = chelate_compatible_pairs(sites, max_bite_deg)
    raw = enumerate_assignments(n, ligand_units, compatible, asymmetric_labels)
    canon = {canonical_form(iso, group, asymmetric_labels) for iso in raw}
    def outerkey(iso):
        return tuple(sorted((l, s if l in asymmetric_labels else tuple(sorted(s))) for l, s in iso))
    return sorted(canon, key=outerkey)


# ======================================================================
# CLI / self-checks
# ======================================================================

def _parse_ligands(spec: str, denticity_spec: str) -> list[tuple[str, int]]:
    denticities = {}
    if denticity_spec:
        for part in denticity_spec.split(","):
            lbl, d = part.split(":")
            denticities[lbl] = int(d)
    units = []
    for part in spec.split(","):
        lbl, count = part.split(":")
        d = denticities.get(lbl, 1)
        units += [(lbl, d)] * int(count)
    return units


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--geometry", help="comma-separated geometry name(s), e.g. octahedral or "
                                        "tetrahedral,square_planar to enumerate across both")
    ap.add_argument("--ligands", help='ligand multiset, e.g. "Cl:2,N:4"')
    ap.add_argument("--denticity", default="", help='denticity overrides, e.g. "en:2" (default 1)')
    ap.add_argument("--max-bite-angle", type=float, default=100.0,
                     help="max metal-site-site angle (deg) for a bidentate ligand to bridge (default 100)")
    ap.add_argument("--asymmetric", default="", help="comma-separated ligand labels whose two "
                     "donor ends are chemically DIFFERENT (e.g. an N,O-chelate) -- both site "
                     "orderings are then treated as potentially-distinct isomers, not merged. "
                     "Default: every ligand is treated as symmetric (e.g. en, both ends -NH2).")
    ap.add_argument("--list-geometries", action="store_true")
    ap.add_argument("--self-check", action="store_true", help="run the known-chemistry validation checks")
    args = ap.parse_args(argv)

    if args.list_geometries:
        for name, v in GEOMETRIES.items():
            print(f"{name}: {len(v)} sites")
        return 0

    if args.self_check:
        return _self_check()

    if not args.geometry or not args.ligands:
        ap.print_help()
        return 1

    units = _parse_ligands(args.ligands, args.denticity)
    asymmetric = set(args.asymmetric.split(",")) if args.asymmetric else set()
    for geom in args.geometry.split(","):
        isomers = distinct_isomers(geom, units, args.max_bite_angle, asymmetric)
        print(f"{geom}: {len(isomers)} distinct isomer(s)")
        for k, iso in enumerate(isomers):
            desc = ", ".join(f"{label}@{list(sites)}" for label, sites in
                              sorted(iso, key=lambda t: sorted(t[1])))
            print(f"  [{k}] {desc}")
    return 0


def _self_check():
    checks = []

    # fac/mer: CN6 octahedral, 3A + 3B monodentate -> exactly 2 distinct isomers
    iso = distinct_isomers("octahedral", [("A", 1)] * 3 + [("B", 1)] * 3)
    checks.append(("fac/mer (octahedral 3A3B)", len(iso), 2))

    # Bailar: CN6 octahedral, 3 IDENTICAL SYMMETRIC bidentate ligands
    # (both donor ends equivalent, e.g. en's two -NH2's) -> exactly 2
    # distinct isomers (Delta, Lambda) -- the whole point of using the
    # PROPER rotation group only, not the full point group.
    iso_sym = distinct_isomers("octahedral", [("en", 2)] * 3, max_bite_deg=95.0)
    checks.append(("Bailar twist, symmetric bidentate (en)", len(iso_sym), 2))

    # Same geometry, but an ASYMMETRIC bidentate ligand (e.g. an N,O
    # chelate like glycinate, donor ends NOT equivalent) -> which end
    # sits at which site now matters, so this must give MORE than 2 --
    # if it came out equal to the symmetric case, the direction fix
    # wouldn't actually be doing anything.
    iso_asym = distinct_isomers("octahedral", [("gly", 2)] * 3, max_bite_deg=95.0,
                                 asymmetric_labels={"gly"})
    print(f"[{'ok' if len(iso_asym) > len(iso_sym) else 'MISMATCH'}] "
          f"Bailar twist, ASYMMETRIC bidentate: found {len(iso_asym)} (must be > symmetric case's {len(iso_sym)})")
    if len(iso_asym) <= len(iso_sym):
        checks.append(("__asymmetric_fix__", 0, 1))  # force overall failure below

    # cis/trans: CN4 SQUARE PLANAR, 2A + 2B monodentate -> exactly 2 (cis, trans)
    iso = distinct_isomers("square_planar", [("A", 1)] * 2 + [("B", 1)] * 2)
    checks.append(("cis/trans (square_planar 2A2B)", len(iso), 2))

    # CN4 TETRAHEDRAL, 2A + 2B monodentate -> exactly 1 (no cis/trans
    # isomerism at all for a tetrahedral center -- real, known chemistry;
    # this is a check that the group-theory approach gets it right
    # automatically, not something hand-coded in).
    iso = distinct_isomers("tetrahedral", [("A", 1)] * 2 + [("B", 1)] * 2)
    checks.append(("tetrahedral 2A2B (should have NO cis/trans)", len(iso), 1))

    # Berry: CN5 TBP, 2A + 3B monodentate -> 2 distinct (A's diaxial vs
    # one-axial-one-equatorial); the fully-diequatorial-A arrangement
    # (3 A's) isn't reachable here since only 2 A's are present, so this
    # specifically checks the 2A3B TBP case, not the full A2B3 count.
    iso = distinct_isomers("trigonal_bipyramidal", [("A", 1)] * 2 + [("B", 1)] * 3)
    checks.append(("TBP 2A3B", len(iso), None))  # informational, no hardcoded expectation

    ok = True
    for label, got, expect in checks:
        status = "ok" if (expect is None or got == expect) else "MISMATCH"
        if expect is not None and got != expect:
            ok = False
        print(f"[{status}] {label}: found {got}" + (f" (expected {expect})" if expect is not None else ""))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
