#!/usr/bin/env python3
"""generate_isomer_structures.py -- given a real starting structure for a
metal complex, detect its ligand fragments and their denticity from
connectivity, enumerate the symmetry-distinct ways those SAME fragments
could be arranged on the coordination polyhedron (via isomer_enumerate),
and for every isomer other than the input's own, rigidly transplant each
ligand fragment (whole fragment, not just its donor atom) onto its new
target site direction -- producing a full set of starting structures for
crest's nci_metal runtype to launch one MTD run from each.

This is the piece that turns an abstract isomer_enumerate.py result
(which site gets which ligand LABEL) into an actual new .xyz file: the
ligand's own internal geometry (backbone, H's, ...) is preserved exactly,
only rigidly rotated+translated so its donor atom(s) land on the new
target direction(s) instead of the old one(s).

Usage:
    python generate_isomer_structures.py input.xyz --metal 1 --outdir isomers/
"""
from __future__ import annotations

import argparse
import itertools
import os

import numpy as np

from isomer_enumerate import GEOMETRIES, distinct_isomers, find_rotation_group, canonical_form

COVALENT_RADIUS = {  # Angstrom, rough single-bond covalent radii
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
}
# common first-row-ish metals; extend as needed
METAL_RADIUS = {"Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Pt": 1.36,
                "Fe": 1.32, "Mn": 1.39, "Cr": 1.39, "Pd": 1.39}


def read_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0])
    els, coords = [], []
    for l in lines[2:2 + n]:
        p = l.split()
        els.append(p[0])
        coords.append([float(p[1]), float(p[2]), float(p[3])])
    return els, np.array(coords)


def write_xyz(path, els, coords, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(els)}\n{comment}\n")
        for el, c in zip(els, coords):
            f.write(f"{el:<2s}  {c[0]:15.8f} {c[1]:15.8f} {c[2]:15.8f}\n")


def radius(el):
    return METAL_RADIUS.get(el, COVALENT_RADIUS.get(el, 0.75))


def detect_bonds(els, coords, scale=1.2):
    n = len(els)
    bonds = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(i + 1, n):
            cutoff = scale * (radius(els[i]) + radius(els[j]))
            if np.linalg.norm(coords[i] - coords[j]) <= cutoff:
                bonds[i, j] = bonds[j, i] = True
    return bonds


def detect_fragments(els, coords, metal_idx, bonds):
    """One fragment per ligand: BFS from each metal-neighbor, excluding
    the metal itself, over the bond graph. Returns a list of dicts with
    'atoms' (indices, metal excluded), 'donors' (the subset directly
    bonded to the metal), each fragment visited exactly once."""
    n = len(els)
    assigned = {metal_idx}
    fragments = []
    for start in range(n):
        if start == metal_idx or not bonds[metal_idx, start] or start in assigned:
            continue
        # BFS over the whole molecule graph EXCLUDING the metal node
        seen = {start}
        queue = [start]
        while queue:
            cur = queue.pop()
            for nxt in range(n):
                if nxt == metal_idx or nxt in seen or not bonds[cur, nxt]:
                    continue
                seen.add(nxt)
                queue.append(nxt)
        donors = sorted(a for a in seen if bonds[metal_idx, a])
        fragments.append({"atoms": sorted(seen), "donors": donors})
        assigned |= seen
    return fragments


def fragment_formula(els, atoms):
    counts = {}
    for a in atoms:
        counts[els[a]] = counts.get(els[a], 0) + 1
    return "".join(f"{el}{counts[el]}" for el in sorted(counts))


def color_refine(atoms, bonds, els, rounds=None):
    """1-WL-style color refinement over the bonding graph restricted to
    `atoms` (element labels + connectivity only -- no coordinates at
    all). Two atoms end up with the SAME final color only if no local
    structural feature, out to `rounds` bonds away, distinguishes them;
    this is invariant under any graph automorphism that respects
    element labels and connectivity. That's exactly what's needed here:
    (a) whether two ligand FRAGMENTS are the same type (compare their
    sorted color multisets) and (b) whether a fragment's own donor
    atoms are chemically interchangeable (compare their individual
    final colors) -- both from topology alone, so ordinary geometric
    noise in a real (as opposed to idealized) structure never affects
    the answer, unlike comparing 3D shapes directly would. Not a fully
    rigorous graph-isomorphism test in every pathological case (1-WL
    famously fails to distinguish a few contrived regular graphs), but
    it's the standard, practical technique behind canonical molecular
    identifiers (e.g. InChI's own canonicalization uses the same idea)
    and is more than sufficient for real ligands."""
    atoms = list(atoms)
    if rounds is None:
        rounds = len(atoms)
    color = {a: els[a] for a in atoms}
    for _ in range(rounds):
        new_sig = {a: (color[a], tuple(sorted(color[b] for b in atoms if bonds[a, b]))) for a in atoms}
        uniq = sorted(set(new_sig.values()), key=str)
        remap = {sig: i for i, sig in enumerate(uniq)}
        color = {a: remap[new_sig[a]] for a in atoms}
    return color


def fragment_topology_label(els, atoms, bonds):
    """A hashable value identifying this fragment's TYPE by bonding
    topology (not just elemental formula, which would wrongly conflate
    two constitutionally different ligands that happen to share an
    atom count -- e.g. two different C2H8N2 connectivities). Pairs the
    human-readable formula with the exact color-refinement signature,
    so two fragments compare equal as the same ligand type iff their
    bonding graphs (element-labelled) are indistinguishable, not just
    their formulas."""
    color = color_refine(atoms, bonds, els)
    return fragment_formula(els, atoms), tuple(sorted(color.values()))


def best_matching_geometry(donor_dirs: np.ndarray) -> str:
    """Which tabulated GEOMETRIES entry (same site count) best matches
    the CURRENT donor-direction set, by comparing sorted pairwise
    angles -- used to pick which idealized site-vector set to treat the
    real structure's positions as an instance of."""
    n = len(donor_dirs)
    def pairwise_angles(v):
        out = []
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                out.append(np.degrees(np.arccos(np.clip(np.dot(v[i], v[j]), -1, 1))))
        return sorted(out)
    target = pairwise_angles(donor_dirs)
    best, best_err = None, None
    for name, sites in GEOMETRIES.items():
        if len(sites) != n:
            continue
        err = np.sum((np.array(pairwise_angles(sites)) - np.array(target)) ** 2)
        if best_err is None or err < best_err:
            best, best_err = name, err
    return best


def rotation_matrix(axis, angle):
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)


def align_vectors(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """The minimal proper rotation mapping unit vector a exactly onto
    unit vector b. Handles the antiparallel case explicitly: a 180-deg
    rotation must be about an axis PERPENDICULAR to a, not about a
    itself (rotating about axis a leaves a fixed -- using a as the axis
    would map a->a, not a->-a)."""
    axis = np.cross(a, b)
    if np.linalg.norm(axis) < 1e-8:
        if np.dot(a, b) > 0:
            return np.eye(3)
        perp = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
        perp = perp - a * np.dot(perp, a)
        perp /= np.linalg.norm(perp)
        return 2 * np.outer(perp, perp) - np.eye(3)
    angle = np.arccos(np.clip(np.dot(a, b), -1, 1))
    return rotation_matrix(axis, angle)


def exact_two_point_rotation(a1, b1, a2, b2) -> np.ndarray:
    """The proper rotation mapping unit vectors a1->b1 AND a2->b2
    exactly, constructed directly rather than via SVD/Kabsch. An exact
    fit is guaranteed to exist here (old and new donor-direction pairs
    both come from the same rigid, idealized polyhedron, so their
    pairwise angle always matches) and is generically UNIQUE (a proper
    rotation is fully determined by its action on two non-collinear
    vectors) -- so this is solved as two robust single-vector
    alignments rather than trusting SVD's proper-rotation sign
    correction in the rank-deficient (only 2 points) regime, where it
    can be numerically unstable (this is what actually produced an
    8.6-deg-off placement via kabsch_rotation during testing).
    Step 1: align a1 -> b1. Step 2: the remaining freedom is a rotation
    about b1 (which leaves b1 fixed) -- resolve it by aligning the
    component of (rotated a2) perpendicular to b1 with the component of
    b2 perpendicular to b1; since both perpendicular components are, by
    construction, orthogonal to b1, their cross product is automatically
    parallel to b1, so this second alignment is automatically a pure
    rotation about b1 and never disturbs the first step's result."""
    R1 = align_vectors(a1, b1)
    a2p = R1 @ a2
    p1 = a2p - b1 * np.dot(a2p, b1)
    p2 = b2 - b1 * np.dot(b2, b1)
    n1, n2 = np.linalg.norm(p1), np.linalg.norm(p2)
    if n1 < 1e-8 or n2 < 1e-8:
        return R1
    p1, p2 = p1 / n1, p2 / n2
    if np.linalg.norm(np.cross(p1, p2)) < 1e-8 and np.dot(p1, p2) < 0:
        # p1, p2 antiparallel: align_vectors' generic fallback would
        # pick SOME axis perpendicular to p1, but that axis is only
        # guaranteed to fix p1->p2 in isolation -- it is NOT guaranteed
        # to leave b1 fixed unless it IS b1. The one rotation that
        # fixes this second donor while never disturbing the first
        # step's (already-correct) b1 alignment is specifically a
        # 180-deg rotation about b1 itself (b1 stays fixed under any
        # rotation about b1, and rotating 180 about b1 sends any vector
        # perpendicular to it to its own antipode within that plane --
        # exactly p1 -> p2 here).
        R2 = rotation_matrix(b1, np.pi)
    else:
        R2 = align_vectors(p1, p2)
    return R2 @ R1


def clash_count(els, coords, group_atoms, other_atoms, slack=0.75):
    """Number of (group, other) atom pairs closer than `slack` times
    their covalent-radius sum -- a real, if crude, steric overlap check.
    Only meant to catch genuine clashes (two atoms nearly on top of each
    other after a rigid transplant), not ordinary close contacts."""
    n = 0
    for a in group_atoms:
        for b in other_atoms:
            cutoff = slack * (COVALENT_RADIUS.get(els[a], 0.75) + COVALENT_RADIUS.get(els[b], 0.75))
            if np.linalg.norm(coords[a] - coords[b]) < cutoff:
                n += 1
    return n


def transplant_fragment(els, coords, metal_pos, fragment, old_donor_dirs, new_donor_dirs,
                         other_atoms=(), expand_by=0.5):
    """Rigidly rotate the WHOLE fragment (all its atoms, not just the
    donor(s)) about the metal so its donor(s) move from old_donor_dirs
    to new_donor_dirs, preserving the fragment's own internal geometry
    exactly, then translate the WHOLE fragment (one shared offset, so
    the move stays rigid) further out by `expand_by` (Angstrom) along
    the average new-donor direction. A pure rotation+translation keeps
    every internal distance in the fragment exactly what it was; giving
    different atoms different offsets is NOT rigid and will distort
    the ligand (confirmed by testing: an earlier version that pushed
    each donor along its own direction while other atoms used the
    shared average visibly stretched/compressed bonds, e.g. N-H
    1.01->0.801 A). The metal-donor bond ends up deliberately too LONG
    by `expand_by` for a multi-donor fragment (the shared average isn't
    parallel to any single donor's own direction) -- that's intentional
    and expected, not a bug: nci_metal's own initial geometry
    optimization, combined with restrain=true (the UFF4MOF flat-bottom
    distance restraint already built for coordination-sphere
    protection), pulls it back to the correct bond length and direction
    as the very first thing that happens, letting the force field
    resolve any remaining close contacts along the way rather than a
    hand-coded search having to get it exactly right up front.

    Returns (new_coords, R, donor_rotation_only): R is the rotation
    used, and donor_rotation_only maps each donor atom index to its
    position BEFORE the expansion offset (i.e. the exact rigid-rotation
    result) -- the caller validates against this, not the final
    (deliberately offset) position, since only the rotation is
    supposed to be exact.

    Denticity 1 leaves a full free rotation about the bond axis
    undetermined (chemically meaningless for the fragment alone, but
    NOT for the rest of the molecule, since that twist decides where
    every other atom in a bulky ligand swings to). Rather than picking
    one twist arbitrarily, sample a set of candidate twists and keep
    whichever has the fewest steric clashes against `other_atoms` (all
    atoms already placed elsewhere in the growing structure) -- not
    guaranteed clash-free for a sufficiently bulky ligand, but a real
    search over the one DOF actually available, and still checked
    afterwards (see clash_count in the caller) rather than hidden.
    Denticity >=2 has no such free DOF: exact_two_point_rotation fully
    determines R, with no search needed."""
    n_donors = len(old_donor_dirs)
    if n_donors == 1:
        a, b = old_donor_dirs[0], new_donor_dirs[0]
        R0 = align_vectors(a, b)
        donor_idx = fragment["donors"][0]
        rotation_only = {donor_idx: metal_pos + R0 @ (coords[donor_idx] - metal_pos)}

        def apply(R):
            out = coords.copy()
            for at in fragment["atoms"]:
                out[at] = metal_pos + expand_by * b + R @ (coords[at] - metal_pos)
            return out

        if not other_atoms:
            return apply(R0), R0, rotation_only

        best_coords, best_clashes = None, None
        for twist_deg in range(0, 360, 15):
            Rtwist = rotation_matrix(b, np.radians(twist_deg))
            candidate = apply(Rtwist @ R0)
            n_clash = clash_count(els, candidate, fragment["atoms"], other_atoms)
            if best_clashes is None or n_clash < best_clashes:
                best_clashes, best_coords = n_clash, candidate
            if best_clashes == 0:
                break
        return best_coords, R0, rotation_only

    R = exact_two_point_rotation(old_donor_dirs[0], new_donor_dirs[0], old_donor_dirs[1], new_donor_dirs[1])
    rotation_only = {d: metal_pos + R @ (coords[d] - metal_pos) for d in fragment["donors"]}
    outward = np.mean(new_donor_dirs, axis=0)
    outward /= np.linalg.norm(outward)
    new_coords = coords.copy()
    for a in fragment["atoms"]:
        rel = coords[a] - metal_pos
        new_coords[a] = metal_pos + expand_by * outward + R @ rel
    return new_coords, R, rotation_only


def generate(input_xyz: str, metal_idx: int, outdir: str, max_bite_angle: float = 100.0):
    els, coords = read_xyz(input_xyz)
    metal_pos = coords[metal_idx]
    bonds = detect_bonds(els, coords)
    fragments = detect_fragments(els, coords, metal_idx, bonds)

    # label each fragment by bonding TOPOLOGY, not just elemental
    # formula -- two fragments are the same ligand type only if their
    # (element-labelled) bonding graphs are indistinguishable, which
    # correctly tells apart two constitutionally different ligands that
    # happen to share an atom count, something formula-only matching
    # would silently conflate. Likewise, whether a bidentate fragment's
    # two donors are chemically interchangeable (symmetric) is now
    # detected directly from its own bonding graph -- do the two donors
    # get the SAME color-refinement color? -- rather than a hardcoded
    # per-formula whitelist that would need updating for every new
    # ligand. Both checks use topology alone, never 3D coordinates, so
    # ordinary geometric noise in a real (non-idealized) structure can
    # never affect the answer.
    for frag in fragments:
        frag["label"] = fragment_topology_label(els, frag["atoms"], bonds)
        frag["formula"] = frag["label"][0]
        frag["denticity"] = len(frag["donors"])
        if frag["denticity"] >= 2:
            donor_colors = color_refine(frag["atoms"], bonds, els)
            frag["asymmetric"] = len({donor_colors[d] for d in frag["donors"]}) > 1
        else:
            frag["asymmetric"] = False

    all_donor_idx = [d for frag in fragments for d in frag["donors"]]
    donor_dirs = np.array([(coords[d] - metal_pos) / np.linalg.norm(coords[d] - metal_pos)
                            for d in all_donor_idx])
    geom_name = best_matching_geometry(donor_dirs)
    geom_sites = GEOMETRIES[geom_name]
    print(f"detected {len(fragments)} ligand fragment(s), best-fit geometry: {geom_name}")
    for frag in fragments:
        print(f"  {frag['formula']} (denticity {frag['denticity']}, "
              f"{'asymmetric' if frag['asymmetric'] else 'symmetric'}) at atoms {frag['atoms']}")

    # match each donor atom to its nearest idealized site index (by angle)
    donor_to_site = {}
    used_sites = set()
    for d, ddir in zip(all_donor_idx, donor_dirs):
        dists = sorted(range(len(geom_sites)), key=lambda s: -np.dot(ddir, geom_sites[s]))
        for s in dists:
            if s not in used_sites:
                donor_to_site[d] = s
                used_sites.add(s)
                break

    ligand_units = [(frag["label"], frag["denticity"]) for frag in fragments]
    asymmetric_labels = {frag["label"] for frag in fragments if frag["asymmetric"]}
    isomers = distinct_isomers(geom_name, ligand_units, max_bite_angle, asymmetric_labels)
    print(f"{len(isomers)} distinct isomer(s) for this ligand set on {geom_name}")

    # figure out which of the enumerated isomers matches the INPUT's own
    # current arrangement, so we can skip writing that one back out.
    # Compare CANONICAL forms directly -- comparing fragments to target
    # groups by iteration order (as an earlier version of this script
    # did) is wrong, because which specific same-label fragment gets
    # matched to which specific same-label target group is arbitrary and
    # doesn't correspond to anything geometrically meaningful; only the
    # canonical form is a well-defined thing to compare against.
    rotation_group = find_rotation_group(geom_sites)
    input_raw = frozenset(
        (frag["label"], tuple(donor_to_site[d] for d in frag["donors"]))
        for frag in fragments
    )
    input_canonical = canonical_form(input_raw, rotation_group, asymmetric_labels)

    os.makedirs(outdir, exist_ok=True)
    written = 0
    for k, iso in enumerate(isomers):
        if iso == input_canonical:
            print(f"  [{k}] matches the input structure's own arrangement -- skipped")
            continue
        # map: label -> list of target site-groups for this isomer. Which
        # specific same-label fragment gets sent to which specific
        # same-label target group is arbitrary and doesn't matter --
        # same-label fragments are chemically interchangeable, so any
        # bijection produces the identical resulting molecule.
        by_label = {}
        for label, sites in iso:
            by_label.setdefault(label, []).append(sites)

        new_coords = coords.copy()
        consumed = {label: 0 for label in by_label}
        for frag in fragments:
            label = frag["label"]
            target_sites = by_label[label][consumed[label]]
            consumed[label] += 1
            old_dirs = [(coords[d] - metal_pos) / np.linalg.norm(coords[d] - metal_pos) for d in frag["donors"]]
            if frag["denticity"] == 1:
                new_dirs = [geom_sites[target_sites[0]]]
            else:
                new_dirs = [geom_sites[s] for s in target_sites]
            other_atoms = [a for a in range(len(els)) if a != metal_idx and a not in frag["atoms"]]
            # adaptive expansion: try the smallest expand_by that avoids a
            # clash for THIS fragment against the structure's current state,
            # rather than always using the full default (0.5 A). A large,
            # unconditional expansion pushes the donor-metal distance past
            # not just the clash it was meant to avoid, but potentially
            # past GFN-FF's OWN bond-perception cutoff too -- confirmed
            # directly: a simple monodentate Cl/NH3 swap (no real clash
            # risk) expanded to Cu-Cl 2.8 A / Cu-N 2.52 A caused GFN-FF to
            # see the whole structure as 5 disconnected fragments ("#
            # fragment in coord: 5"), which no restraint can fix after the
            # fact (restrain=true is an added biasing force, not a change
            # to GFN-FF's own once-only topology perception). Starting
            # from expand_by=0 and growing only as far as actually needed
            # keeps every donor at or near its real bonding distance
            # whenever no clash exists to avoid in the first place.
            for trial_expand in (0.0, 0.15, 0.3, 0.5):
                trial_coords, _, trial_rotation_only = transplant_fragment(
                    els, new_coords, metal_pos, frag, old_dirs, new_dirs,
                    other_atoms=other_atoms, expand_by=trial_expand)
                if clash_count(els, trial_coords, frag["atoms"], other_atoms) == 0:
                    break
            new_coords, rotation_only = trial_coords, trial_rotation_only
            # hard check: the ROTATION must have put every donor exactly
            # on its intended target direction -- checked on the
            # pre-expansion position (rotation_only), not the final one,
            # since the deliberate expansion offset is only guaranteed
            # to keep a SINGLE donor (denticity 1) exactly aligned; for
            # denticity >=2 it intentionally does not (see
            # transplant_fragment's docstring). This must never fail on
            # the rotation itself -- if it does, something upstream
            # (site assignment, or a numerical-precision issue in the
            # rotation-group finder for a geometry with real symmetry)
            # produced an inconsistent target, and writing the structure
            # anyway would silently hand crest a broken starting point
            # instead of surfacing the bug.
            for donor, target_site in zip(frag["donors"], target_sites):
                rel = rotation_only[donor] - metal_pos
                got_dir = rel / np.linalg.norm(rel)
                want_dir = geom_sites[target_site]
                if not np.allclose(got_dir, want_dir, atol=1e-2):
                    raise RuntimeError(
                        f"isomer {k}: donor atom {donor} (fragment {frag['formula']}) rotation-only position "
                        f"landed at direction {got_dir}, not its intended target site {target_site} "
                        f"({want_dir}) -- aborting rather than writing a structure with an "
                        f"inconsistent placement")

        # final, whole-structure clash report -- catches anything the
        # per-fragment twist search above couldn't avoid (always true for
        # denticity>=2 fragments, which have no free DOF to search over;
        # possible even for denticity 1 if every sampled twist still
        # clashes with something, e.g. a genuinely bulky ligand jammed
        # into a tight coordination environment).
        n_clash = 0
        for frag in fragments:
            other_atoms = [a for a in range(len(els)) if a != metal_idx and a not in frag["atoms"]]
            n_clash += clash_count(els, new_coords, frag["atoms"], other_atoms)
        n_clash //= 2  # each cross-fragment pair gets counted from both sides above
        if n_clash > 0:
            print(f"  [{k}] **WARNING** {n_clash} steric clash(es) detected in this generated "
                  f"structure -- inspect before using it as an MTD starting point (a bulky/flexible "
                  f"ligand's one free twist wasn't enough to avoid every overlap here)")

        outpath = os.path.join(outdir, f"isomer_{written}.xyz")
        write_xyz(outpath, els, new_coords, comment=f"generated isomer {k} of {input_xyz}")
        print(f"  [{k}] -> {outpath}")
        written += 1
    print(f"wrote {written} new starting structure(s) (of {len(isomers)} total distinct isomers)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_xyz")
    ap.add_argument("--metal", type=int, required=True, help="1-based index of the metal atom")
    ap.add_argument("--outdir", default="isomers")
    ap.add_argument("--max-bite-angle", type=float, default=100.0)
    args = ap.parse_args(argv)
    generate(args.input_xyz, args.metal - 1, args.outdir, args.max_bite_angle)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
