"""Single-molecule coordinate templates pulled from an existing structure.

A template is one residue's atoms, centred on its geometric centre and rotated
so the long molecular axis lies along z.  Taking the conformation from a real
equilibrated system means the molecule we stamp into the core is already
physical -- no internal geometry is ever generated here.
"""

import numpy as np
import MDAnalysis as mda


def _principal_axis(positions):
    """Unit eigenvector of the largest moment of a coordinate set."""
    centred = positions - positions.mean(axis=0)
    evals, evecs = np.linalg.eigh(np.cov(centred.T))
    return evecs[:, np.argmax(evals)]


def _rotation_between(a, b):
    """Rotation matrix mapping unit vector ``a`` onto unit vector ``b``."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


class Template:
    """One molecule's atoms, ready to be stamped into the core."""

    def __init__(self, resname, names, positions):
        self.resname = resname
        self.names = list(names)
        self.positions = np.asarray(positions, dtype=float)

    @property
    def n_atoms(self):
        return len(self.names)

    @property
    def span(self):
        """Bounding-box size (x, y, z) in angstrom."""
        return self.positions.max(axis=0) - self.positions.min(axis=0)

    def randomly_rotated(self, rng):
        """Coordinates under a uniformly random 3-D rotation."""
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        w, x, y, z = q
        rot = np.array([
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ])
        return self.positions @ rot.T


def extract_template(reference, resname, max_candidates=200):
    """Return a :class:`Template` for ``resname`` taken from ``reference``.

    Every copy of the residue (up to ``max_candidates``) is stood upright and
    the most compact one wins, so grid placement starts with the fewest clashes.
    """
    u = reference if isinstance(reference, mda.Universe) else mda.Universe(reference)
    sel = u.select_atoms(f"resname {resname}")
    if not len(sel):
        raise ValueError(f"reference structure has no residue named {resname!r}")

    best = None
    for res in sel.residues[:max_candidates]:
        ag = res.atoms
        pos = ag.positions - ag.positions.mean(axis=0)
        pos = pos @ _rotation_between(_principal_axis(pos), np.array([0.0, 0.0, 1.0])).T
        footprint = float((pos[:, :2].max(axis=0) - pos[:, :2].min(axis=0)).max())
        if best is None or footprint < best[0]:
            best = (footprint, ag.names, pos)

    _, names, pos = best
    return Template(resname, names, pos)
