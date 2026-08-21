"""ld_lite -- build a lipid droplet from a bilayer and equilibrate it.

A neutral-lipid (triolein) core inserted between the two leaflets of a solvated
CHARMM-GUI bilayer, plus the standard CHARMM-GUI minimization and six-stage
equilibration to settle it.

    >>> from ld_lite import build_ld, write_top, write_index
    >>> ld = build_ld("bilayer.gro", n_trio=150)
    >>> ld.atoms.write("system.gro")
    >>> write_top(ld, "topol.top")
    >>> write_index(ld, "index.ndx")

or from the command line::

    python -m ld_lite build --bilayer bil.gro --n-trio 150 --out runs/ld1
    python -m ld_lite equil runs/ld1
"""

__version__ = "0.1.0"

from .build import build_ld, make_whole
from .gmx import molecule_counts, write_top, write_index
from .templates import Template, extract_template
from .relax import relax_clashes

__all__ = [
    "build_ld",
    "make_whole",
    "molecule_counts",
    "write_top",
    "write_index",
    "Template",
    "extract_template",
    "relax_clashes",
]
