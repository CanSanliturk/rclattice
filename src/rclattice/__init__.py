"""rclattice — lattice modelling of reinforced concrete, analysed with OpenSees.

The Problem (problem.py), the FE Model (model.py), the mesher (mesh.py) and the builders
(builders.py) are backend-agnostic and do NOT import openseespy. All OpenSees calls live in
opensees.py.
"""

from .builders import build_continuum, build_lattice, build_lattice_rc, select_nodes
from .calibration import CalibrationTargets, calibrate_lattice, continuum_targets
from .materials import (
    concrete_nd_elastic,
    concrete_uniaxial_elastic,
    concrete_uniaxial_nonlinear,
    concrete_uniaxial_regularized,
    steel_uniaxial,
)
from .model import Element, Load, Model, NDMaterial, Node, Support, UniaxialMaterial
from .problem import (
    BoxLoad,
    BoxSupport,
    CompoundRectangles,
    ConcreteGrade,
    EdgeLoad,
    EdgeSupport,
    Problem,
    Rebar,
    RectangleDomain,
    SteelGrade,
    portal_frame,
)

__all__ = [
    "Problem",
    "ConcreteGrade",
    "RectangleDomain",
    "CompoundRectangles",
    "portal_frame",
    "EdgeSupport",
    "EdgeLoad",
    "BoxSupport",
    "BoxLoad",
    "build_lattice",
    "build_continuum",
    "select_nodes",
    "calibrate_lattice",
    "continuum_targets",
    "CalibrationTargets",
    "Model",
    "Node",
    "Element",
    "UniaxialMaterial",
    "NDMaterial",
    "Support",
    "Load",
]
__version__ = "0.0.1"
