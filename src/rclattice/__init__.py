"""rclattice — lattice modelling of reinforced concrete, analysed with OpenSees.

The Problem (problem.py), the FE Model (model.py), the mesher (mesh.py) and the builders
(builders.py) are backend-agnostic and do NOT import openseespy. All OpenSees calls live in
opensees.py.
"""

from .builders import (
    build_continuum,
    build_continuum_rc,
    build_lattice,
    build_lattice_rc,
    select_nodes,
)
from .calibration import CalibrationTargets, calibrate_lattice, continuum_targets
from .materials import (
    concrete_nd_elastic,
    concrete_nd_nonlinear,
    concrete_uniaxial_elastic,
    concrete_uniaxial_nonlinear,
    concrete_uniaxial_regularized,
    steel_uniaxial,
    steel_uniaxial_elastic,
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
    # problem / geometry / materials grades
    "Problem",
    "ConcreteGrade",
    "SteelGrade",
    "Rebar",
    "RectangleDomain",
    "CompoundRectangles",
    "portal_frame",
    "EdgeSupport",
    "EdgeLoad",
    "BoxSupport",
    "BoxLoad",
    # builders
    "build_lattice",
    "build_lattice_rc",
    "build_continuum",
    "build_continuum_rc",
    "select_nodes",
    # calibration
    "calibrate_lattice",
    "continuum_targets",
    "CalibrationTargets",
    # material mappings (grade -> OpenSees material)
    "concrete_uniaxial_elastic",
    "concrete_uniaxial_nonlinear",
    "concrete_uniaxial_regularized",
    "concrete_nd_elastic",
    "concrete_nd_nonlinear",
    "steel_uniaxial",
    "steel_uniaxial_elastic",
    # FE model objects
    "Model",
    "Node",
    "Element",
    "UniaxialMaterial",
    "NDMaterial",
    "Support",
    "Load",
]
__version__ = "0.0.1"
