"""FOLDTRACE Guided: an undergraduate research workflow that paces a structure-first
active-site investigation into observe -> predict -> compute -> compare -> decide steps.

It orchestrates the released FOLDTRACE v0.2.0 functions (foldtrace.search, foldtrace.mapping)
without changing their scientific behaviour, thresholds, or output definitions. FOLDTRACE
tracks what has happened to *prespecified* chemistry-defining positions from an experimentally
characterized reference; it does not discover catalytic residues de novo.
"""
from .project import (
    GuidedProject, GuidedError, CheckpointError, PredictionLockedError,
    Stage, STAGES, STAGE_IDS, PREDICTABLE_STATES, TM_GATE, OFFSET_A,
)
__all__ = ["GuidedProject", "GuidedError", "CheckpointError", "PredictionLockedError",
           "Stage", "STAGES", "STAGE_IDS", "PREDICTABLE_STATES", "TM_GATE", "OFFSET_A"]
