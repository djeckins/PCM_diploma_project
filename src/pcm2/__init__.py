"""pcm2 — per-frame prediction of ion-channel permeation readiness from MD trajectories.

Layers and boundaries: trajectory reading → features → models and evaluation → output.
The MD reading library (MDAnalysis) lives only in the reading layer and in the feature
modules that need PBC distances; the ML library appears in no feature module; the report
and figures read artifacts. The boundary is enforced by grep (tests/test_layer_boundaries.py).
"""

__version__ = "0.1.0"
