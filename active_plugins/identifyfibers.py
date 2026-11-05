"""CellProfiler discovery entry point for the IdentifyFibers plugin.

Keep this file and the adjacent ``identifyfibers_engine`` package directory together
inside the directory selected in CellProfiler Preferences as the plugins
directory.
"""

from identifyfibers_engine.cellprofiler_module import IdentifyFibers

__all__ = ["IdentifyFibers"]

