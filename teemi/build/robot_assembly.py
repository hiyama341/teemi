"""Deprecated compatibility wrapper for legacy robot assembly helpers."""

from warnings import warn

warn(
    "teemi.build.robot_assembly is legacy and will move out of the main build "
    "API in a future release. Import from teemi.legacy.build.robot_assembly "
    "instead.",
    DeprecationWarning,
    stacklevel=2,
)

from teemi.legacy.build.robot_assembly import *  # noqa: F401,F403
