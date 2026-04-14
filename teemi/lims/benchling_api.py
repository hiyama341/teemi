"""Deprecated compatibility wrapper for legacy Benchling helpers."""

from warnings import warn

warn(
    "teemi.lims.benchling_api is legacy and will move out of the main API in a "
    "future release. Import from teemi.legacy.lims.benchling_api instead.",
    DeprecationWarning,
    stacklevel=2,
)

from teemi.legacy.lims.benchling_api import *  # noqa: F401,F403
