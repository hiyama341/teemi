"""Deprecated compatibility wrapper for legacy CSV database helpers."""

from warnings import warn

warn(
    "teemi.lims.csv_database is legacy and will move out of the main API in a "
    "future release. Import from teemi.legacy.lims.csv_database instead.",
    DeprecationWarning,
    stacklevel=2,
)

from teemi.legacy.lims.csv_database import *  # noqa: F401,F403
