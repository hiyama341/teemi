"""Deprecated compatibility wrapper for legacy picklist helpers."""

from warnings import warn

warn(
    "teemi.build.containers_wells_picklists is legacy and will move out of the "
    "main build API in a future release. Import from "
    "teemi.legacy.build.containers_wells_picklists instead.",
    DeprecationWarning,
    stacklevel=2,
)

from teemi.legacy.build.containers_wells_picklists import *  # noqa: F401,F403
