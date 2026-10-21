#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
"""A poisoned socket, shared by the ``test_combinatorial_*`` modules.

Not a ``conftest.py``: ``tests/test_PCR.py`` and ``tests/test_gibson_cloning.py``
live in the same directory and are allowed to reach NEB, so the guard is opted
into per module rather than imposed on the directory.

Module-scoped, because a module-scoped design fixture would otherwise be built
before a function-scoped guard could cover it; autouse fixtures are set up
ahead of other fixtures of the same scope, so this one wins.

``pydna.tm.requests`` does not exist -- ``tm_neb`` does ``import requests``
inside its own body -- so the verb on the ``requests`` module is what has to be
replaced, not an attribute of pydna.
"""

from __future__ import annotations

import pytest
import requests

import teemi.build.PCR


@pytest.fixture(autouse=True, scope="module")
def no_network():
    def forbidden(*args, **kwargs):
        raise AssertionError(
            "a combinatorial test reached the network; stub the reply instead"
        )

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(requests, "get", forbidden)
        patch.setattr(requests, "post", forbidden)
        patch.setattr(requests.Session, "request", forbidden)
        patch.setattr(teemi.build.PCR, "_post_neb_tm_api", forbidden)
        yield
