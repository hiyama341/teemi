#!/usr/bin/env python

"""Keep Versioneer's build hooks; package metadata lives in pyproject.toml."""
import sys
from pathlib import Path

# Isolated builds must be able to import the bundled Versioneer module.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import versioneer
from setuptools import setup

setup(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
