#!/usr/bin/env python

"""The setup script."""
import versioneer
from setuptools import setup, find_packages

# Readme & History
with open("README.rst") as readme_file:
    readme = readme_file.read()
with open("HISTORY.rst") as history_file:
    history = history_file.read()

# Requirements
requirements = []
with open("requirements.txt", "r") as fh:
    for line in fh:
        line = line.strip()
        if line and not line.startswith("#"):
            requirements.append(line)

extra_requirements = {
    "dev": [
        "pytest>=9.0.3",  # CVE-2025-71176: tmpdir handling
        "pylint>=2.13.9",
        "black>=26.3.1",  # CVE-2024-21503, CVE-2026-32274, CVE-2026-31900
        "pytest-cov>=7.1.0",
    ],
    "docs": ["sphinx", "pydata-sphinx-theme", "myst-parser"],
    # teemi.learn.auto_ml drives H2O, which pulls in a JVM-backed server;
    # it is only needed for the LEARN notebooks.
    "learn": ["h2o>=3.38"],
}

setup(
    author="Lucas Levassor",
    author_email="lucaslevassor@gmail.com",
    python_requires=">=3.11",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    description="A Python package for constructing microbial strains",
    license="MIT license",
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/x-rst",  
    include_package_data=True,
    keywords="teemi",
    name="teemi",
    packages=find_packages(include=["teemi", "teemi.*"]),
    url="https://github.com/hiyama341/teemi",
    project_urls={
        "Documentation": "https://teemi.readthedocs.io",
        "Source": "https://github.com/hiyama341/teemi",
        "Paper": "https://doi.org/10.1371/journal.pcbi.1011929",
    },
    install_requires=requirements,
    extras_require=extra_requirements,
    ### Change version and put tag to release on PYPI
    # First time making a package, then comment this section and write: version= "0.0.1"
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
