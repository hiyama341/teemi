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
with open('requirements.txt', 'r') as fh:
    for line in fh:
        requirements.append(line.strip())

extra_requirements={
    "dev": ["pytest==7.1.2", "pylint==2.13.9", "black==22.3.0", "pytest-cov==4.1.0"]
}

# Main setup scripts
setup(
    author="Lucas Levassor",
    author_email="lucaslevassor@gmail.com",
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    description="A Python package for constructing microbial strains",
    entry_points={
        "console_scripts": [
            "teemi=teemi.cli:main",
        ],
    },

    tests_require=extra_requirements["dev"],
    license="MIT license",
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/x-rst",  
    include_package_data=True,
    keywords="teemi",
    name="teemi",
    packages=find_packages(include=["teemi", "teemi.*"]),
    test_suite="tests",
    url="https://github.com/hiyama341/teemi",
    install_requires = requirements,
    ### Change version and put tag to release on PYPI
    # First time making a package, then comment this section and write: version= "0.0.1"
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
)
