History
-------

Unreleased
~~~~~~~~~~

- Repaired and modernized the documentation with a Sphinx + PyData Sphinx Theme setup.
- Added a cleaner front page, quickstart example, workflow gallery refresh, and improved Colab notebook links.

1.0.5 (2026-03-15)
~~~~~~~~~~~~~~~~~~

- Added the Aspergillus CRISPR workflow notebook.
- Introduced the ``teemi.legacy.build`` namespace to preserve older build helpers while simplifying the main build modules.
- Added tests around CRISPR guide logic and sequence fetching.

1.0.4 (2026-03-14)
~~~~~~~~~~~~~~~~~~

- Updated the PyPI publishing workflow.

1.0.3 (2026-03-14)
~~~~~~~~~~~~~~~~~~

- Added roadmap and release documentation to support project planning and maintenance.

1.0.2 (2025-11-14)
~~~~~~~~~~~~~~~~~~

- Fixed a small packaging/versioning issue in the release process.

1.0.1 (2025-02-15)
~~~~~~~~~~~~~~~~~~

- Fixed a PyPI packaging syntax issue in the README metadata used during release.

1.0.0 (2025-02-15)
~~~~~~~~~~~~~~~~~~

- Updated CI/test configuration around supported Python versions for the 1.0 release line.

0.5.1 (2024-10-01)
~~~~~~~~~~~~~~~~~~

- Updated packaging so project dependencies are included directly in the PyPI package.

0.5.0 (2023-11-08)
~~~~~~~~~~~~~~~~~~

- Expanded the README with more examples.
- Fixed installation of development extras such as ``pip install teemi[dev]``.

0.4.0 (2023-01-08)
~~~~~~~~~~~~~~~~~~

- Refactored the ``DesignAssembly`` class.
- Simplified redundant methods.
- Added support for multiple pads when constructing assemblies.

0.3.1 (2023-07-31)
~~~~~~~~~~~~~~~~~~

- Release attempt affected by a README bug.

0.3.0 (2023-06-22)
~~~~~~~~~~~~~~~~~~

- Added the ``gibson_cloning`` submodule for simple Gibson cloning workflows.

0.2.0 (2023-05-31)
~~~~~~~~~~~~~~~~~~

- Added the ``CRISPRsequencecutter`` and ``sequence_finder`` submodules.

0.1.0 (2023-02-01)
~~~~~~~~~~~~~~~~~~

- First release on PyPI.
