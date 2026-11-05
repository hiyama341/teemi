History
-------

Unreleased
~~~~~~~~~~

- Repaired and modernized the documentation with a Sphinx + PyData Sphinx Theme setup.
- Added a cleaner front page, quickstart example, workflow gallery refresh, and improved Colab notebook links.
- Fixed the source distribution, which could not be built or installed because
  ``setup.py`` reads ``requirements.txt`` and ``MANIFEST.in`` did not ship it.
- Registered the ``dev``, ``docs`` and new ``learn`` extras. ``extras_require``
  was never passed to ``setup()``, so ``pip install teemi[dev]`` silently
  installed nothing; the 0.5.0 entry claiming this was fixed has been removed.
- Removed the ``teemi`` console script. It pointed at ``teemi.cli:main``, a
  module that has never existed, so the command only ever raised
  ``ModuleNotFoundError``.
- Declared ``requests``, ``seaborn`` and ``scipy``, which the package imports
  but never declared, and dropped the docs-only ``nbsphinx`` pin and the
  test-only ``pytest-cov`` pin from the runtime requirements.
- Corrected dependency floors that could not be satisfied together: ``pydna``
  requires ``numpy>1.26``, and ``pandas`` requires ``openpyxl>=3.1.0``.
- Fixed CRISPR base editing on minus-strand genes. ``process_base_editing``
  read ``sgrna_strand`` on its own, but it and ``gene_strand`` are both
  genome-relative, so a guide reads along the CDS only when the two agree.
  Every minus-strand edit was placed on the wrong residue.
- Fixed the Cas9 off-target seed register. ``find_off_target_hits`` took the
  seed as the bases ending at the PAM match, but the ``(?=GG)`` pattern matches
  at the first G of NGG, so the seed included the PAM's N and dropped the
  protospacer's last base, one base out of step with the slice the finders look
  the count up by. Cas9 and CRISPRi off-target counts therefore change, and the
  ``off_target_upper`` filter admits a different set of guides.
- Changed the default C-to-T editing window in ``identify_base_editing_sites``
  from 3-10 to 4-10, the window measured for CRISPR-cBEST. Position 3 has no
  reported conversion; pass ``editing_window_start=3`` for the old behaviour.
- ``process_base_editing`` now requires the ``gene_strand`` and ``sgrna``
  columns, and raises instead of printing a warning when a protospacer is not
  at its stated offset or an edited position is not a cytosine.
- ``find_off_target_hits`` accepts ``dcas9`` as an alias for ``cas9`` and names
  an unknown ``cas_type`` instead of raising ``KeyError``; ``SgRNAargs``
  rejects an ``off_target_seed`` outside ``1..protospacer_len``; and
  ``filter_guides`` builds its GC and off-target masks from the frame being
  filtered, fixing an ``IndexError`` on a multi-contig guide table.
- Added ``tests/test_crispr_validation.py``, which checks guide finding and
  base editing against CRISPy-web output and a hand curation of a minus-strand
  gene.
- Stopped the sgRNA finders reporting guides that the molecule does not carry.
  ``feature.extract()`` concatenates the parts of a compound location, so a
  window spanning two of them read a protospacer and PAM assembled across the
  gap -- across an intron on a spliced eukaryotic gene, or across the 10 bp
  overlap in ``SCO0022``, the one compound CDS in the *S. coelicolor*
  chromosome teemi validates against, which alone yielded 10 cas9 and 2 cas3
  guides occurring nowhere in the genome on either strand. Such guides also
  scored ``off_target_count`` 0, ranking them above real ones. The finders now
  skip windows that span a genomic discontinuity; parts that are *adjacent* on
  the DNA, such as a CDS crossing the origin of a circular record, are
  unaffected. Guide counts drop for spliced genes (2.95% to 0% phantoms on
  *A. oryzae* chromosome 1); one-part CDSs, which is every CDS in a bacterial
  genome bar the above, are untouched.

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
