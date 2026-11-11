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
- Collapsed the three sgRNA finders into one scan. ``find_sgrna_hits_cas9``,
  ``find_sgrna_hits_cas12a`` and ``find_sgrna_hits_cas3`` were near-copies of
  one another -- the cas12a and cas3 bodies were 96% byte-identical, differing
  only in two constants and whether a ``print`` was commented out -- so a fix
  had to be made three times and in practice was not: a mutation audit found
  that 6 of 8 deliberate defects in the cas12a copy went undetected by the test
  suite, against 1 of 6 in cas3, which has golden-row tests. The three now
  delegate to one implementation driven by a ``_Nuclease`` record of everything
  that differs: PAM pattern and length, protospacer length, which side of the
  protospacer the PAM sits on, the strand sign, and the diagnostic messages.
  Their names, signatures and argument order are unchanged, since ``crispri``
  and the tests call them positionally.
- No change in what the finders return. Verified byte-identical over 43,663
  guide rows in 14 configurations: the *S. coelicolor* chromosome across 134
  loci for all three nucleases, two cas12a protospacer lengths, the ``filter``
  step, and synthetic circular, spliced and opposite-strand records. Preserved
  deliberately, because this change is not the place to alter them: the printed
  messages, including which nucleases stay silent; the cas12a and cas3
  ``sgrna_loc``, which the position formula misplaces (fixed separately,
  below); and the ``TypeError`` raised on a CDS whose parts lie on opposite
  strands, which Biopython gives ``strand=None``.
- A differential search over 2.95 million inputs -- truncated windows and
  each guard path, every compound-location arrangement, degenerate sequence
  content, parameter extremes, printed output, and every CDS of the
  *S. coelicolor* and *A. oryzae* genomes -- found three edge paths where the
  collapse had shifted behaviour, all since restored: an empty protospacer
  prints cas9's out-of-border message again, that nuclease having had no
  empty branch of its own; the window arithmetic stays inside the match loop,
  so a ``protospacer_len`` it cannot add still fails at the first match
  rather than at the call; and the strand multiplication stays ahead of the
  seed slice, which decides which ``TypeError`` a strandless CDS raises.
- One deliberate difference: ``find_sgrna_hits_cas9`` calls its ``revcomp``
  argument once per candidate window rather than twice, taking the guide and
  its PAM from one reverse complement instead of two. Any pure ``revcomp``,
  the module's own included, yields the same rows; only one that counts, logs
  or caches its calls can tell the difference.
- Dropped, in the course of that collapse: ``gene_counter`` and
  ``sgrna_counter``, incremented in three functions and never read; three
  commented-out ``print`` statements; and a docstring documenting a ``filepath``
  parameter the function did not take. The scan-flag branch was an ``if``/``elif``
  with no ``else``, leaving three names reachable-unbound, and is now an
  ``if``/``else``. The module goes from 941 lines to 699, and from 246 statements
  to 176.

- Fixed ``sgrna_loc`` for Cas12a and Cas3. The column is one past the
  protospacer's end, in coordinates of the CDS as extracted 5'->3'. It is how
  ``crispr_best.process_base_editing`` finds a guide, by slicing
  ``[sgrna_loc - len(sgrna) : sgrna_loc]``, so a value out by one base edits the
  wrong residue. Two errors with separate causes were both inert for Cas9 and so
  went unseen. The forward scan added a literal 3 to ``match.end()``, which for a
  5' PAM already covers the PAM -- Cas9's PAM pattern is a zero-width lookahead,
  so there the 3 supplies the PAM length rather than duplicating it -- and the
  reverse scan subtracted that same literal 3 where the nuclease's own
  ``pam_len`` belonged. Cas12a positions were therefore 3 nt out on the forward
  scan and 1 nt out on the reverse, and Cas3's 3 nt out on the forward scan only,
  its PAM being 3 nt. Base editing a Cas3 guide table raised ``RuntimeError`` in
  consequence -- the offset guard in ``process_base_editing`` catching the
  protospacer not being where the table said -- and now completes. A Cas12a
  table still fails there, but for an unrelated reason that this change does not
  touch: on SCO5087 it yields no editable guides at all, and
  ``process_base_editing`` raises ``ValueError`` on an empty frame rather than
  returning one.
- Checked by reading the slice back out of the CDS over the whole
  *S. coelicolor* chromosome: all 1,725,388 Cas9 guides already sat where they
  claimed, against none of the 8,124 Cas12a guides and none of the 75,555
  forward-scan Cas3 guides. Every guide of all three now does, and Cas9's eleven
  columns are unchanged across those 1,725,388 rows, so only the Cas12a and Cas3
  ``sgrna_loc`` moves. ``tests/test_crispr_validation.py`` asserts this per
  nuclease and requires both scan directions to occur, so that neither position
  formula can go unexercised.
- This is a point where teemi departs from StreptoCAD, which its guide finders
  were ported from. StreptoCAD computes the same two expressions, but every
  nuclease it supports has a 3 nt PAM -- its Cas12a is experimental and uses a
  3 nt ``TTN`` -- so the hardcoded 3 is never wrong on its reverse scan. teemi
  brought that half of the defect to light by adopting Cas12a's 4 nt ``TTTV``
  PAM. The forward-scan error is present in both, and misplaces Cas3 guides
  there as it did here.
- Pinned each nuclease's PAM and spacer length to the enzyme it was measured
  for. The three parameter sets were bare literals, so nothing in the module
  recorded which enzyme ``cas_type="cas3"`` meant, or where a 34 nt spacer and a
  ``TTC`` PAM had come from. Each entry now names its enzyme and cites a source
  that states both of its constants: ``cas9`` is SpCas9 (Jinek et al. 2012),
  ``cas12a`` is As/LbCas12a (Toth et al. 2018), and ``cas3`` is the minimal
  Type I-C CASCADE-Cas3 of pCRISPR-Cas3 (Whitford et al. 2025), which is this
  group's own system rather than a guess at one. All six constants were checked
  against the primary literature and none of them changed.
- Two things that surfaced in doing so, both worth recording. Cas3's subtype is
  I-C, and that does not follow from its PAM: Type I-B uses ``TTC`` as well, so
  the subtype comes from the system being modelled, while canonical Type I-E
  takes ``AAG`` instead. And Type I-C cuts bi-directionally from the target
  site, so the extent of a Cas3 deletion cannot be predicted to run one way.
- Left unpinned deliberately, both because they change what the finders return.
  ``off_target_seed`` defaults to 13 for all three nucleases, which is a Cas9
  figure, where Whitford et al. 2025 score Cas3 off-targets on the PAM plus the
  8 PAM-proximal bases. And ``SgRNAargs.protospacer_len`` defaults to 20 and is
  forwarded only to the cas12a finder, so Cas12a guides are 23 nt when that
  finder is called directly and 20 nt on the default ``extract_sgRNAs`` path.

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
