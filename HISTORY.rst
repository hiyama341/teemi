History
-------

Unreleased
~~~~~~~~~~

1.2.0 (2026-11-19)
~~~~~~~~~~~~~~~~~~

.. warning::

   **Upgrading changes CRISPR results, and nothing raises to say so.**
   Guide tables, off-target counts and base-edit calls produced with 1.1.0
   or earlier do not match what this release produces from the same inputs.
   Regenerate anything you mean to order, publish or build on. The changes
   that move results silently are: every base edit on a minus-strand gene,
   which was placed on the wrong residue; the default C-to-T editing window,
   from 3-10 to 4-10; Cas9 and CRISPRi ``off_target_count`` wherever a guide
   has off-targets at all, the seed having been read one base out of step;
   Cas12a spacers, from 20 nt to 23 nt on a 4 nt ``TTTV`` PAM; Cas12a and
   Cas3 ``sgrna_loc``, by up to 3 nt; ``editing_context``, and so which rows
   survive ``editing_context=True``, for a window starting at position 1,
   which used to read the guide's LAST base as the context for its first; and
   guide tables for genes with a compound
   location, which shrink, because a guide spanning one of its gaps is no
   longer reported. Base editing a Cas12a or Cas3 guide table raises
   instead, where it used to return invented amino-acid changes.

- Fixed CRISPR base editing on minus-strand genes. ``process_base_editing``
  read ``sgrna_strand`` on its own, but it and ``gene_strand`` are both
  genome-relative, so a guide reads along the CDS only when the two agree.
  Every minus-strand edit was placed on the wrong residue.
- ``process_base_editing`` now refuses a guide table that has not been through
  ``filter_sgrnas_for_base_editing``, instead of failing on it. The column
  guard's message named the finders and ``identify_base_editing_sites`` but not
  that filter, which is not optional: ``identify_base_editing_sites`` leaves
  ``editable_cytosines`` empty for a guide with no cytosine in the window --
  319 of the 7,651 cas9 guides over the first 50 single-part CDSs of the
  *S. coelicolor* validation chromosome, SCO0001 to SCO0058, which is 4.2% of
  them and reaches 45 of the 50 genes -- and following the message literally
  reached ``int('')``
  and raised pandas' ``invalid literal for int() with base 10: ''``, naming
  neither the column nor the step that had been left out. The message now names
  all three steps, and an empty ``editable_cytosines`` is refused by name and
  count.
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
  ``if``/``else``. Measured against the 1.1.0 tag, the collapse takes
  ``crispr_cas.py`` from 216 ``ast`` statements to 190, while the file grows
  from 824 lines to 856 -- the sourcing work further down this list is longer
  than the code it documents, and the two numbers move opposite ways for that
  reason rather than by accident.
- ``process_base_editing`` and ``identify_base_editing_sites`` now refuse a
  guide table from any nuclease but Cas9. CRISPR-BEST is an SpCas9 system:
  its editors are rAPOBEC1 (cBEST) and ecTadA (aBEST) fused to a
  *Streptomyces* codon-optimised SpCas9 D10A nickase, and the default 4-10
  editing window is the one measured for cBEST on that 20 nt geometry (Tong
  et al. 2019, doi:10.1073/pnas.1913493116). A guide table carries no
  nuclease label -- the finders' frame is eleven columns and ``cas_type`` is
  not one of them -- so the decision is made on the ``pam`` column, which
  must end in ``GG``. Every SpCas9 PAM does, and neither other nuclease can
  produce one that does: Cas12a's pattern is ``TTT[ACG]`` and Cas3's is
  ``TTC``. ``process_base_editing`` therefore now also requires the ``pam``
  column, and ``editing_context`` when the flag of that name is left True.
  The check is strict there, where an amino-acid change gets asserted, and
  advisory in ``identify_base_editing_sites``, which only reports where the
  cytosines in a spacer are and so checks the PAMs when present without
  demanding them. Two limits are deliberate: a ``pam`` cell is an ordinary
  string, so a caller who overwrites one is believed, and the PAM fixes the
  nuclease but not the geometry -- a table found with ``cas_type="cas9"``
  and ``protospacer_len`` overridden carries real NGG PAMs and still has the
  4-10 window read off a spacer that window was never measured on. An empty
  guide table carries no PAM to judge and is allowed through; it now comes
  back empty rather than raising the bare pandas ``ValueError`` described
  below.
- Fixed ``sgrna_loc`` for Cas12a and Cas3. The column is one past the
  protospacer's end, in coordinates of the CDS as extracted 5'->3'. It is
  how ``crispr_best.process_base_editing`` finds a guide, by slicing
  ``[sgrna_loc - len(sgrna) : sgrna_loc]``, so a value out by one base edits
  the wrong residue. Two errors with separate causes were both inert for
  Cas9 and so went unseen. The forward scan added a literal 3 to
  ``match.end()``, which for a 5' PAM already covers the PAM -- Cas9's PAM
  pattern is a zero-width lookahead, so there the 3 supplies the PAM length
  rather than duplicating it -- and the reverse scan subtracted that same
  literal 3 where the nuclease's own ``pam_len`` belonged. Cas12a positions
  were therefore 3 nt out on the forward scan and 1 nt out on the reverse,
  and Cas3's 3 nt out on the forward scan only, its PAM being 3 nt. Base
  editing a Cas3 guide table raised ``RuntimeError`` in consequence -- the
  offset guard in ``process_base_editing`` catching the protospacer not
  being where the table said -- and ran to completion once the positions
  were right. That was a regression this change introduced, recorded here
  because the ``RuntimeError`` had accidentally been the only thing stopping
  a scope error: with the positions correct a Cas3 table returned base edits
  instead of refusing. On SCO5087 all 31 of its guides were reported as
  carrying an editable cytosine, and ``process_base_editing`` returned five
  rows of invented amino-acid changes at its default arguments -- 28 with
  ``editing_context=False`` -- every row carrying a 34 nt spacer and a
  ``TTC`` PAM, the C-to-T window having been read off positions 4 to 10 of a
  34-mer. All of it was fiction, and the nuclease guard above now refuses
  it. A Cas12a table failed instead of lying, for a reason this change does
  not touch and which is not about base editing at all: SCO5087 yields no
  Cas12a **guides**. Its 1,404 nt CDS is 72.1% GC and carries not one
  ``TTT`` on either strand, so there is no TTTV site to find, and none would
  appear if the PAM were widened to TTTN. What reached
  ``process_base_editing`` was an empty frame, and the pandas assignment
  raised a bare ``ValueError`` reading "Cannot set a DataFrame with multiple
  columns to the single column mutated_sequence" -- not the module's own
  column guard above it, so the message named nothing a caller could act
  on. The entry above returns an empty frame instead.
- Checked by reading the slice back out of the CDS over the whole
  *S. coelicolor* chromosome: all 1,725,388 Cas9 guides already sat where they
  claimed, against none of the 8,124 Cas12a guides and none of the 75,555
  forward-scan Cas3 guides. Every guide of all three now does, and Cas9's eleven
  columns are unchanged across those 1,725,388 rows, so only the Cas12a and Cas3
  ``sgrna_loc`` moves. ``tests/test_crispr_validation.py`` asserts this per
  nuclease and requires both scan directions to occur, so that neither position
  formula can go unexercised.
- The same two expressions are in StreptoCAD, which teemi's guide finders were
  ported from, and the forward-scan error was live there too. StreptoCAD has
  since dropped the ``+ 3`` from its Cas3 finder, so the two now agree on every
  Cas3 ``sgrna_loc`` -- on SCO5087, 31 guides each, same spacers and PAMs,
  positions equal on all 31. Its reverse scan still subtracts a literal 3, and
  still cannot be wrong there, because every nuclease it supports has a 3 nt
  PAM. That is the half teemi had to fix: teemi keeps the 4 nt ``TTTV`` Cas12a
  that both projects inherited from the same pre-split finder and that
  StreptoCAD has since deleted.
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
- Fixed Cas12a guides coming back 20 nt long from ``extract_sgRNAs``. Each
  parameter set carries the spacer length its enzyme was measured for, but
  ``SgRNAargs.protospacer_len`` defaulted to 20 and was forwarded to the cas12a
  finder alone, so the default path overrode the As/LbCas12a 23 with a Cas9
  number while cas9 and cas3 kept their own. ``protospacer_len`` now defaults to
  ``None``, meaning whatever the chosen ``cas_type`` was measured for, and an
  explicit value applies to whichever nuclease is in use rather than to cas12a
  only. Over 40 *S. coelicolor* CDSs this moves Cas12a from 89 guides of 20 nt
  to 88 of 23 nt, which is what calling the finder directly already returned.
  Notebook 10 had been passing ``protospacer_len=23`` by hand to get it.
- Fixed the ``off_target_seed`` ceiling refusing legal values. It was checked
  against ``protospacer_len``, which was 20 whatever the nuclease, so
  ``cas_type="cas3", off_target_seed=25`` was rejected although Cas3's
  protospacer is 34 nt long and the seed is sliced out of it. The ceiling is now
  the protospacer length the ``cas_type`` resolves to, and the message names the
  nuclease that set it.
- Replaced the three substring guards in ``extract_sgRNAs`` with one lookup in
  the parameter table. ``if "cas9" in args.cas_type`` is what had made ``dcas9``
  work, it matched a nuclease named anywhere in the string, and a ``cas_type``
  none of the three recognised fell through all of them and left ``sgrna_df``
  unbound. An unknown ``cas_type`` is now rejected when ``SgRNAargs`` is
  constructed, since resolving the protospacer length is what validates it.
- Moved the off-target scan's PAM patterns into that same table as
  ``off_target_pam_pattern``. They were a second dictionary holding Cas9's PAM
  in the opposite orientation -- ``(?=GG)`` against ``(?=CC)`` -- with nothing
  recording that the two scans read their seed off different strands and so need
  the PAM in different senses. Both scans read both strands, so the seeds
  collected were always the same; the counters are unchanged, verified identical
  across four ``cas_type`` values and four seed lengths over both strands of the
  chromosome and a set of degenerate short inputs.
- Settled ``off_target_seed`` by leaving all three nucleases on 13. It is not a
  stray Cas9 constant: it is CRISPy-web's own window, and Blin et al. 2020 sort
  candidates by off-target hits "with 0, 1, or 2 mismatches in the 13 bp seed
  region right upstream of the PAM". CRISPy-web holds that 13 fixed and varies
  only the PAM between nucleases, so one seed length for all three is upstream's
  design rather than something overlooked here. Whitford et al. 2025's Cas3
  figure is real and ``off_target_seed=8`` reproduces it exactly, but as a
  default it costs more than it buys: over the first 500 CDSs of NC_003888.3,
  Cas3 at 13 keeps all 8712 candidate guides under the default
  ``off_target_upper=10`` with a median count of 0, while at 8 it keeps 54.2%,
  the median count is 9, the 90th percentile is 32 and three genes are left with
  no guide at all. The seed and that cutoff are a coupled pair, and Whitford's
  script applies no numeric cutoff, so moving one of the two alone does not
  reproduce their method. The parameter stays overridable and its provenance is
  now recorded beside the table.
- Recorded the gap that matters more than the seed length. ``off_target_count``
  counts exact matches, where CRISPy-web counts near-misses in the same window.
  This repository's own fixture carries both: guide CY00000196 of SCO5087 has 0
  off-targets at 0 mismatches, 30 at one and 591 at two, and teemi reproduces
  the first of those three columns only. So ``off_target_upper`` rejects almost
  nothing at the shipped defaults -- 0.2% of Cas9 candidates over those 500
  CDSs -- and a guide reported as unique may have hundreds of near-misses.
  Nothing in the test suite guards the default either, since every call site
  passes ``off_target_seed`` explicitly.
- Corrected the ``off_target_upper`` docstring, which called it a "maximum
  number of seed nucleotides" where ``filter_guides`` compares it against a
  guide's count of off-target sites.
- Added the first external reference for Cas12a. The nuclease had none: the
  CRISPy-web export validates Cas9, and SCO5087, the gene it covers, carries no
  TTTV site on either strand of this GC-rich chromosome. Two *S. cerevisiae*
  genes designed with LbCas12a in GenScript's tool fill the gap -- CAN1
  (YEL063C, minus strand) and LYP1 (YNL268W, plus strand), 55 guides between
  them -- with each gene's own slice of its chromosome taken from RefSeq so the
  comparison runs offline. teemi finds every one of those guides, with the same
  PAM, at the same chromosome coordinate. This checks what no test written
  against teemi's own parameter table can: GenScript encodes the TTTV PAM and
  the 23 nt spacer itself, so agreeing with it tests the two constants rather
  than only the code that reads them.
- Recorded that excluding TTTT from the Cas12a PAM is load-bearing twice over.
  It is the measured As/LbCas12a consensus, and it is also what keeps the
  non-overlapping ``re.finditer`` scan complete: no two TTTV sites can start
  within 4 nt of each other, so none is skipped. TTTN has no such property --
  over this repository's chromosome it finds 7,996 sites where an overlapping
  scan finds 9,338, and over the two yeast genes it loses 14 of the guides
  GenScript reports.
- Added the first external reference for Cas3. No design service offers
  Type I-C: GenScript has no Type I system, and CRISPy-web 3.0 added CRISPRi
  and TnpB but is still Cas9-only among the cutters. The authors of
  pCRISPR-Cas3 hit the same wall -- "the majority of existing computational
  tools do not support type I/type I-C CRISPR systems" -- and published the
  notebook they picked their own spacers with, `NBChub/Type1_Protospacers
  <https://github.com/NBChub/Type1_Protospacers>`_. That makes it the
  reference implementation of the system modelled here, and the same source
  ``TTC`` and ``34`` are taken from. Over SCO5087 (plus strand) and SCO7048
  (minus strand) it reports 31 spacers each and teemi reports exactly the same
  sets -- same sequence, same chromosome coordinate and same strand, guide for
  guide, with equality rather than containment because the notebook ranks but
  does not truncate. Over the first 400 single-part CDSs of NC_003888.3 all
  6,924 guides agree on both. It is the same laboratory as teemi and as
  CRISPy-web, so it does not make the biology independent; what it pins
  independently is the mechanics, which side of the protospacer the PAM sits
  on, how far the spacer runs, how the reverse strand is reached and where the
  result lands -- the four things that were wrong for cas12a and cas3 before.
- Confirmed the ``off_target_seed=8`` claim against that notebook rather than
  against its paper's prose. teemi's ``off_target_count`` at 8 is the
  notebook's seed count minus one for every guide of both genes, up to counts
  of 53. The offset is the on-target site, which teemi subtracts and the
  notebook leaves in, as its README says: "a match count of 1 indicates that
  no off-target sites were found". A count that reached only one strand, or
  anchored the seed a base out, would diverge by more than a constant.
- Added the spacers pCRISPR-Cas3 was actually cloned with, which is the only
  CRISPR reference here made of guides known to cut rather than guides a scan
  proposed. The paper's data record (DTU Data, doi:10.11583/DTU.22786157)
  lists ten protospacers targeting this chromosome; seven lie inside a
  single-part CDS and so are reachable by a locus tag, and teemi reports every
  one of them at the coordinate and strand the record gives, across six genes
  from 92 kb to 7.4 Mb, four on the minus strand and three on the plus. The
  other three are intergenic, which ``extract_sgRNAs`` cannot reach because it
  takes locus tags and nothing else -- Cas3 deletes processively, so the paper
  also aims it between genes. This cannot validate ``TTC`` or ``34``, since
  those spacers were chosen by searching for that PAM at that length, and the
  test says so; what it validates is the scanner, over far more of the
  chromosome than two genes reach.
- Recorded a source for the Cas3 PAM from outside teemi's own authors, which
  the parameter table had been missing. Leon et al. 2021 (NAR 49:2114,
  doi:10.1093/nar/gkab006) determined it in the system's native host and by
  the laboratory that built this compact I-C system rather than the one that
  ported it: "Alignments of 3163 protospacers with upstream and downstream
  regions revealed the consensus PAM ... to be 5'-TTC-3'". The same paper puts
  34 nt inside a measured range instead of at an optimum -- "Because Type I-C
  spacer length ranges from 32-37 nt, contrary to Type I-F spacers
  consistently measuring 32 nt" -- so 34 is the pCRISPR-Cas3 cloning
  convention chosen within a functional window. The notebook credits Csorgo
  et al. 2020 for its defaults; that paper introduced the system, and the PAM
  determination for it is Leon.
- Noted one published number that does not reproduce, so that it is not
  mistaken for an oracle later. Whitford et al. report "158 341 5'-TTC-3' PAM
  sequences" in *S. coelicolor*, "averaging one PAM for every 54.7 bp", where
  scanning NC_003888.3 on both strands finds 162,382 and collapsing identical
  spacers gives 159,792. Neither is their figure and the 2.5% gap is
  undocumented.
- Re-executed the two *Aspergillus oryzae* CRISPR notebooks, 09 (Cas9) and 10
  (Cas12a), whose committed output predated all of the above. Both ran against
  the eight local RIB40 chromosomes, 37,848,727 bp and 60,033 features. Two of
  the changes here move them: the splice fix replaced one of
  ``AO090009000002``'s two chosen Cas9 guides, ``GTATTAGATGGGTATTTGGA`` at
  ``sgrna_loc`` 3596 giving way to ``CGATTATCATTCAACCAGAT`` at 3504, and the
  Cas12a 20 -> 23 nt correction reaches notebook 10's whole guide table.
- Notebook 10's committed output had never been produced by notebook 10. All
  twelve of its code cells' outputs were byte-identical to notebook 09's, the
  two whose source differs included, so it reported 20 nt spacers under the
  label ``sgRNA spacer`` although its own code passes ``protospacer_len=23``
  and sets ``item_type='Cas12a crRNA spacer'``; both notebooks also carried the
  same execution counts, 17 to 28, which is the tell. It now carries its own
  output: 23 nt spacers on ``TTTA``, ``TTTC`` and ``TTTG`` PAMs. This was a
  pre-existing defect rather than one the fixes above caused.
- Fixed an ``in_tail`` part taking the whole library down with it. Such a
  part has two independent limits: ``in_tail_max_bp``, what a primer tail
  will physically carry, which is pydna's ``maxlink``; and
  ``homology_length``, what the overlap has to *spend*, because an
  ``in_tail`` part sits inside the homology rather than being added to it.
  Only the first was checked, and the two defaults sit close together -- 40
  against a ``homology_length`` of 35 -- so a part between 36 and 40 bp was
  admitted by ``plan_layout`` and then refused by ``design_junction``,
  which raised ``JunctionError`` and returned no variants, no PCRs and
  nothing designed. ``plan_layout`` now demotes such a part to its own
  fragment, as it already did for the physical cap: one more PCR rather
  than no library. The two ``JunctionError`` cases are reported apart now,
  because their remedies differ -- a part longer than the whole homology is
  a number in the spec, whereas two parts too short to supply the budget is
  the choice of parts, and printing the part lengths against a negative
  budget pointed at the parts for a failure that was neither of theirs.
- A refused plan no longer returns its stop alone. ``_stopped`` built exactly
  one ``Diagnostic``, so a ``no_footprint`` or ``impossible_junction`` result
  had a one-element ``diagnostics`` and ``diagnostics[0]`` was the stop; the
  warnings already settled before it were dropped, which reported a part that
  could not be used and said nothing about the ones that had already been
  demoted. Both stops now carry what was established first, so the critical
  entry is LAST rather than only. A consumer that keyed on "a refused design
  has one diagnostic, and it is the stop" has to read the last one instead.
- ``plan_layout``'s ``in_tail`` fix above emits a second diagnostic code,
  ``in_tail_over_homology``, beside the existing ``in_tail_too_long``. The
  string is named here because ``Diagnostic.code`` crosses the boundary through
  ``DesignResult.to_dict``, so anything filtering or grouping on codes meets a
  value it has not seen.
- ``filter_guides`` refuses a pattern filter on a column the finders do not
  produce -- ``ValueError`` naming the column, where it used to reach
  ``frame[column]`` and raise ``KeyError``. ``downstream_remove`` and
  ``upstream_remove`` have always been such columns: ``_HIT_COLUMNS`` is
  eleven and neither is among them.
- ``SgRNAargs.protospacer_len`` is resolved on construction, not just
  defaulted. The entry above says it "defaults to ``None``", which describes
  the parameter; the ATTRIBUTE is the nuclease's own length by the time
  anything reads it, so ``SgRNAargs(..., cas_type="cas3").protospacer_len`` is
  34 and cas12a's is 23, where 1.1.0 gave 20 for both. Anything that read it
  to slice, label or report reads a different number.
- The ``Development Status`` classifier goes from ``2 - Pre-Alpha`` to
  ``5 - Production/Stable``. It is metadata rather than behaviour, but it is
  what PyPI shows, and this is the release that moves it.
- ``teemi.design.combinatorial.flow`` imported ``Annotated`` from
  ``typing_extensions``, which teemi does not declare and had only because
  pydantic depends on it. It comes from ``typing`` now, where it has been
  since Python 3.9, against a floor of 3.11, so nothing need be declared.
- Added ``teemi.design.combinatorial.idt``: the primer order in the two shapes
  IDT takes it in. ``tube_frame`` is its bulk tube form, ``plate_frames`` one
  frame per 96-well plate with the wells filled **down the columns** -- an
  8-channel pipette spans exactly one column, so a part-full plate's oligos
  sit in whole columns instead of being smeared across every row -- and
  ``primers_fasta`` the same oligos for whatever reads no CSV. Both frames
  read ``result.primers`` in order, so the plate is a re-shaping of the order
  sheet and the two cannot disagree about what is in the order. The module is
  pure; ``export.py`` remains the only file-writing edge of the package.
- ``oligo_name`` puts the primer's label before its id, which is what makes
  the length cap safe: a label is unique within one plan by construction, so
  two long ids that truncate alike still give two distinct names. It does not
  collapse runs of underscores: ``ID_SEPARATOR`` is ``__``, so collapsing
  doubles WOULD make ``A__B`` and ``A_B`` the same oligo. Conditional because
  no release had it -- a sanitiser ported from elsewhere collapsed them, and
  the module's own parametrised test caught it before any of this shipped.
- The bundle gained six files and a provenance stamp:
  ``02_primers/oligo_order_idt.csv``, ``02_primers/oligo_order_idt_plate.xlsx``
  and ``02_primers/primers.fasta``; ``03_pcrs/worklist_by_ta.csv``;
  ``05_constructs/constructs.csv``, which is ``variant_overview_dataframe`` --
  built on every plan and, until now, never written, so the archive carried a
  map per construct and no table of what was being built; and
  ``05_constructs/00_all_constructs.gb``. The xlsx is the only new binary and
  needs no new dependency, ``openpyxl`` already being one.
- ``write_bundle``, ``write_bundle_zip`` and ``build_summary`` take ``run=``, a
  mapping the caller supplies, written out as ``09_run.json`` beside the seven
  things teemi can see for itself: ``teemi``, ``pydna``, ``primer3``,
  ``pandas`` and ``Bio``, plus ``python`` and ``platform``. Same argument
  ``part_info`` already makes: provenance this package cannot know. The keyword
  is defaulted and the file conditional, so every EXISTING caller is
  unaffected, and ``00_LAYOUT.txt``, which is built from what was actually
  written, cannot promise a file that is not there.
- Two things a NEW caller of ``run=`` has to know, since the keyword is not as
  opaque as it looks. Those seven names are reserved: ``run_json`` raises
  ``ValueError`` on a key that collides with one, because a bundle naming a
  teemi version its own code did not produce cannot be checked, and the clash
  is likelier a mistake than an intention. And ``build_summary`` reads one key
  by name -- ``built_at``, which it prints under the rule at the top of
  ``01_SUMMARY.txt`` -- so a caller spelling it ``timestamp`` or ``date`` gets
  no line and no complaint. Nothing here reads the clock: the timestamp is the
  caller's, because it is the one fact that has to agree with their own
  filenames and log lines.
- ``bundle_layout``'s description column starts at 40 characters instead of
  30, because the plate workbook's path is 37 and the old gutter wrapped it
  into the text. Every line of ``00_LAYOUT.txt`` therefore moves, which is
  worth saying out loud because that file is the one in the bundle meant to be
  read rather than parsed. It also gained a keyword-only ``omit``, which drops
  individual lines from an entry that was written -- additive, and public,
  ``bundle_layout`` being in the package surface.
- ``build_summary`` states how many plates the oligos fill and, under 24 of
  them, that the order is below IDT's floor for a 96-well plate. The observed
  count against the limit, the way a ``Violation`` reads, rather than advice.
- Added ``pcr_worklist_dataframe``, the reactions grouped by the temperature
  they anneal at. It groups on the **whole degree with ties going down**, which
  is the rule ``pcrs.csv`` already prints Ta by, now extracted as
  ``report.whole_ties_down`` so both readers share one owner: grouping on the
  raw float would turn 61.2 and 61.4 into two blocks of one reaction each,
  since an offline plan gets its Tm from pydna in tenths, and could label a
  block a degree away from the Ta printed for that same reaction. It reports no
  thermocycler count, because a gradient block covers a range.
- Fixed ``constructs_genbank`` holding one record. pydna ends a formatted
  record at ``//`` with no trailing newline, so joining them directly gave
  ``//LOCUS`` on one line and a file that holds a single record as far as any
  parser is concerned. Counting ``LOCUS`` lines is what missed it; the test
  parses the file back.
- Fixed LOCUS name collisions in the same file. ``simulate`` names a product
  from ``variant.id[:16]``, which is not distinct across variants -- 8 designs
  gave 4 distinct names, 12 gave 8 -- so the combined file re-stamps each
  record from its label and id on a copy, leaving the caller's records alone.
- Added ``construct_maps(folder)``, and this is the one thing in the bundle
  work that an existing caller has to know about: ``glob("*.gb")`` over
  ``05_constructs`` now matches the combined file too, which is a ``.gb`` that
  is not a construct map. Anything counting those files counts one too many,
  and anything handing each to ``SeqIO.read`` gets a multi-record file --
  which is exactly what four of teemi's own tests did. Use
  ``construct_maps(folder)`` rather than the glob.

1.1.0 (2026-10-30)
~~~~~~~~~~~~~~~~~~

- Added the ``teemi.design.combinatorial`` subpackage, the JSON-spec driven
  combinatorial design backend, in `#47
  <https://github.com/hiyama341/teemi/pull/47>`_.
- Extended support to Python 3.12, 3.13 and 3.14.
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
