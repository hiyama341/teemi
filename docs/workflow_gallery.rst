Workflow Gallery
================

`teemi` is most useful when it is approached through workflows rather than isolated functions.
This page highlights example notebooks that can be adapted to new hosts, targets, and design goals.

CAD Workflows
-------------

These notebooks live in ``teemi_cad_workflows/inspiration_notebooks/notebooks`` and are meant to be compact,
project-oriented examples.

CRISPR and Strain Engineering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``05_terminator_library_for_expression_tuning.ipynb``
  Outline for designing terminator libraries for expression tuning.
- ``06_crispr_multiplex_editing_planner.ipynb``
  Planning scaffold for multiplex CRISPR workflows.
- ``08_promoter_swap_for_overexpression_library_in_pichia.ipynb``
  Planning scaffold for promoter swap workflows in *Pichia*.
- ``09_crispr_cas9_deletion_example_in_aspergillus_oryzae.ipynb``
  Worked Cas9 deletion example with combined-genome off-target scanning and repair oligo design.
- ``10_crispr_cas12a_deletion_example_in_aspergillus_oryzae.ipynb``
  Worked Cas12a deletion example with the same Aspergillus workflow structure.

Library and Protein Engineering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``00_Designing_primers_for_kozak_library.ipynb``
  Example for designing primer sets for a Kozak library.
- ``01_HT_cazyme_primer_design.ipynb``
  High-throughput primer design workflow for CAZyme targets.
- ``02_primer_directed_mutagenesis.ipynb``
  Planning scaffold for primer-directed mutagenesis.
- ``03_investigate_plastic_degrading_enzymes.ipynb``
  Planning scaffold for enzyme discovery and prioritization.
- ``04_promoter_library_design_in_s_cerevisiae.ipynb``
  Planning scaffold for promoter library design in yeast.
- ``07_homolog_mining_to_codon_optimized_expression_panel.ipynb``
  Planning scaffold for homolog mining and expression panel generation.

Legacy Colab Notebooks
----------------------

The ``colab_notebooks`` directory contains larger, older, end-to-end case-study notebooks.
They are still useful, but the long-term direction of the repo is to expand the lighter workflow-oriented
examples in ``teemi_cad_workflows/inspiration_notebooks``.

Recommended Next Examples
-------------------------

The next useful workflow additions for the repo are:

- Cas12a multiplex deletion planning in *Aspergillus oryzae*
- genotyping planner for CRISPR edits
- Cas9 vs Cas12a guide comparison notebook
- signal peptide swapping workflow
- homolog mining to expression panel
