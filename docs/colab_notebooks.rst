Colab Notebooks
===============

These notebooks capture the larger, end-to-end strictosidine case study workflows used with
``teemi``. They are the best starting point if you want to run the published DBTL examples in
Google Colab with as little setup as possible.

Each notebook links to both Google Colab and the GitHub source file.

Strictosidine Case: First DBTL Cycle
------------------------------------

Design
~~~~~~

.. list-table::
   :widths: 7 53 40
   :header-rows: 1

   * - ID
     - Focus
     - Links
   * - 00
     - Automatically fetch homologs from NCBI in a standardizable and repeatable way.
     - `Open in Colab <https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/00_1_DESIGN_Homologs.ipynb>`__ | `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/00_1_DESIGN_Homologs.ipynb>`__
   * - 01
     - Select promoters from RNA-seq data and fetch them from online databases with quality checks.
     - `Open in Colab <https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/01_1_DESIGN_Promoters.ipynb>`__ | `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/01_1_DESIGN_Promoters.ipynb>`__
   * - 02
     - Generate a combinatorial library with ``DesignAssembly`` and prepare it for downstream build planning.
     - `Open in Colab <https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/02_1_DESIGN_Combinatorial_library.ipynb>`__ | `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/02_1_DESIGN_Combinatorial_library.ipynb>`__

Build
~~~~~

.. list-table::
   :widths: 7 53 40
   :header-rows: 1

   * - ID
     - Focus
     - Links
   * - 03
     - Assemble a CRISPR plasmid with USER cloning.
     - `Open in Colab <https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/03_1_BUILD_gRNA_plasmid.ipynb>`__ | `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/03_1_BUILD_gRNA_plasmid.ipynb>`__
   * - 04
     - Construct the background strain by knocking out G8H and CPR.
     - `Open in Colab <https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/04_1_BUILD_Background_strain.ipynb>`__ | `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/04_1_BUILD_Background_strain.ipynb>`__
   * - 05
     - Generate the first combinatorial library covering 1280 possible combinations.
     - `Open in Colab <https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/05_1_BUILD_Combinatorial_library.ipynb>`__ | `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/05_1_BUILD_Combinatorial_library.ipynb>`__

Test
~~~~

.. list-table::
   :widths: 7 53 40
   :header-rows: 1

   * - ID
     - Focus
     - Links
   * - 06
     - Process LC-MS data and genotype the generated strains.
     - `Open in Colab <https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/06_1_TEST_Library_characterisation.ipynb>`__ | `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/06_1_TEST_Library_characterisation.ipynb>`__

Learn
~~~~~

.. list-table::
   :widths: 7 53 40
   :header-rows: 1

   * - ID
     - Focus
     - Links
   * - 07
     - Use AutoML to predict the best combinations for a targeted second round of library construction.
     - `Open in Colab <https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/07_1_LEARN_Modelling_and_predictions.ipynb>`__ | `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/07_1_LEARN_Modelling_and_predictions.ipynb>`__

Strictosidine Case: Second DBTL Cycle
-------------------------------------

Design
~~~~~~

.. list-table::
   :widths: 7 53 40
   :header-rows: 1

   * - ID
     - Focus
     - Links
   * - 08
     - Translate machine-learning results into a targeted library of strains.
     - `Open in Colab <https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/08_2_DESIGN_Model_recommended_combinatiorial_library.ipynb>`__ | `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/08_2_DESIGN_Model_recommended_combinatiorial_library.ipynb>`__

Build
~~~~~

.. list-table::
   :widths: 7 53 40
   :header-rows: 1

   * - ID
     - Focus
     - Links
   * - 09
     - Construct the targeted second-round library of strains.
     - `Open in Colab <https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/09_2_BUILD_Combinatorial_library.ipynb>`__ | `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/09_2_BUILD_Combinatorial_library.ipynb>`__

Test
~~~~

.. list-table::
   :widths: 7 53 40
   :header-rows: 1

   * - ID
     - Focus
     - Links
   * - 10
     - Process second-cycle LC-MS characterization data.
     - `Open in Colab <https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/10_2_TEST_Library_characterization.ipynb>`__ | `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/10_2_TEST_Library_characterization.ipynb>`__

Learn
~~~~~

.. list-table::
   :widths: 7 53 40
   :header-rows: 1

   * - ID
     - Focus
     - Links
   * - 11
     - Run the second ML cycle and inspect improved performance and saturation of top strains.
     - `Open in Colab <https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/11_2_LEARN_Modelling_and_predictions.ipynb>`__ | `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/11_2_LEARN_Modelling_and_predictions.ipynb>`__

Publication
-----------

For the full scientific context, see the
`teemi publication <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011929>`__.
