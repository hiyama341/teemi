Colab Notebooks
===============

These notebooks capture the larger, end-to-end strictosidine case study workflows used with
``teemi``. They are the best starting point if you want to run the published DBTL examples in
Google Colab with as little setup as possible.

Each notebook includes a direct Colab launch badge plus a GitHub source link.

.. |colab00| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open notebook 00 in Google Colab
   :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/00_1_DESIGN_Homologs.ipynb

.. |colab01| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open notebook 01 in Google Colab
   :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/01_1_DESIGN_Promoters.ipynb

.. |colab02| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open notebook 02 in Google Colab
   :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/02_1_DESIGN_Combinatorial_library.ipynb

.. |colab03| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open notebook 03 in Google Colab
   :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/03_1_BUILD_gRNA_plasmid.ipynb

.. |colab04| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open notebook 04 in Google Colab
   :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/04_1_BUILD_Background_strain.ipynb

.. |colab05| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open notebook 05 in Google Colab
   :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/05_1_BUILD_Combinatorial_library.ipynb

.. |colab06| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open notebook 06 in Google Colab
   :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/06_1_TEST_Library_characterisation.ipynb

.. |colab07| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open notebook 07 in Google Colab
   :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/07_1_LEARN_Modelling_and_predictions.ipynb

.. |colab08| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open notebook 08 in Google Colab
   :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/08_2_DESIGN_Model_recommended_combinatiorial_library.ipynb

.. |colab09| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open notebook 09 in Google Colab
   :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/09_2_BUILD_Combinatorial_library.ipynb

.. |colab10| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open notebook 10 in Google Colab
   :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/10_2_TEST_Library_characterization.ipynb

.. |colab11| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open notebook 11 in Google Colab
   :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/11_2_LEARN_Modelling_and_predictions.ipynb

Strictosidine Case: First DBTL Cycle
------------------------------------

Design
~~~~~~

.. list-table::
   :widths: 7 43 25 25
   :header-rows: 1

   * - ID
     - Focus
     - Launch
     - Source
   * - 00
     - Automatically fetch homologs from NCBI in a standardizable and repeatable way.
     - |colab00|
     - `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/00_1_DESIGN_Homologs.ipynb>`__
   * - 01
     - Select promoters from RNA-seq data and fetch them from online databases with quality checks.
     - |colab01|
     - `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/01_1_DESIGN_Promoters.ipynb>`__
   * - 02
     - Generate a combinatorial library with ``DesignAssembly`` and prepare it for downstream build planning.
     - |colab02|
     - `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/02_1_DESIGN_Combinatorial_library.ipynb>`__

Build
~~~~~

.. list-table::
   :widths: 7 43 25 25
   :header-rows: 1

   * - ID
     - Focus
     - Launch
     - Source
   * - 03
     - Assemble a CRISPR plasmid with USER cloning.
     - |colab03|
     - `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/03_1_BUILD_gRNA_plasmid.ipynb>`__
   * - 04
     - Construct the background strain by knocking out G8H and CPR.
     - |colab04|
     - `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/04_1_BUILD_Background_strain.ipynb>`__
   * - 05
     - Generate the first combinatorial library covering 1280 possible combinations.
     - |colab05|
     - `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/05_1_BUILD_Combinatorial_library.ipynb>`__

Test
~~~~

.. list-table::
   :widths: 7 43 25 25
   :header-rows: 1

   * - ID
     - Focus
     - Launch
     - Source
   * - 06
     - Process LC-MS data and genotype the generated strains.
     - |colab06|
     - `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/06_1_TEST_Library_characterisation.ipynb>`__

Learn
~~~~~

.. list-table::
   :widths: 7 43 25 25
   :header-rows: 1

   * - ID
     - Focus
     - Launch
     - Source
   * - 07
     - Use AutoML to predict the best combinations for a targeted second round of library construction.
     - |colab07|
     - `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/07_1_LEARN_Modelling_and_predictions.ipynb>`__

Strictosidine Case: Second DBTL Cycle
-------------------------------------

Design
~~~~~~

.. list-table::
   :widths: 7 43 25 25
   :header-rows: 1

   * - ID
     - Focus
     - Launch
     - Source
   * - 08
     - Translate machine-learning results into a targeted library of strains.
     - |colab08|
     - `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/08_2_DESIGN_Model_recommended_combinatiorial_library.ipynb>`__

Build
~~~~~

.. list-table::
   :widths: 7 43 25 25
   :header-rows: 1

   * - ID
     - Focus
     - Launch
     - Source
   * - 09
     - Construct the targeted second-round library of strains.
     - |colab09|
     - `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/09_2_BUILD_Combinatorial_library.ipynb>`__

Test
~~~~

.. list-table::
   :widths: 7 43 25 25
   :header-rows: 1

   * - ID
     - Focus
     - Launch
     - Source
   * - 10
     - Process second-cycle LC-MS characterization data.
     - |colab10|
     - `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/10_2_TEST_Library_characterization.ipynb>`__

Learn
~~~~~

.. list-table::
   :widths: 7 43 25 25
   :header-rows: 1

   * - ID
     - Focus
     - Launch
     - Source
   * - 11
     - Run the second ML cycle and inspect improved performance and saturation of top strains.
     - |colab11|
     - `GitHub <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/11_2_LEARN_Modelling_and_predictions.ipynb>`__

Publication
-----------

For the full scientific context, see the
`teemi publication <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011929>`__.
