.. _HomologSelection: https://github.com/hiyama341/teemi/blob/main/colab_notebooks/00_DESIGN_HomologSelection.ipynb
.. _PromoterSelection: https://github.com/hiyama341/teemi/blob/main/colab_notebooks/01_DESIGN_PromoterSelection_colab.ipynb
.. _Combinatorial_library_with_DesignAssembly: https://github.com/hiyama341/teemi/blob/main/colab_notebooks/02_Combinatorial_library_with_DesignAssembly.ipynb
.. _USER_gRNA_plasmid: https://github.com/hiyama341/teemi/blob/main/colab_notebooks/03_BUILD_USER_gRNA_plasmid.ipynb
.. _Background_strain: https://github.com/hiyama341/teemi/blob/main/colab_notebooks/04_BUILD_Background_strain.ipynb
.. _croStrains_test_build: https://github.com/hiyama341/teemi/blob/main/colab_notebooks/05_BUILD_croStrains_test_build.ipynb
.. _CombinatorialLibrary: https://github.com/hiyama341/teemi/blob/main/colab_notebooks/06_BUILD_CombinatorialLibrary.ipynb
.. _LibraryCharacterisation: https://github.com/hiyama341/teemi/blob/main/colab_notebooks/07_TEST_LibraryCharacterisation.ipynb 
.. _DataAnalysisML: https://github.com/hiyama341/teemi/blob/main/colab_notebooks/08_LEARN_DataAnalysis.ipynb


Colab notebooks
===============

Welcome to the Colab notebook section. Here we will show how teemi can be used in the DBTL cycle. 



DESIGN:
-------
1. The `HomologSelection`_ notebook describes how we automatically can fetch homologs from NCBI from a query in a standardizable and repeatable way. 

2. The `PromoterSelection`_ notebook describes how promoters can be selected from RNAseq data and fetched from online database with various quality measurements implemented.

3. The `Combinatorial_library_with_DesignAssembly`_ notebook describes how a combinatorial library can be generated with the DesignAssembly class along with robot executable intructions. 

BUILD:
------
4. The `USER_gRNA_plasmid`_ notebook describes the assembly of a CRISPR plasmid with USER cloning.

5. The `Background_strain`_ notebook describes the construction of the background strain by K/O of G8H and CPR 

6. The `croStrains_test_build`_ notebook describes how a smaller combinatorial library was build and test the method before moving on to a larger library.  

7. The `CombinatorialLibrary`_ notebook shows how the first large combinatorial library was generated for 1440 possible combinations. 


TEST:
-----

8. The `LibraryCharacterisation`_ notebook describes Con data processing of LC-MS data and genotyping of the generated strains. 

LEARN:
------

9. The `DataAnalysisML`_ notebook describes how we use AutoML to predict the best combinations for a targeted second round of library construction.


For more information head over to our `publication <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/08_LEARN_DataAnalysis.ipynb>`__ describing the use of teemi in a literate programming context. 

