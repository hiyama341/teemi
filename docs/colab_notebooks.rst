

Colab notebooks for enhancing biological workflows
--------------------------------------------------

Teemi is a software tool designed to simplify the process of
creating and running scientific workflows. By integrating 
literate programming into the workflow, the code is written 
in a way that is easy to understand, with detailed 
documentation included within the code itself. 
This approach enhances the reproducibility and 
transparency of the workflow, making it easier 
for others to understand, use, and build upon.

In the specific case of the study, teemi was used with 
literate programming to create a 
simulation-guided, iterative, and evolution-guided 
laboratory workflow for optimizing strictosidine production
in yeast. By using this approach, the researchers were able
to streamline the workflow, making it more efficient and 
effective at achieving the desired outcome.

The notebooks provided in the study are an excellent resource
for others interested in applying teemi and literate 
programming to their own bioengineering workflows.
By examining these notebooks, users can gain a better 
understanding of how to implement these tools into their 
own workflow and how to optimize the workflow for their 
specific application.

In summary, the study demonstrates the benefits of using 
teemi and literate programming in bioengineering workflows, 
providing an excellent starting point for those interested 
in optimizing their own workflows.


Strictosidine case : First DBTL cycle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**DESIGN:**

..  |Notebook 00| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Notebook 00
    :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/00_1_DESIGN_Homologs.ipynb 

..  |Notebook 01| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Notebook 01
    :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/01_1_DESIGN_Promoters.ipynb

..  |Notebook 02| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Notebook 02
    :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/02_1_DESIGN_Combinatorial_library.ipynb
    

00. Describes how we can automatically fetch homologs from NCBI from a query in a standardizable and repeatable way |Notebook 00|. 


01. Describes how promoters can be selected from RNAseq data and fetched from online database with various quality measurements implemented |Notebook 01|.



02. Describes how a combinatorial library can be generated with the DesignAssembly class along with robot executable intructions |Notebook 02|. 




**BUILD:**

..  |Notebook 03| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Notebook 03
    :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/03_1_BUILD_USER_gRNA_plasmid.ipynb


..  |Notebook 04| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Notebook 04
    :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/04_1_BUILD_Background_strain.ipynb


..  |Notebook 05| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Notebook 05
    :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/05_1_BUILD_CombinatorialLibrary_AllStrains.ipynb


03. Describes the assembly of a CRISPR plasmid with USER cloning |Notebook 03|.

04. Describes the construction of the background strain by K/O of G8H and CPR |Notebook 04|.

05. Shows how the first combinatorial library was generated for 1280 possible combinations |Notebook 05|. 



**TEST:**


..  |Notebook 06| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Notebook 06
    :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/06_1_TEST_LibraryCharacterisation.ipynb


06. Describes data processing of LC-MS data and genotyping of the generated strains |Notebook 06|.  


**LEARN:**

..  |Notebook 07| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Notebook 07
    :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/07_1_LEARN_DataAnalysis.ipynb


07. Describes how we use AutoML to predict the best combinations for a targeted second round of library construction |Notebook 07|.



Strictosidine case : Second DBTL cycle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



**DESIGN:**

..  |Notebook 08| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Notebook 08
    :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/09_BUILD_Library2.ipynb

08. Shows how results from the ML can be translated into making a target library of strains |Notebook 08|. 



**BUILD:**


..  |Notebook 09| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Notebook 09
    :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/09_2_BUILD_CombinatorialLibrary.ipynb


09. Shows the construction of a targeted library of strains |Notebook 09|. 




**TEST:**

..  |Notebook 10| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Notebook 10
    :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/10_2_TEST_Library_Characterization.ipynb



10. Describes the data processing of LC-MS data like in notebook 7 |Notebook 10|.


**LEARN:**

..  |Notebook 11| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Notebook 11
    :target: https://colab.research.google.com/github/hiyama341/teemi/blob/main/colab_notebooks/11_2_LEARN_DataAnalysisML.ipynb


11. Second ML cycle of ML showing how the model increased performance and saturation of best performing strains |Notebook 11|. 



For more information head over to our `publication <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/08_LEARN_DataAnalysis.ipynb>`__ describing the use of teemi in a literate programming context. 

