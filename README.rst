
.. image:: https://raw.githubusercontent.com/hiyama341/teemi/main/pictures/teemi1.svg?token=GHSAT0AAAAAABTYCY2LIAJ7B7ZTWBR4NMHIYZUJQ4A
  :width: 400
  :alt: teemi logo 

teemi: Literate programming can streamline bioengineering workflows
-----------------------------------------------------------------------

.. summary-start

.. image:: https://badge.fury.io/py/teemi.svg
        :target: https://badge.fury.io/py/teemi

.. image:: https://github.com/hiyama341/teemi/actions/workflows/main.yml/badge.svg
        :target: https://github.com/hiyama341/teemi/actions

.. image:: https://readthedocs.org/projects/teemi/badge/?version=latest
        :target: https://teemi.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

.. image:: https://img.shields.io/github/license/hiyama341/teemi
        :target: https://github.com/hiyama341/teemi/blob/main/LICENSE

.. image:: https://img.shields.io/pypi/pyversions/teemi.svg
        :target: https://pypi.org/project/teemi/
        :alt: Supported Python Versions

.. image:: https://codecov.io/gh/hiyama341/teemi/branch/main/graph/badge.svg?token=P4457QACUY 
        :target: https://codecov.io/gh/hiyama341/teemi

.. image:: https://img.shields.io/badge/code%20style-black-black
        :target: https://black.readthedocs.io/en/stable/


What is teemi?
~~~~~~~~~~~~~~~~~~

**teemi**, named after the Greek goddess of fairness, is a python package designed
to make microbial strain construction reproducible and FAIR (Findable, Accessible, 
Interoperable, and Reusable). With teemi, you can simulate the steps of 
a strain construction cycle, from generating genetic parts to designing 
a combinatorial library and keeping track of samples through a commercial
Benchling API and a low-level CSV file database. 
This tool can be used in literate programming to 
increase efficiency and speed in metabolic engineering tasks. 
To try teemi, visit our Google Colab notebooks. 


teemi not only simplifies the strain construction process but also offers the 
flexibility to adapt to different experimental workflows through its open-source
Python platform. This allows for efficient automation of repetitive tasks and
a faster pace in metabolic engineering.

Our demonstration of teemi in a complex machine learning-guided
metabolic engineering task showcases its efficiency 
and speed by debottlenecking a crucial step in the strictosidine pathway. 
This highlights the versatility and usefulness of this tool in various  
biological applications. To see teemi in action, check out our Google 
Colab notebooks.

Curious about how you can build strains easier and faster? 
Head over to our `Google Colab notebooks <https://github.com/hiyama341/teemi/tree/main/colab_notebooks>`__
and give it a try.

Please cite our paper (link tba) if you've used teemi in a scientific publication.

.. summary-end


Features
--------

* Combinatorial library generation
* HT cloning and transformation workflows
* Flowbot One instructions
* CSV-based LIMS system as well as integration to Benchling
* Genotyping of microbial strains
* Advanced Machine Learning of biological datasets with the AutoML `H2O <https://docs.h2o.ai/h2o/latest-stable/h2o-docs/automl.html>`__
* Workflows for selecting enzyme homologs
* Promoter selection workflows from RNA-seq datasets
* Data analysis of large LC-MS datasets along with workflows for analysis


..  image:: https://raw.githubusercontent.com/hiyama341/teemi/main/pictures/DBTL_teemi_with_modules1.svg
    :width: 800
    :alt: Overview of teemi's features throughout the DBTL cycle. 

Getting started
~~~~~~~~~~~~~~~
To get started with making microbial strains in an HT manner please follow the steps below: 

1. Install teemi. You will find the necessary information below for installation.

2. Check out our `notebooks <https://github.com/hiyama341/teemi/tree/main/colab_notebooks>`__ for inspiration to make HT strain construction with teemi.

3. You can start making your own workflows by importing teemi into either Google colab or Jupyter lab/notebooks.


Colab notebooks
---------------
As a proof of concept we show how teemi and literate programming can be used to streamline bioengineering workflows.
These workflows should serve as a guide or a help to build your own workflows and thereby harnessing the power of literate programming with teemi. 

Specifically, in this study we present how teemi and literate programming to build simulation-guided, iterative,
and evolution-guided laboratory workflows for optimizing strictosidine production in yeast.

Below you can find all the notebooks developed in this work. 
Just click the Google colab badge to start the workflows. 

First DBTL cycle
~~~~~~~~~~~~~~~~



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



Second DBTL cycle
~~~~~~~~~~~~~~~~~



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



Installation
~~~~~~~~~~~~

.. installation-start

Use pip to install teemi from `PyPI <https://pypi.org/project/teemi/>`__.

::

    $ pip install teemi


If you want to develop or if you cloned the repository from our `GitHub <https://github.com/hiyama341/teemi/>`__
you can install teemi in the following way.

::

    $ pip install -e <path-to-teemi-repo>  


You might need to run these commands with administrative
privileges if you're not using a virtual environment (using ``sudo`` for example).
Please check the `documentation <https://teemi.readthedocs.io/en/latest/installation.html#>`__
for further details.

.. installation-end

Documentation and Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~

Documentation is available on through numerous Google Colab notebooks with
examples on how to use teemi and how we use these notebooks for strain
construnction. The Colab notebooks can be found here 
`teemi.notebooks <https://github.com/hiyama341/teemi/tree/main/colab_notebooks>`__. 

* Documentation: https://teemi.readthedocs.io.


Contributions
~~~~~~~~~~~~~

Contributions are very welcome! Check our `guidelines <https://teemi.readthedocs.io/en/latest/contributing.html>`__ for instructions how to contribute.


License
~~~~~~~
* Free software: MIT license

Credits
-------
- This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter

.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

- teemis logo was made by Jonas Krogh Fischer. Check out his `website <http://jkfischerproductions.com/kea/portfolio/index.html>`__. 