=============================================
Data repository for G8H:CPR lirbary 
=============================================

This repository is used in relation with the Python package teemi 
for constructing microbial strains to de-bottleneck the strictosidine pathway. 
The notebooks can be found  `here <https://github.com/hiyama341/ConStrain/tree/main/colab_notebooks>`__.



teemi: Literate programming can streamline bioengineering workflows
-----------------------------------------------------------------------

.. summary-start

.. image:: https://badge.fury.io/py/teemi.svg
    :target: https://badge.fury.io/py/ConStrain

.. image:: https://github.com/hiyama341/ConStrain/actions/workflows/main.yml/badge.svg
        :target: https://github.com/hiyama341/ConStrain/actions

.. image:: https://readthedocs.org/projects/constrain/badge/?version=latest
        :target: https://teemi.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

.. image:: https://img.shields.io/github/license/hiyama341/ConStrain
        :target: https://github.com/hiyama341/ConStrain/blob/main/LICENSE

.. image:: https://img.shields.io/pypi/pyversions/teemi.svg
        :target: https://pypi.org/project/ConStrain/
        :alt: Supported Python Versions


What is Teemi?
~~~~~~~~~~~~~~~~~~

**teemi** is an easy-to-use python package with functions that
can be used in literate programming to simulate steps of a strain 
construction cycle from generating genetic parts, to designing a 
combinatorial library along with instructions for the assembly. 
A fully integrated LIMS system is presented to keep track of samples 
and allocation through both a commercial Benchling API and a low-level CSV file database. 

Here, we demonstrate the use of teemi in a complex machine learning-guided
metabolic engineering task. We envision that literate programming for biology 
can be adapted for any experimental workflow and be mixed and matched for the 
benefit of the user. As this tool is built to be flexible through its open-source
Python platform, future repetitive tasks can be automated and thus increase 
the speed at which we engineer biology. 

Curious about how you can build strains easier and faster? Head over to our `Google Colab notebooks <https://github.com/hiyama341/ConStrain/tree/main/colab_notebooks>`__
and give it a try.

Please cite our paper (link tba) if you've used teemi in a scientific publication.

.. summary-end

This repository is used in relation with the Python package teemi for constructing microbial strains strains on google colab. 


* Free software: MIT license
* Documentation: https://teemi.readthedocs.io.


Features
--------

* Combinatorial lirbary generation
* Cloning and transformation workflows
* Flowbot One instructions
* CSV based LIMS system as well as integration to Benchling
* Genotyping of microbial strains


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
