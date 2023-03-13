.. image:: https://raw.githubusercontent.com/hiyama341/teemi/main/pictures/teemi_logo.svg
  :width: 400
  :alt: teemi logo 

======================================
Welcome to teemi's documentation!
======================================

Our mission is to revolutionize the way in which we make biology.
How do we do that? We pursue reproducible high-throughput strain construction by adopting literate programming with teemi.
We hope you want to be part of the revolution.
Below you will find documentation for how to install and use teemi.


Suppose you are interested in examples, head over to our `notebooks <https://github.com/hiyama341/teemi/tree/main/colab_notebooks>`__. 
Feel free to reach out if something is unclear.


Features of teemi
---------------------

- teemi/

    - design/
        - combinatorial_design.py
        - teselagen_helpers.py
        - cloning.py
        - retrieve_gene_homologs.py
        - fetch_sequences.py
    
    - build/
        - transformation.py
        - containers_wells_picklists.py
        - robot_assembly.py
        - PCR.py

    - test/
        - genotyping.py

    - learn/
        - plotting.py
        - auto_ml.py

    - lims/
        - benchling_api.py
        - csv_database.py  

    - utils.py


    
Getting started
~~~~~~~~~~~~~~~
To get started with making microbial strains in an HT manner please follow the steps below: 

1. Install teemi. You will find the necessary information below under installation.

2. Check out our `notebooks <https://github.com/hiyama341/teemi/tree/main/colab_notebooks>`__ for inspiration to make HT strain construction with teemi.

3. You can start making your own workflows by importing teemi into either Google colab or Jupyter lab/notebooks.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   readme
   installation
   usage
   contributing
   authors
   history
   modules
   colab_notebooks
   teemi.design
   teemi.build
   teemi.test
   teemi.learn
   teemi.lims
   api

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
