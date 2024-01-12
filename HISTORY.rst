History
-------
0.3.4 (2024-10-01)

From issues #16 - Including dependencies in Pypi package
- Have updated setup.py file and dependencies are included in Pypi directly. 

Thanks to @manulera for raising this issue.   

0.3.3 (2023-11-08)

- Updated readme file with extra examples

- Fixed setup.py file for installation of development packages like: pip install teemi[dev]


0.3.2 (2023-01-08)

This release features a re-factored DesignAssembly class with: 

- Simplified methods i.e. redundant methods have been removed.

- The ability to add more than one pad, which can be used to make constructs with overlapping ends for for plasmid cloning.
 

0.3.1 (2023-31-07)
- This release failed due to a bug in the readme file.


0.3.0 (2023-22-06)
~~~~~~~~~~~~~~~~~~

* New submodules: gibson_cloning

This module is used to perform simple Gibson cloning workflows. 
While the addition of the "gibson_cloning" submodule is an exciting development, this module is still a work in progress.
Next, a golden gate module. Keep posted on the progress. 


0.2.0 (2023-31-05)
~~~~~~~~~~~~~~~~~~

* New submodules: CRISPRsequencecutter, sequence_finder. 

CRISPRSequenceCutter is a dataclass that is used to cut DNA through CRISPR-cas9 double-stranded break.
SequenceFinder is a dataclass that finds upstream and downstream sequences from a sequence input, annotates them and saves them.

0.1.0 (2023-01-02)
~~~~~~~~~~~~~~~~~~

* First release on PyPI.


