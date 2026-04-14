# Quickstart

`teemi` works best in notebooks or Python scripts where design, build, test, and learn steps live together.

## Import the package

```python
import teemi
```

## Explore the main namespaces

- `teemi.design` for sequence design, cloning, and workflow planning helpers.
- `teemi.build` for assembly, PCR, and transformation-oriented utilities.
- `teemi.test` for genotyping and data wrangling helpers.
- `teemi.learn` for analysis and plotting utilities.
- `teemi.lims` for sample and sequence tracking integrations.

## Example: combinatorial library planning

One of the quickest ways to get started is to assemble a small combinatorial design workflow
from sequence files and inspect the generated primer plan.

```python
from teemi.design.fetch_sequences import read_genbank_files
from teemi.design.combinatorial_design import DesignAssembly
from pydna.dseqrecord import Dseqrecord

path = "../data/genetic_parts/G8H_CYP_CPR_PARTS/"
pCPR_sites = [Dseqrecord(seq) for seq in read_genbank_files(path + "CPR_promoters.gb")]
CPR_sites = [Dseqrecord(seq) for seq in read_genbank_files(path + "CPR_tCYC1.gb")]

library = DesignAssembly(
    [pCPR_sites, CPR_sites],
    pad=[Dseqrecord("TCGGTC")],
    position_of_pads=[1],
    overlap=35,
    target_tm=55,
)

library.primer_list_to_dataframe()
```

This example is adapted from the repository README and shows the basic flow:
load sequence parts, define the assembly configuration, and inspect the resulting primer table.

Example output:

```text
      id anneals to sequence  annealing temperature  length  price(DKK)
0   P001      pMLS1      ...                  56.11      20        36.0
1   P002      pMLS1      ...                  56.18      49        88.2
2    ...        ...      ...                    ...     ...         ...
```

## Continue with examples

- Start with the [workflow gallery](workflow_gallery) for concise notebook-style examples.
- Browse [Colab notebooks](colab_notebooks) for larger end-to-end case studies.
- Use the [API reference](api/index.md) when you want package-level documentation.
- For a fuller walkthrough, see the combinatorial library example in the repository [README](https://github.com/hiyama341/teemi/blob/main/README.rst).
