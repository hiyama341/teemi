# Roadmap

## Vision

`teemi` should become a practical home for synthetic biology and strain engineering workflows:

- a Python package for designing, planning, and documenting DBTL workflows
- a notebook and workflow gallery for inspiration
- a place where users can adapt real examples for their own host, pathway, or library

The long-term goal is to make `teemi` useful both for:

- computational workflow design
- practical lab planning for build and validation work

## Near-Term Priorities

### 1. Strengthen the workflow-facing side of the repo

- Expand `teemi_cad_workflows` and `colab_notebooks` as first-class examples
- Create a workflow gallery in the docs
- Add more lightweight project notebooks that show how to plan real strain engineering projects
- Keep examples biologically motivated and easy to adapt

### 2. Improve documentation structure

- Keep the README focused on onboarding
- Use dedicated docs/pages for roadmap, releasing, and workflow discovery
- Make the docs less API-first and more workflow-first
- Ensure docs build reliably on Read the Docs

### 3. Add protein engineering support

Add an initial protein engineering layer under `teemi.design` focused on practical sequence design rather than structure-heavy methods.

Planned capabilities:

- homolog mining and panel selection
- codon-optimized expression panel generation
- saturation mutagenesis library design
- signal peptide swapping
- tag and linker design
- sequence domestication and basic synthesis QC

### 4. Modernize build and lab-planning functions

Refactor older build helpers into clearer, more reusable planning interfaces.

Focus areas:

- standardized primer outputs
- standardized edit plans
- standardized validation plans
- cleaner assembly and plate-planning outputs
- clearer separation between planning logic and export/integration logic

### 5. Reduce coupling to specific external systems

Over time, keep the core package less dependent on one lab stack.

Likely direction:

- core planning code in `teemi.design`, `teemi.build`, `teemi.test`, and `teemi.learn`
- external system adapters moved toward integration-specific modules

## Proposed Package Direction

This is directional, not a strict immediate refactor.

```text
teemi/
  design/
    protein_engineering.py
    variant_library.py
    sequence_qc.py
  build/
    assembly.py
    plate_layout.py
    liquid_handling.py
  test/
    validation.py
    screening.py
  learn/
    ranking.py
    analysis.py
  integrations/
    benchling.py
    vendors.py
    lims.py
```

## Planned Workflow Themes

### Strain engineering

- promoter swap overexpression libraries
- multiplex CRISPR editing planners
- terminator libraries for expression tuning
- reporter knock-in workflows

### Protein engineering

- homolog mining to codon-optimized expression panels
- saturation mutagenesis libraries
- signal peptide swapping
- tagged enzyme library design

### Build and validation

- PCR planning for multi-part assemblies
- assembly planning and plate maps
- genotyping planners
- expected amplicon validation workflows

## UV Adoption

The repo should keep `pip` support for users, but `uv` should become the preferred developer workflow.

Planned steps:

- document `uv` setup in the README
- add `uv` examples for installing and running tests
- optionally migrate CI install/test steps to `uv` later

## Release Hygiene

The release process should stay simple and documented:

- test before release
- tag with plain version numbers, e.g. `1.0.3`
- let GitHub Actions publish tagged releases to PyPI

See [RELEASING.md](./RELEASING.md).
