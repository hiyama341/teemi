# teemi strictosidine notebooks

This directory hosts the **Jupyter notebooks** used to generate the results in the Petersen & Levassor et al. _teemi_ paper.  
They exemplify how **literate programming** can enable full Design → Build → Test → Learn (DBTL) workflows in synthetic biology.

---

## 📖 Notebook overview

Below is a summary of the notebooks, ordered by DBTL round and stage, with a short description of each.

| DBTL Round | Notebook                                              | Stage  | Description                                                                  |
| ---------- | ----------------------------------------------------- | ------ | ---------------------------------------------------------------------------- |
| 1          | `00_1_DESIGN_Homologs`                                | DESIGN | Automatic selection of homologs from NCBI (repeatable, programmatic)         |
| 1          | `01_1_DESIGN_Promoters`                               | DESIGN | Promoter selection from RNA-seq / database with quality filters              |
| 1          | `02_1_DESIGN_Combinatorial_library`                   | DESIGN | Generation of combinatorial library via `DesignAssembly`, robot instructions |
| 1          | `03_1_BUILD_gRNA_plasmid`                             | BUILD  | Assembly of CRISPR plasmid using USER cloning                                |
| 1          | `04_1_BUILD_Background_strain`                        | BUILD  | Construction of background strain via genomic knockouts (G8H, CPR)           |
| 1          | `05_1_BUILD_Combinatorial_library`                    | BUILD  | Build of large combinatorial library (1,280 combinations)                    |
| 1          | `06_1_TEST_Library_characterisation`                  | TEST   | Data processing, LC-MS phenotyping, genotype–phenotype mapping               |
| 1          | `07_1_LEARN_Modelling_and_predictions`                | LEARN  | Use of AutoML to model and predict promising combinations                    |
| 2          | `08_2_DESIGN_Model_recommended_combinatorial_library` | DESIGN | Use model predictions to design second focused library                       |
| 2          | `09_2_BUILD_Combinatorial_library`                    | BUILD  | Construct the focused second-round library                                   |
| 2          | `10_2_TEST_Library_characterization`                  | TEST   | Process data for second-round library (LC-MS, genotyping)                    |
| 2          | `11_2_LEARN_Modelling_and_predictions`                | LEARN  | Modeling and performance comparison across cycles                            |

---

## 🛠 Usage notes

- Run the notebooks **in sequence** (from early DESIGN to final LEARN) to replicate the analytical flow.
- The notebooks assume specific datasets (genotype, phenotype) and module usage consistent with **teemi**.
- You may link each notebook filename to its corresponding file in your repo (e.g. via GitHub links) for ease of navigation.
- Use this as a **template or scaffold** for your own literate programming–based bioengineering workflows.

---

## 📌 Citation & provenance

These notebooks are integral to the reproducibility and transparency of the teemi paper. If you reuse or extend them, please cite:

> Petersen SD, Levassor L, Pedersen CM, Madsen J, Hansen LG, Zhang J, et al. (2024) _teemi: An open-source literate programming approach for iterative design-build-test-learn cycles in bioengineering_. PLoS Comput Biol 20(3): e1011929. :contentReference[oaicite:2]{index=2}

All associated data and code are available via the main teemi repository and linked data archives.
