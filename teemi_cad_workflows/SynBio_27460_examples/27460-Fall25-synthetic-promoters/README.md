# Exam project for 27460 Synthetic Biology Fall 2025 - Group 1
This repository contains the code, workflow and data used to design and generate synthetic promoters and constructs for **in vivo assembly** in ***Coprinopsis cinerea*** using teemi! This work supports the project:

**"Expression of omega-3 DHA PUFA synthesis in non-conventional host: *Coprinopsis cinerea*."**

---

## Summary
In this project a library of candidate promoters and terminators were created by cross-referencing commonly used elements (e.g., TEF1) from model organisms with homologous regions in the *C.cinerea* genome. For each region, 1000 bp were extracted and analyzed as potential promoter and terminator sequences.

The repository includes:
- **Notebooks** for:
    - Promoter extraction and synthetic promoter generation
    - Promoter and terminator screening using dTomato
    - Construct assembly and cloning
    - Primer design for PCR amplification
- **Data** containing the generated sequences
- **src** containing helper functions and models used througout the notebooks 


For the full methodology please refer to the report: **"Expression of omega-3 DHA PUFA synthesis in non-conventional host: *Coprinopsis cinerea*."**

---

## Repository Structure
```text
/data/                                  # Raw and processed sequence files  
    /constructs                         # Finalized data for full insert
    /generated_promoters                # Synthetically generated promoters
    /promoter_terminator_library        # Full list of promoter and terminators, both gathered and generated
    /insert_sequences                   # Used genetic sequences for inserts

/src/                                   # Utility functions designed for notebook use  
/notebooks/                             # Analysis and exploration  
```

---

### Usage
To explore the workflow, **open the notebooks in the `/notebooks/` directory** and follow them in order:

1. **Synthetic promoter generation**  
2. **Combinatorial library and assembly**  
3. **Construct assembly**
3. **Primer design**

---

### Authors
Alexander Videbæk and Rasmus Tøffner-Clausen - (Group 1 - 27460 Synthetic Biology, DTU (Fall 2025))