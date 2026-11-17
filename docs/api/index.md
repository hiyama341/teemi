# API reference

This reference focuses on stable package entry points that are safe to document during a lightweight docs build.
The active lab-workflow surface is centered on `teemi.build`, while `teemi.lims` remains available as a legacy
compatibility namespace for older notebook workflows and integrations.

```{eval-rst}
.. autosummary::
   :toctree: generated

   teemi
   teemi.build
   teemi.design
   teemi.learn
   teemi.lims
   teemi.test
   teemi.utils
```

## CRISPR design

The guide-finding and base-editing modules behind the Cas9, Cas12a, Cas3,
CRISPRi and CRISPR-BEST workflows. These carry their reference documentation in
the modules themselves, so they are pulled in directly rather than through the
summary table above.

### Guide finding and off-targets

```{eval-rst}
.. automodule:: teemi.design.crispr_cas
```

### CRISPRi

```{eval-rst}
.. automodule:: teemi.design.crispri
```

### Base editing (CRISPR-BEST)

```{eval-rst}
.. automodule:: teemi.design.crispr_best
```
