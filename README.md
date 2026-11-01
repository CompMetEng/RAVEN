<picture>
  <source media="(prefers-color-scheme: dark)" srcset="img/raven-logo-white.png">
  <img src="img/raven-logo.png" alt="RAVEN" width="400px">
</picture>

[![Release](https://img.shields.io/github/v/release/CompMetEng/RAVEN?include_prereleases)](https://github.com/CompMetEng/RAVEN/releases/)
[![Documentation](https://readthedocs.org/projects/raven-docs/badge/?version=latest)](https://raven-docs.readthedocs.io/en/latest/)
[![GitHub Discussions](https://img.shields.io/github/discussions-search?query=repo%3ACompMetEng%2FRAVEN&label=GitHub%20Discussions)](https://github.com/CompMetEng/RAVEN/discussions)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.3689518.svg)](https://doi.org/10.5281/zenodo.3689518)

**[Documentation](https://raven-docs.readthedocs.io/en/latest/)** · [Installation](https://raven-docs.readthedocs.io/en/latest/installation/index.html) · [Guide](https://raven-docs.readthedocs.io/en/latest/guide/index.html) · [API reference](https://raven-docs.readthedocs.io/en/latest/api/index.html) · [Migrate from RAVEN 2](https://raven-docs.readthedocs.io/en/latest/migrate.html) · [Citations](https://raven-docs.readthedocs.io/en/latest/references.html)

## About RAVEN

The **RAVEN** (Reconstruction, Analysis and Visualization of Metabolic Networks) Toolbox is a software suite for MATLAB for the semi-automated reconstruction of genome-scale metabolic models (GEMs). It reconstructs from published models or from the KEGG database, with gap-filling and quality control, and it simulates and analyses models and visualises simulation results and omics data. For context-specific model extraction, `ftINIT` is the method to use; `getINITModel` and `runINIT` are the original tINIT implementation, kept to reproduce models built with it.

RAVEN has a Python port, [raven-toolbox](https://github.com/CompMetEng/raven-toolbox), built on [cobrapy](https://github.com/opencobra/cobrapy). The documentation covers both toolboxes, and [RAVEN 3 and raven-toolbox](https://raven-docs.readthedocs.io/en/latest/raven3-vs-raven-toolbox.html) lists the functions that exist in one of them and the functions that return a different result in each.

This repository is where RAVEN 3 is developed. RAVEN 3 is a breaking release: the code is reorganised into functional folders, optional function arguments are passed by position, by name or both, external binaries and KEGG reference data are fetched on demand, several MATLAB toolbox dependencies are dropped, and low-value and duplicate functions are removed. The code can change until 3.0.0 is tagged. The last RAVEN 2 release is 2.11.4, maintained in [SysBioChalmers/RAVEN](https://github.com/SysBioChalmers/RAVEN). [Migrating from RAVEN 2](https://raven-docs.readthedocs.io/en/latest/migrate.html) lists every change for users who upgrade.

## Get started

1. Install RAVEN as described in the [installation instructions](https://raven-docs.readthedocs.io/en/latest/installation/raven.html). The Add-Ons manager of MATLAB installs the latest stable release, which is RAVEN 2; a RAVEN 3 pre-release is installed from this repository, and the instructions describe how. `checkRaven` verifies the installation and the solver.
2. Follow [Getting started](https://raven-docs.readthedocs.io/en/latest/guide/getting-started.html) in the guide, which has MATLAB and Python examples side by side, or open a script in the [`tutorials`](https://github.com/CompMetEng/RAVEN/tree/develop/tutorials) folder. The tutorials come from the original RAVEN paper and are described in [Legacy tutorials](https://raven-docs.readthedocs.io/en/latest/tutorials/index.html).
3. The [GEM reconstruction protocol](https://raven-docs.readthedocs.io/en/latest/protocol/index.html) reconstructs a model for *Hansenula polymorpha* from start to finish.

## Documentation and support

The installation instructions, the guide, the protocol, the tuned parameter defaults and the API reference for MATLAB and Python are at [raven-docs](https://raven-docs.readthedocs.io/en/latest/). The function help texts in this repository generate the MATLAB part of the API reference. Use [GitHub Discussions](https://github.com/CompMetEng/RAVEN/discussions) for questions and comments, and [issues](https://github.com/CompMetEng/RAVEN/issues) for bug reports.

## Cite us

If you use RAVEN in your work, please cite:

> Wang H, Marcišauskas S, Sánchez BJ, Domenzain I, Hermansson D, Agren R, Nielsen J, Kerkhoven EJ. (2018) RAVEN 2.0: A versatile toolbox for metabolic network reconstruction and a case study on _Streptomyces coelicolor_. PLoS Comput Biol 14(10): e1006541. doi:[10.1371/journal.pcbi.1006541](https://doi.org/10.1371/journal.pcbi.1006541).

All releases are archived on [Zenodo](https://doi.org/10.5281/zenodo.3689518), which allows citing the version used. If you use ftINIT, please also cite:

> Gustafsson J, Anton M, Roshanzamir F, Jörnsten R, Kerkhoven EJ, Robinson JL, Nielsen J. (2023) Generation and analysis of context-specific genome-scale metabolic models derived from single-cell RNA-Seq data. Proc Natl Acad Sci 120(6): e2217868120. doi:[10.1073/pnas.2217868120](https://doi.org/10.1073/pnas.2217868120)

[Citations](https://raven-docs.readthedocs.io/en/latest/references.html) lists the papers behind individual methods, among them the first RAVEN paper (2013), random sampling and the INIT method.

## Contributing

Contributions are welcome. Read the [contributing guidelines](https://github.com/CompMetEng/RAVEN/blob/develop/.github/CONTRIBUTING.md) first. RAVEN 3 is developed on the `develop` branch, and RAVEN 2 is maintained in [SysBioChalmers/RAVEN](https://github.com/SysBioChalmers/RAVEN).
