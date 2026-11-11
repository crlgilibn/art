# De Sitter Horizon Galactic Gravity

*Parameter-Free Planck–Bose Horizon Mechanics, QUMOND Action, and the SPARC Radial Acceleration Relation*

---

## 1. Repository Overview

This repository provides a reproducible, zero-parameter benchmark of the galactic Radial Acceleration Relation (RAR) across the SPARC database (Lelli, McGaugh & Schombert 2016, 2017), deriving the de Sitter horizon acceleration scale $a _ 0 = \frac{c\thinspace H _ 0}{2\pi}{}$, the thermal Planck–Bose transition $g _ {\text{eff}}(g _ {\text{bar}}) = g _ {\text{bar}} / \bigl(1 - \exp(-\sqrt{g _ {\text{bar}}/a _ 0})\bigr){}$, and its exact closed-form polylogarithm QUMOND variational principle $\mathcal{F}(y){}$.

| File / Directory | Role | Description |
| :--- | :--- | :--- |
| **[`SPARC.md`](SPARC.md)** | **Technical Article** | Full derivation of $a _ 0 = \frac{c H _ 0}{2\pi}{}$, the closed-form polylogarithm QUMOND action $\mathcal{F}(y){}$ ($\text{Li} _ 2, \text{Li} _ 3, \text{Li} _ 4$), Solar-System Cassini immunity ($10^{-342}{}$), empirical SPARC benchmark ($\chi^2 _ \nu = 1{,}6203$, $\Delta\text{BIC} = -45{,}6$), and discussion of relativistic lensing, External Field Effect (EFE), and cluster limits. |
| **[`study/sparc_octonion_test.py`](study/sparc_octonion_test.py)** | **Verification Script** | Self-contained Python 3 benchmark (standard library only, `0` external dependencies) verifying all RAR bins, unbinned points with official error bars $(\sigma _ {\log g _ {\text{bar}}}, \sigma _ {\log g _ {\text{obs}}}){}$, 175 rotation curves, QUMOND Euler–Lagrange derivatives, Cassini bounds, weak lensing, and `RC100` $\text{high-}z$ disks. |
| **`study/SPARC_RARbins.mrt`** | **Dataset ($14$ bins)** | Binned Radial Acceleration Relation data from McGaugh, Lelli & Schombert (2016, *Phys. Rev. Lett.* **117**, 201101). |
| **`study/SPARC_RAR_all2693.mrt`** | **Dataset ($2\thinspace 693$ pts)** | Spatially resolved acceleration pairs ${}(g _ {\text{bar}}, g _ {\text{obs}}){}$ and uncertainties across $153$ SPARC galaxies (Lelli et al. 2017, *Astrophys. J.* **836**, 152). |
| **`study/SPARC_MassModels.mrt`** | **Dataset ($175$ galaxies)** | Full radial mass models $V _ {\text{obs}}(R), V _ {\text{gas}}(R), V _ {\text{disk}}(R), V _ {\text{bul}}(R){}$ across $3\thinspace 391$ points (Lelli et al. 2016, *Astron. J.* **152**, 157). |
| **[`study/fig1_sparc_rar_and_rotcurves.png`](study/fig1_sparc_rar_and_rotcurves.png)** | **Figure** | Visual comparison of the parameter-free RAR curve against Newton/GR and four representative SPARC galaxy rotation curves. |

---

## 2. Benchmark Summary

All models below are evaluated on the unmodified SPARC datasets with fixed universal stellar mass-to-light ratios $\Upsilon _ {\text{disk}} = 0{,}50\thinspace M _ \odot/L _ \odot$ and $\Upsilon _ {\text{bul}} = 0{,}70\thinspace M _ \odot/L _ \odot$ and official observational uncertainties:

| Model | Free Fit Parameters $k$ | 14-Bin RMS ${}[\text{dex}]{}$ | All $2\thinspace 693$ Points RMS ${}[\text{dex}]{}$ | Mean Bias ${}[\text{dex}]{}$ | Official $\chi^2 _ \nu$ ($2\thinspace 693$ pts) | Cassini Saturn Bound $\lvert \delta g / g _ N\rvert < 2\times 10^{-9}{}$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Pure Newton / Einstein GR** *(baryons only)* | **0** | 0.4746 | 0.5140 | -0.4305 | 21.6322 | Satisfied ($0$) |
| **Milgrom Standard MOND** (*1983*) | **1** (*fitted* $a _ 0$) | 0.0658 | 0.1484 | -0.0582 | 2.0574 | Satisfied ($1{,}3 \times 10^{-12}{}$) |
| **Famaey–Binney Simple MOND** (*2005*) | **1** (*fitted* $a _ 0$) | 0.0227 | 0.1331 | +0.0094 | 1.6324 | **Excluded ($1{,}6 \times 10^{-6}{}$)** |
| **Verlinde Emergent Gravity** (*2017,* $c H _ 0/6$) | **0** | 0.0553 | 0.1439 | -0.0513 | 1.9146 | Requires extra screening |
| **Planck–Bose Horizon (Planck $H _ 0 = 67{,}4$)** | **0** | **0.0256** | **0.1339** | -0.0179 | **1.6348** | **Satisfied ($10^{-342}{}$)** |
| **Planck–Bose Horizon (SH0ES $H _ 0 = 73{,}04$)** | **0** | **0.0186** | **0.1329** | **-0.0057** | **1.6203** | **Satisfied ($10^{-342}{}$)** |
| *Empirical 1-Parameter Fit* (*McGaugh 2016*) | **1** (*fitted* $g _ \dagger$) | *0.0172* | *0.1329* | *+0.0036* | *1.6343* | *Satisfied (*$10^{-318}{}$) |

---

## 3. Reproducible Verification

Run the self-contained verification suite directly from the repository root:

```bash
python3 study/sparc_octonion_test.py --check
```

For the full mathematical derivation, observational comparisons, and discussion of open theoretical questions, see **[`SPARC.md`](SPARC.md)**.
