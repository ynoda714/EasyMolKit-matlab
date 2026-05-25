# Third-Party Notices

This project uses or references the following third-party resources:

## Algorithms

| Algorithm | Source | Usage |
|---|---|---|
| Wildman-Crippen LogP | Wildman, S.A.; Crippen, G.M. (1999). "Prediction of Physicochemical Parameters by Atomic Contributions". *J. Chem. Inf. Comput. Sci.* 39(5), 868-873. DOI: 10.1021/ci990307l | `emk.descriptor.calculate()` — LogP field |
| Ertl TPSA / HBA / HBD | Ertl, P.; Rohde, B.; Selzer, P. (2000). "Fast Calculation of Molecular Polar Surface Area as a Sum of Fragment-Based Contributions". *J. Med. Chem.* 43(20), 3714-3717. DOI: 10.1021/jm000942e | `emk.descriptor.calculate()` — TPSA, NumHAcceptors, NumHDonors fields |
| Lovering FractionCSP3 | Lovering, F.; Bikker, J.; Humblet, C. (2009). "Escape from Flatland: Increasing Saturation as an Approach to Improving Clinical Success". *J. Med. Chem.* 52(21), 6752-6756. DOI: 10.1021/jm901241e | `emk.descriptor.calculate()` — FractionCSP3 field |
| Morgan Fingerprint | Morgan, H.L. (1965). "The Generation of a Unique Machine Description for Chemical Structures". *J. Chem. Doc.* 5(2), 107-113. DOI: 10.1021/c160017a018 | `emk.fingerprint.morgan()` |
| MACCS Keys | MDL Information Systems. 166 public keys. | `emk.fingerprint.maccs()` |
| Tanimoto Coefficient | Tanimoto, T.T. (1958). "An Elementary Mathematical Theory of Classification and Prediction". IBM Internal Report. | `emk.similarity.tanimoto()` |
| Dice Coefficient | Dice, L.R. (1945). "Measures of the Amount of Ecologic Association Between Species". *Ecology* 26(3), 297-302. | `emk.similarity.dice()` |
| Lipinski Rule of Five | Lipinski, C.A.; Lombardo, F.; Dominy, B.W.; Feeney, P.J. (1997). "Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings". *Adv. Drug Deliv. Rev.* 23(1-3), 3-25. DOI: 10.1016/s0169-409x(96)00423-1 | `emk.filter.lipinski()` |
| Bemis-Murcko Scaffold | Bemis, G.W.; Murcko, M.A. (1996). "The Properties of Known Drugs. 1. Molecular Frameworks". *J. Med. Chem.* 39(15), 2887-2893. DOI: 10.1021/jm9602928 | `emk.mol.scaffold()` |

## Software Libraries

| Library | License | Usage |
|---|---|---|
| RDKit | BSD-3-Clause | Cheminformatics core: molecular parsing, descriptor calculation, fingerprints, similarity |
| Python | PSF License | Runtime for RDKit via MATLAB `pyenv` bridge |
| MATLAB | MathWorks Commercial | Primary development and execution environment |
| PyTorch (CPU-only) | BSD-3-Clause | Deep learning runtime for R08/R09; installed via `installExtra("torch")` |
| PyTorch Geometric (PyG) | MIT | Graph neural network framework for R08; installed via `installExtra("torch_geometric")` |
| torch_scatter / torch_sparse / torch_cluster | MIT | PyG companion packages required for message passing and sparse ops; installed alongside torch_geometric |
| HuggingFace Transformers | Apache-2.0 | Transformer model framework for R09 (ChemBERTa); installed via `installExtra("transformers")` |

## Data Sources

| Source | License | Files | Usage |
|---|---|---|---|
| PubChem (NIH/NCBI) | CC0 (Public Domain) | `data/list/everyday_chemicals.csv` | 30 household chemical SMILES |
| ChEMBL (EMBL-EBI) | CC-BY-SA 3.0 | `data/list/fda_drugs.csv`, `data/list/forensic_challenge.csv` | 200 FDA-approved drug SMILES + properties; forensic challenge set |
| RDKit Data/Pains/wehi_pains.csv | BSD-3-Clause | `data/list/pains.csv` | 480 PAINS structural alert SMARTS patterns |

### ChEMBL Attribution (CC-BY-SA 3.0)

The files `data/list/fda_drugs.csv` and `data/list/forensic_challenge.csv` contain data derived
from the ChEMBL database (https://www.ebi.ac.uk/chembl/), provided by EMBL-EBI under the
Creative Commons Attribution-ShareAlike 3.0 Unported License (CC-BY-SA 3.0).

These files are distributed under CC-BY-SA 3.0. Recipients must:
- Give appropriate credit to ChEMBL / EMBL-EBI
- Distribute any derivative of these files under the same CC-BY-SA 3.0 license

ChEMBL license: https://chembl.gitbook.io/chembl-interface-documentation/about#data-licensing

### RDKit PAINS Data Attribution (BSD-3-Clause)

The file `data/list/pains.csv` is derived from
https://github.com/rdkit/rdkit/blob/master/Data/Pains/wehi_pains.csv (BSD-3-Clause).
Original PAINS filters: Baell JB, Holloway GA. J Med Chem 53 (2010) 2719-2740.
DOI: 10.1021/jm901137j

---

All trademarks are the property of their respective owners.
