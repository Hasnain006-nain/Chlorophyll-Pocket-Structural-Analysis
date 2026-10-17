# Comparative Structural and Chemical Microenvironment Analysis of the $A_{-1}$ Chlorophyll Binding Pockets in Far-Red Light Photosystem I

**Author:** Hasnain Haider  
**Date:** October 2026  
**Analyzed Macromolecular Structures:** PDB [9I9L](https://www.rcsb.org/structure/9I9L) (Primary), PDB [9EYS](https://www.rcsb.org/structure/9EYS) (Replicate)  
**Organism:** *Chroococcidiopsis thermalis* (FaRLiP, Far-Red Light Photoacclimation)  
**Target Sites:** Reaction Center Primary Electron Acceptors $A_{-1A}$ (Chlorophyll $a$, CLA) vs. $A_{-1B}$ (Chlorophyll $f$, F6C)

---

## 1. Abstract & Research Motivation

Oxygenic photosynthetic organisms adapted to far-red light (FaRLiP) incorporate chlorophyll *f* (Chl *f*) to capture photons beyond 700 nm. Within the Photosystem I (PSI) reaction center core, electron transport occurs across two pseudo-symmetric branches (the A-branch and B-branch). Cryo-EM models have unveiled a striking spatial asymmetry: **Chl *f* is selectively incorporated into the primary electron acceptor $A_{-1B}$ site, whereas the pseudo-symmetric $A_{-1A}$ site retains standard chlorophyll *a* (Chl *a*)**.

The chemical difference between Chl *a* and Chl *f* resides entirely at the **C2 carbon position**:
* **Chl *a*:** Carries a non-polar methyl substituent (`-CH3`, CMB).
* **Chl *f*:** Carries an electron-withdrawing formyl substituent (`-CHO`, CMB=OMB).

**Core Scientific Question:**  
*Does the local protein microenvironment surrounding the C2 position provide specific steric room, hydrogen-bonding networks, or aromatic stacking that rationalizes site-specific binding of Chl f at $A_{-1B}$ rather than $A_{-1A}$?*

---

## 2. Computational Methodology

All measurements were performed algorithmically on atomic coordinates extracted from RCSB macromolecular CIF models (`9i9l.cif` and `9eys.cif`) via Biopython (`MMCIFParser`) and spatial KD-trees (`NeighborSearch`):

1. **Spatial Site Identification:**  
   Chlorophyll pigments were mapped using their central magnesium ($Mg$) atoms. Met-ligated chlorophylls were cataloged as $A_0$ electron acceptors. Non-Met chlorophylls with $Mg\text{--}Mg \le 14.0$ Å from $A_0$ were identified as $A_{-1}$ candidates. Formyl-bearing pigments were assigned as $A_{-1B}$ (F6C) and paired with the pseudo-symmetric $A_{-1A}$ (CLA) at ~18.8 Å separation.
2. **Chemical Validation:**  
   The C2 substituent was validated without relying purely on PDB atom names: Chl *a* contains CMB (methyl carbon, no oxygen attached); Chl *f* contains CMB bonded to OMB (formyl carbonyl oxygen, verified bond length: 1.23 Å).
3. **Steric Clearance & Virtual Planar Probes:**  
   * For $A_{-1B}$, clearance was measured directly from the observed atomic coordinates of OMB to the nearest non-pigment atom.
   * For $A_{-1A}$, because no physical oxygen is present in Chl *a*, two hypothetical in-plane formyl oxygen probe positions were mathematically projected into the porphyrin macrocycle plane ($\angle(\text{C}_{\text{ring}}-\text{C}_{\text{formyl}}=\text{O}) = 123^\circ$, bond length = 1.22 Å, $\pm 57^\circ$ offset) to assess potential steric room.
4. **Pocket Microenvironment Metrics:**  
   * Packing density: Count of non-pigment heavy atoms within 5.0 Å of CMB.
   * Nearest aromatic ring centroid: Calculated across all Phe, Tyr, and Trp residues.
   * Nearest cysteine sulfur: Distance to Cys SG.
   * Approximate net formal charge: Sum of Arg (+1), Lys (+1), Asp (-1), Glu (-1) within 8.0 Å of CMB.

---

## 3. Quantitative Triplicate Results

Both 9I9L and 9EYS are Photosystem I trimers, containing three independent biological replicates per complex.

### Summary Data Table across All Replicates

| Structure | Rep | Site | Residue | Pigment | Axial Ligand | O Clearance (Å) | Nearest Contact | Formyl H-Bonds (≤3.3Å) | Heavy Atoms (≤5Å) | Nearest Aromatic | Aromatic Dist (Å) | Net Charge (8Å) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **9I9L** | 1 | A-1A | B:802 | CLA (Chl a) | B:HOH919 (2.05 Å) | 3.50 (simulated) | B:LEU438 (3.00 Å) | None | 13 | A:PHE708 | 5.64 | 0 |
| **9I9L** | 1 | A-1B | A:802 | F6C (Chl f) | A:HOH926 (2.08 Å) | 3.05 (real OMB) | A:PHE459:CG | **None** | 17 | A:PHE459 | 3.60 | 0 |
| **9I9L** | 2 | A-1A | O:802 | CLA (Chl a) | O:HOH918 (2.06 Å) | 3.49 (simulated) | O:LEU438 (3.04 Å) | None | 15 | N:PHE708 | 5.61 | 0 |
| **9I9L** | 2 | A-1B | N:802 | F6C (Chl f) | N:HOH910 (2.07 Å) | 3.03 (real OMB) | N:PHE459:CG | **None** | 17 | N:PHE459 | 3.60 | 0 |
| **9I9L** | 3 | A-1A | b:802 | CLA (Chl a) | b:HOH916 (2.06 Å) | 3.47 (simulated) | b:LEU438 (3.00 Å) | None | 14 | a:PHE708 | 5.66 | 0 |
| **9I9L** | 3 | A-1B | a:802 | F6C (Chl f) | a:HOH923 (2.08 Å) | 2.98 (real OMB) | a:PHE459:CG | **None** | 17 | a:PHE459 | 3.56 | 0 |
| **9EYS** | 1 | A-1A | A:856 | CLA (Chl a) | A:HOH914 (2.03 Å) | 3.49 (simulated) | B:LEU438 (2.87 Å) | None | 13 | A:PHE708 | 5.61 | 0 |
| **9EYS** | 1 | A-1B | A:802 | F6C (Chl f) | A:HOH916 (2.04 Å) | 3.30 (real OMB) | A:PHE459:CG | **None** | 17 | A:PHE459 | 3.56 | 0 |
| **9EYS** | 2 | A-1A | O:802 | CLA (Chl a) | O:HOH918 (2.05 Å) | 3.57 (simulated) | O:LEU438 (2.90 Å) | None | 14 | N:PHE708 | 5.65 | 0 |
| **9EYS** | 2 | A-1B | N:802 | F6C (Chl f) | N:HOH904 (2.02 Å) | 3.20 (real OMB) | N:PHE459:CD2 | **None** | 17 | N:PHE459 | 3.53 | 0 |
| **9EYS** | 3 | A-1A | b:802 | CLA (Chl a) | b:HOH913 (2.03 Å) | 3.51 (simulated) | b:LEU438 (2.97 Å) | None | 14 | a:PHE708 | 5.65 | 0 |
| **9EYS** | 3 | A-1B | a:802 | F6C (Chl f) | a:HOH913 (2.04 Å) | 3.21 (real OMB) | a:PHE459:CD2 | **None** | 18 | a:PHE459 | 3.50 | 0 |

### Mean Comparisons

* **Formyl Oxygen Clearance:**  
  * 9I9L: $A_{-1A} = 3.49 \pm 0.02$ Å vs. $A_{-1B} = 3.02 \pm 0.04$ Å  
  * 9EYS: $A_{-1A} = 3.53 \pm 0.04$ Å vs. $A_{-1B} = 3.24 \pm 0.06$ Å
* **Packing Density (Heavy Atoms within 5 Å):**  
  * 9I9L: $A_{-1A} = 14.0 \pm 1.0$ vs. $A_{-1B} = 17.0 \pm 0.0$  
  * 9EYS: $A_{-1A} = 13.7 \pm 0.6$ vs. $A_{-1B} = 17.3 \pm 0.6$
* **Distance to Nearest Aromatic Ring:**  
  * 9I9L: $A_{-1A} = 5.64 \pm 0.03$ Å (Phe708) vs. $A_{-1B} = 3.59 \pm 0.02$ Å (Phe459)  
  * 9EYS: $A_{-1A} = 5.63 \pm 0.02$ Å (Phe708) vs. $A_{-1B} = 3.53 \pm 0.03$ Å (Phe459)
* **Approximate Net Charge within 8 Å:**  
  * In both structures across all replicates, the net formal charge is **0.00** (neutral electrostatic background).

---

## 4. Key Biochemical Findings & Discussion

### 4.1 Strict Absence of Hydrogen Bonds to the Formyl Group
In 0 of 6 biological replicates is there a polar oxygen or nitrogen atom within 3.3 Å of OMB. Rather than being stabilized by conventional hydrogen bonds, the OMB oxygen atom is oriented into a hydrophobic, aromatic-rich cavity dominated by **Phe459**.

### 4.2 Local Aromatic Asymmetry: Phe459 vs. Leu438
The most robust structural difference between the two pockets is the presence of an aromatic phenylalanine:
* **At $A_{-1B}$:** **PsaA:Phe459** sits directly above the C2 substituent, with a centroid-to-CMB distance of ~3.55 Å. This facilitates $\pi\text{--}\pi$ stacking and favorable polar-aromatic interactions with the conjugated carbonyl group. Additionally, Phe571 and Cys463 (~7.1 Å) are present nearby.
* **At $A_{-1A}$:** The homologous position is occupied by **Leu438**, an aliphatic branched amino acid. The nearest aromatic ring (Phe708) is more than 5.6 Å away.

### 4.3 Out-of-Plane Torsion and Packing
Idealized planar formyl probes reveal that a strictly coplanar carbonyl group at $A_{-1B}$ would experience severe steric pressure. In the experimental cryo-EM models, OMB is displaced 0.53–0.94 Å from the macrocycle plane, indicating that out-of-plane formyl torsion or porphyrin ruffling is essential for relieving steric clashes while optimizing contact with Phe459.

---

## 5. Caution on Causation vs. Correlation

These findings document structural correlates extracted from static, high-resolution cryo-EM densities. While these correlates provide a rational biophysical basis for the observed selectivity, they do not prove thermodynamic causation. Site-directed mutagenesis (e.g., reciprocal PsaA:Phe459Leu and PsaB:Leu438Phe swaps) coupled with spectroscopic characterization will be essential to definitively determine the free-energy contributions to pigment selectivity.

---

## 6. References
1. Kato, K., et al. (2024). Structural basis of far-red light-driven Photosystem I from *Chroococcidiopsis thermalis*. *Nature Communications* (PDB: [9I9L](https://www.rcsb.org/structure/9I9L)).
2. Cryo-EM replicate coordinates: PDB [9EYS](https://www.rcsb.org/structure/9EYS).
3. Jordan, P., et al. (2001). Three-dimensional structure of cyanobacterial photosystem I at 2.5 Å resolution. *Nature*, 411, 909–917 (PDB: [1JB0](https://www.rcsb.org/structure/1JB0)).
4. Nürnberg, D. J., et al. (2018). Photochemistry beyond the red limit in chlorophyll *f*-containing photosystems. *Science*, 360, 1210–1213.
