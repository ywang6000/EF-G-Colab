# EF-G Conformation Folding and Ancestral Reconstruction

This repository contains the computational workflows and analysis code supporting the paper titled:  
*Evolutionary Insights into Elongation Factor G Using AlphaFold and Ancestral Analysis*.  
**[Computers in Biology and Medicine]**, Accepted.

## Overview

This project integrates structural prediction, phylogenetic analysis, and functional validation to explore EF-G conformational evolution. Specifically, it covers:

- Identification of EF-G conformations (con1 vs. con2) using AlphaFold predictions.
- Clustering of multiple sequence alignments (MSAs) based on AlphaFold-predicted structures.
- Phylogenetic reconstruction and ancestral sequence inference.
- Functional and biophysical characterization of two reconstructed ancestral EF-G proteins.

## Contents

- `EF-G_folding.ipynb`  
  → Main notebook for folding, structural and sequence-based clustering, and visualization.

- `SEQ/wt.fasta`  
  → Wild-type EF-G query sequence.

- `SEQ/`  
  → will contain initial ColabFold outputs and parent `.a3m` alignment files generated via MMseqs2.

- `Cluster/`  
  → will Contain clustered `.a3m` files from AFCluster, with subfolders for predicted structures and analysis.

- `Cluster/folding/`  
  → will contain PDB files generated from the clustered `.a3m` files.

- `Cluster/folding/sel_pdb/`  
  → will contain PDB files that passed the pLDDT confidence threshold.

- `Cluster/folding/sel_pdb/1st_1-286_RMSD/`  
  → will contain structure clustering and PSE-based selection data and further analysis.

- `ClusterMSA_min3.py`  
  → Modified sequence clustering script using a DBSCAN-based method.

- `utils.py`  
  → Supporting utility functions adapted from:  
    *H. Wayment-Steele et al., Nature 625, 832–839 (2024), doi:10.1038/s41586-023-06832-9.*
- `a3m_file`
	→ example a3m files of con1 and con2 for WT
- `WT_list1/2.txt`
	→ example lists of a3m after structure clustering
- `one .pse`
	→ example pse file to identify con1 and con2 for WT
- `two .treefile and one .state file`
	→ example treefiles and ancestor sequence file.

## Key Tools and Dependencies
- requirement.txt
- **AlphaFold / ColabFold** — for structure prediction  
- **IQ-TREE 2** — for maximum likelihood phylogenetic reconstruction  
- **PyMOL** — for structural visualization and comparison  
- **Biopython**, **NumPy**, **Matplotlib** — for parsing, analysis, and plotting

## modularity
Most cells function as modular units, capable of carrying out specific tasks independently. 

## Reproducibility
To reproduce the computing environment:
```bash
conda env create -f myenv.yml
conda activate myenv  
```
To reproduce the full analysis pipeline:

1. Set your working directory (`default_path`).
2. Start from **cell #2** in `EF-G_folding.ipynb`, adjusting the path to your local ColabFold installation.
3. Run **cell #3** to generate the `Cluster/` directory.
4. Follow through cells #2–8 to reproduce Figures 1–2.  
   Cells #9–16 generate Figure 3, and #18–22 produce Figure 5.  
5. Cell #23 (standalone) uses Logomaker for Figure 4.  
6. Cell #17 produces Supplementary Figure S4A.

## Citation

If you use this code or find it helpful, please cite:

> Shawonur Rahaman,Jacob H. Steele,Yi Zeng,Shoujun Xu,Yuhong Wang1 (2025).  
> *Evolutionary Insights into Elongation Factor G Using AlphaFold and Ancestral Analysis*.  
> **Computers in Biology and Medicine**, Accepted.

## Contact

For questions or collaboration inquiries, please contact:  
**Yuhong Wang**  
📧 ywang60@uh.edu

