# Persistent spectral based ensemble learning (PerSpect-EL) for protein-protein binding affinity prediction

This manual is for the code implementation of paper "Persistent spectral based ensemble learning (PerSpect-EL) for protein-protein binding affinity prediction".

# Code Requirements
---
        Platform: Python>=3.6, MATLAB 2016B
        Python Packages needed: math, numpy>=1.19.5, scipy>=1.4.1, scikit-learn>=0.20.3, GUDHI 3.0.0
        
# Pipeline of PerSpect EL Models
---
![image](https://user-images.githubusercontent.com/32187437/147686016-da8e2f85-7de9-47c1-b7a5-8cabb1a4635b.png)

# Details about each step

## Protein-protein complex
Before the representation, the atom coordinates from each protein-protein complex needs to extracted accordingly to its element specific combinations e.g. mut_binding_C_C.pts.
The atom coordinates in *.pts can be generated from the original TopNetTree code 10.24433/CO.0537487.v1.

## Persistent Spectral (PerSpect) representation and Feature generation

For each protein-protein complex, the relevant element specific atom groups are used to construct the simplicial complexes to generate the Hodge Laplacians. 
```python
def hodge(args):
    # this function generates the eigenvalues from the L0 hodge laplacians in wild and mutation types from Vietoris Rips Complexes. 
    # parameter "args" is foldername format (PDBID, chainid, wildtype, resid, mutanttype) e.g. (1AK4 D A 488 G)
    # Runs the "binding_hodge.m" and "mutation_hodge.m" to generate the hodge laplacians for each complex. 
    
def alpha(args):
    # this function generates the eigenvalues from the L1 and L2 hodge laplacians in wild and mutation types from alpha complexes.
    # parameter "args" is foldername format (PDBID, chainid, wildtype, resid, mutanttype) e.g. (1AK4 D A 488 G)
    # Outputs numpy array of eigenvalues which can be used to generate persistent attributes.
 
def read_eig(args):
    # this function takes the eigenvalues from "alpha" and "hodge" functions and compute the persistent spectral attributes. 
    # parameter "args" is foldername format (PDBID, chainid, wildtype, resid, mutanttype) e.g. (1AK4 D A 488 G)
    # Outputs the persistent attributes in a numpy array after passing through the subfunction "compute_stat"

def compute_stat(dist):
    # this function takes a distribution of eigenvalues of a particular filtration parameter and computes the persistent spectral attributes. 
    # Outputs the persistent attributes of that filtration parameter. 
```
The above functions can be found in feature/skempi_feature.py and feature/run_all.py.

```python
binding_hodge.m --> Computes L0 Hodge Laplacian by constructing Vietoris Rips Complex from the atom coordinates between the binding sites. 
mutation_hodge.m --> Computes L0 Hodge Laplacian by constructing Vietoris Rips Complex from the atom coordinates between the mutation site and its neighborhood.
computeVRcomplex.m --> Constructs the VR complex from atom coordinates.
```

## Base Learner Training
---
Each persistent attribute is put into a base learner (1D-CNN) for training. 

  spectcnn_ab_prelu_v6 --> Base Learner for AB-Bind S645 (with non-binders)
  
  spectcnn_ab_prelu_Nout_v6 --> Base Learner for AB-Bind S645 (without non-binders)
  
  spectcnn_skempi_prelu_v6 --> Base Learner for SKEMPI 1131
  
  spectcnn_homology_test --> Base Learner for AB-Bind S645 Blind Homology Test Prediction
  
## Meta Learner Training 
--- 
The outputs from Base Learners are then combined and input into the Meta Learner for final prediction. 

  spectnettree_ab_prelu_v6 --> Meta Learner for AB-Bind S645 (with non-binders)
  
  spectnettree_ab_prelu_Nout_v6 --> Meta Learner for AB-Bind S645 (without non-binders)
  
  spectnettree_skempi_prelu_v6 --> Meta Learner for SKEMPI 1131
  
  spectnettree_homology_test --> Meta Learner for AB-Bind S645 Blind Homology Test Prediction

# Auxiliary Features
---     
        Software: Jackal, PDB2PQR, SPIDER2, MIBPB
        Relevant code and software can be found in https://doi.org/10.24433/CO.0537487.v1. 
        The auxiliary features for AB-Bind S645 is in models/X_ab_aux.npy.
        The auxiliary features for SKEMPI-1131 is in models/X_skempi_aux.npy.

## Cite
If you use this code in your research, please cite our paper:

* JunJie Wee and Kelin Xia, "Persistent spectral based ensemble learning (PerSpect-EL) for protein-protein binding affinity prediction", Briefings In Bioinformatics, bbac024 (2022)
