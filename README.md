# Persistent spectral based ensemble learning (PerSpect-EL) for protein-protein binding affinity prediction

This manual is for the code implementation of paper "Persistent spectral based ensemble learning (PerSpect-EL) for protein-protein binding affinity prediction".

# Code Requirements
---
        Platform: Python>=3.6
        Packages needed: math, numpy>=1.18.1, scipy>=1.4.1, scikit-learn>=0.22.1
        
# Flow of PerSpect (EL) models
---
Protein-protein interaction complex  -->  PerSpect representation and Feature generation -->  Base Learner Training --> Meta Learner

# Details about each step

## Protein-protein interaction complex
Before the representation, the atom coordinates from each protein-protein complex needs to extracted accordingly.  

## Persistent Spectral (PerSpect) representation and Feature generation
For each protein-protein complex, the relevant element specific atom groups are used to construct the simplicial complexes to generate the Hodge Laplacians. 
```python
def hodge(args):
    # this function generates the eigenvalues from the L0 hodge laplacians in wild and mutation types from Vietoris Rips Complexes. 
    # parameter "args" is foldername format (PDBID, chainid, wildtype, resid, mutanttype) e.g. (1AK4 D A 488 G)
    
def alpha(args):
    # this function generates the eigenvalues from the L1 and L2 hodge laplacians in wild and mutation types from alpha complexes.
    # parameter "args" is foldername format (PDBID, chainid, wildtype, resid, mutanttype) e.g. (1AK4 D A 488 G)
    # Outputs numpy array of features which can be splitted into the L1 and L2 pers
```

# Auxiliary Features Generation
---     
        Software: Jackal, PDB2PQR, SPIDER2, MIBPB
        Relevant code and software can be found in https://doi.org/10.24433/CO.0537487.v1. 

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

