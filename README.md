***  
**This is a code sample**
***  


### About Library

#### Introduction

This is a library of 24 electrostatic potential based 2D molecular descriptors.  
This work is based on **[Hunter, C. A. (2004). Angewandte Chemie - International Edition, 43(40), 5310–5324](https://doi.org/10.1002/anie.200301739)**    
 
**Code author/developer:** Andre Frade (with support of Dr. Patrick McCabe and Prof. Richard Cooper)   

#### Content
The library:
- contains a molecule charge neutralization method
- detects up to 46 different functional groups
- calculates up to 24 different electrostatic potential based molecular descriptors:
    - max_alpha, max_beta: maximum alpha and beta values in the molecule
    - min_alpha, min_beta: minimum alpha and beta values in the molecule
    - total_alpha, total_beta: total sum of all alpha and beta values in a molecule
    - average_alpha, average_beta: total sum of all alpha and beta values in a molecule, averaged by nr of identified functional groups
    - alpha_mw_norm, beta_mw_norm: total sum of all alpha and beta values in a molecule normalised by MolWt
    - alpha_vsa_norm, beta_vsa_norm: total sum of all alpha and beta values in a molecule normalised by vdW surface area
    - alpha_logp_norm, beta_logp_norm: total sum of all alpha and beta values in a molecule normalised by LogP
    - a_VSA0, a_VSA1, a_VSA2, a_VSA3, a_VSA4: fraction of surface area with a value of alpha between x and y 
    - b_VSA0, b_VSA1, b_VSA2, b_VSA3, b_VSA4: fraction of surface area with a value of alpha between x and y

MolWt calculation: the average molecular weight of the molecule
Surface Area calculation: Labute. P. J. Mol. Graph. Mod. _18_ 464-477 (2000)
LogP calculation: S. A. Wildman and G. M. Crippen JCICS 39 868-873 (1999)
     
#### Dependencies
In order to use this library with no problems, you run it on Python 3.7+ and must have the following dependencies installed:
- pickle
- numpy
- pandas
- rdkit
 
#### Usage
One just needs to:
- Ensure all dependencies are installed
- Change the ht_filepath and complexgroup_filepath variables with the absolute path of the ABvalues.csv and Complex_groups respectively (in _init)
- Copy and past the commented code in the class method's __init__
- Run it   

When using this library for descriptor calculation be aware:   
1. It is highly recommeded that your input molecules are in their netral charge form. The charge neutralization method works well, but it may fail. In this scenario, charged functional groups won't be accounted for and possibly, even an error might be thrown.  
2. The code will only the functional groups that are listed in ABvalues.csv. If your molecule has a group that is not on this list, the code will run, but the group won't be accounted for. Report to the original paper for instructions on how to calculate alpha and beta values, if you would like to add a new functional group to the list  



