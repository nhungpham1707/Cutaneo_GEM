**Genome-scale constraint-based metabolic model of* Cutaneotrichosporon curvatum*  **

**Folder structure description: **
A genome-scale constraint-based metabolic model of* Cutaneotrichosporon curvatum*, iNP636, was constructed using *Yarrowia lipolytica* model iNL895 as a template. The model was then used to study lipid synthesis in *C.curvatus*.

Start date: Nov-2016

The data folder contains data to draft, curate and validate the model. 
The script folder contains codes to build model (construction folder), Validation (same name folder), and Simulation (same name folder). 
In the construction folder, BuildCryptoGEM is the master script to build *C.curvatum *model. This function calls all other functions in the current folder to build the model.

More details are explained in the script. 

Nhung, 07th Feb 2019

**Instruction to use model**

**Required software **


*  Matlab2015b and above (some functions in COBRA toolbox had not been designed for older Matlab version)
*  CobraToolbox version 2.0.6
*  Or other software or toolboxes that support genome-scale metabolic models in SBML format 

**Model files**
iNP636.xml and iNP636.mat. 

**Scripts**
