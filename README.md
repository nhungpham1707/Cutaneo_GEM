# Project: Constructing Genome-scale constraint-based metabolic model of *Cutaneotrichosporon oleaginosus* 
A genome-scale constraint-based metabolic model of* Cutaneotrichosporon oleaginosus*, iNP636_Coleaginosus_ATCC20509, was constructed using *Yarrowia lipolytica* model iNL895 as a template. 
The model was then used to study lipid synthesis in *C. oleaginosus*.
* Start date: 20-Nov-2016
* Complete : 20-Dec-2019
* Nhung, 07th Feb 2019

The model focus on lipid metabolism. A curated lipid metabolic pathway is constructed as in figure below

<p>
<img src="https://github.com/nhungpham1707/Cutaneo_GEM/blob/master/Results_Figures/lipid_v01.pdf" width="400" height="300" alt>
</p>
<p>
    <em>Curated lipid metabolic pathway in *C. oleaginosus*<em>.
        </p>


## Folder structure description
The data folder contains data to draft, curate and validate the model. 
The script folder contains codes to build model (construction folder), Validation (same name folder), and Simulation (same name folder). 
In the construction folder, BuildCutaneoGEM is the master script to build *C. oleaginosus *model. This function calls all other functions in the current folder to build the model.

More details are explained in the script. 

## Instruction to use model 

**Required software**


*  Matlab2015b and above (some functions in COBRA toolbox had not been designed for older Matlab version)
*  CobraToolbox version 2.0.6
*  Or other software or toolboxes that support genome-scale metabolic models in SBML format

**Model files**
iNP636_Coleaginosus_ATCC20509.xml and iNP636_Coleaginosus_ATCC20509.mat. 

**Scripts**
- Clone this repository to local computer and add it to malab path
* To build model  - call BuildCutaneoGEM - the function will generate model from the template and carry all curation 
* To simulate - there are ready-to-use scripts with documentation for how to use and what simulation it performs, i.e. model validation, compare with process model, lipid simulation at different C/N ratio
* To generate biomass for different C/N concentration
* - When simulate lipid production at different C/N concentration and N-limit, use the function BiomassInDifConditions to generate the corresponding biomass reaction 
