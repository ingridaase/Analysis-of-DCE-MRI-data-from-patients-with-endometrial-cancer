# Master's thesis repository
This repository includes the code for the master's thesis "Analysis of DCE-MRI data from patients with endometrial cancer". In this project the pharmacokinetic model extended Tofts' method (ETM) was used to analyse DCE-MRI data from patients with endometrial cancer for quantitative perfusion parameters. 

Most of the algorithms used for running ETM on DCE-MRI data are included in this repository. Two images are required for running ETM, including whole-volume tumor masks and the DCE-MRI images. In addition, an arterial input function is required. Examples of how the resulting parameters were analysed and visualized are also included in this repository. 

## Contents of the repository 
This repository is organized in three folders: 

**The arterial input function:** This folder contains three different methods for estimating an arterial input function, including manual annotation of the AIF, an atuomatic algorithm for estimating the AIF and a population-based AIF. The population-based AIF is based on manually annotated AIFs, whereas the two other methodologies requires an input DCE-MRI image in DICOM format. 

**Extended Tofts method:** This folder contains algorithms for running ETM. An arterial input function is required for ETM. This folder contains scripts for running ETM inside the research information system, in addition to scripts for running ETM on a local machine.

**Visualization and statistics:** This folder contains examples of how the resulting model parameters were analysed and visualized. Additionally, clinical parameters were imported and included in the analysis. 

## Requirements 
The requirements for the algorithms in this project are included in the "environment.yml" file. 

## Contact
I can be contacted at ingrid-aase.cowles@outlook.com for quiestions regarding this project. 
