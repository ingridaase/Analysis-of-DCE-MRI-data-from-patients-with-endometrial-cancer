# Master's thesis repository
This repository includes the code for the master's thesis "Analysis of DCE-MRI data from patients with endometrial cancer". In this project the pharmacokinetic model extended Tofts' method (ETM) was used to analyse DCE-MRI data from patients with endometrial cancer for quantitative perfusion parameters. 

Most of the algorithms used for running ETM on DCE-MRI data are included in this repository. Two images are required for running ETM, including whole-volume tumor masks and the DCE-MRI images. In addition, an arterial input function is required. Examples of how the resulting parameters were analysed and visualized are also included in this repository. 

## Contents of the repository 
This repository is organized in three folders: 

**The arterial input function:** This folder contains three different methods for estimating an arterial input function. An AIF ROI can be annotated from a DCE-MRI image, using the "annotate_AIF.py" file. A population-based AIF can be calculated based on manually annotated AIFs using the "population_AIF.py". An automatic AIF can be estimated from a DCE-MRI image using the "AIF_deterministic.py". An additional jupyter notebook is added to the folder to show the different steps of the automatic AIF algorithm. 

**Extended Tofts method:** This folder contains algorithms for running ETM. An arterial input function is required for ETM. This folder contains code for running ETM inside the research information system, where the main file for running ETM is "stub.py". This folder also contains code for running ETM on a local machine, where "runtofts_script.py" is the main file for running ETM. In this file comments are included to illustrate which changes are required when different AIF methods are used.

**Visualization and statistics:** This folder contains examples of how the resulting model parameters were analysed and visualized. Examples of how the model parameters obatined from local modeling were analysed are included in the "local_modeling_results.ipynb" notebook. Examples of how the model paramters from the research environment were analysed and how clinical parameters were included in the analysis are shown in the "research_environment_results.ipynb" notebook. Additional, the notebook "parameter_maps.ipynb" shows how the parameters maps were created from the voxelwise obtained model parameters.

## Requirements 
The requirements for the algorithms in this project are included in the "environment.yml" file. 

## Contact
I can be contacted at ingrid-aase.cowles@outlook.com for quiestions regarding this project. 
