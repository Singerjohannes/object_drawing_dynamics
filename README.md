# object_drawing_dynamics
This repository containes code for the paper "The spatiotemporal neural dynamics of object recognition for natural images and line drawings". With the material contained in this repository all of the group level results and plots in the paper can be reproduced. All the first level results needed for the group level analyses are provided on OSF (Link: ). In addition, code for the first-level analyses is provided along with preprocessed data from sample subjects, which can be retrieved from OSF (Link: ). Link to preprint: 


## Requirements: 

To run the code in this repository you will need the following toolboxes: 

- SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- The Decoding Toolbox 3.999 or newer (https://sites.google.com/site/tdtdecodingtoolbox/) 


## Group-level analyses: 

All group-level results and the statistics in the paper can be reproduced with the code provided in this repository.  
For this, you first need to download the data from OSF (Link:), unzip the data and put the folders in the following folder structure: 

data/  
data/fmri/  
data/fmri/decoding/  
data/fmri/crossdecoding/  
data/fmri/preproc/  
data/fmri/rsa/  
data/meg/decoding/  
data/meg/crossdecoding/  
data/meg/preproc/ data/meg/rsa/  
data/meg/temporal_generalization/

With this folder structure all scripts can be executed without changes.  
To reproduce the results for the MEG-data execute the script meg_decoding_group_wrapper.m  
To reproduce the results for the fMRI-ROI-data execute the script fmri_group_roi_wrapper.m  
To reproduce the results for the fMRI-volume-data execute the script fmri_group_searchlight_wrapper.m  
To reproduce the MEG-fMRI fusion results execute the script MEG_fmri_fusion_wrapper.m 

## First-level analyses:

We provide exemplary preprocessed MEG and fMRI data for a single subject to demonstrate how the first-level results are computed.  
Make sure to unzip the files in the /data/meg/preproc/ and data/fmri/preproc/ folders before running the first-level scripts.  

To run the decoding, cross-decoding, temporal generalization and representational similarity analyses for the MEG-data run meg_first_level_wrapper.m and specify in the script which part (e.g. decoding, cross-decoding etc.) you want to run. 
Depending on the type of analyses this might be time intensive (up to 15-20 hours for the temporal generalization analysis). 

To run the decoding (ROI or searchlight), cross-decoding (ROI or searchlight) and representational similarity analyses (only ROI) for the fMRI-data run fmri_first_level_wrapper.m and specify in the script which part (e.g. decoding, cross-decoding etc.) you want to run. 
Depending on the type of analyses this might be time and/or memory intensive (up to 3-4 hours for the searchlight analyses; with parallelization and the amount of cores used the memory needed increases). 


