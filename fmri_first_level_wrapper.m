%% Script for running all the first level analyses of the fMRI data
%
% possible analysis steps are:
% Category Decoding - ROI and searchlight
% Category Crossdecoding - ROI and searchlight
% RSA - ROI

clear all
clc

%setup paths

path = pwd;
figure_path = fullfile(path,'figures');

% add utils

addpath(fullfile(path,'utils'));

% add first level functions

addpath(fullfile(path,'first_level','fmri'));

% add the decoding toolbox

addpath(fullfile(path,'tdt_3.999','decoding_toolbox'));

%initialize decoding toolbox
decoding_defaults;

% set plot defaults

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3, 'DefaultAxesFontName', 'Helvetica','DefaultTextFontName', 'Helvetica')

% get colormap
cmap = colormap('redblueTecplot');
close all

% specify which steps to compute

cfg.do.decoding_ROI = 1; % runs the category decoding for a single subject for photos, drawings and sketches separately in 2 ROIs: EVC and LOC
cfg.do.decoding_searchlight = 0; % runs the category decoding for a single subject for photos, drawings and sketches separately across the whole brain
cfg.do.cross_decoding_ROI = 0; % runs the category crossdecoding for a single subject for three comparisons: photo-drawing, drawing-sketch and photo-sketch separately in 2 ROIs: EVC and LOC
cfg.do.cross_decoding_searchlight = 0; % runs the category crossdecoding for a single subject for three comparisons: photo-drawing, drawing-sketch and photo-sketch separately across the whole brain
cfg.do.RSA = 0; % computes representational dissimilarity matrices (RDMs) for a single subject for two ROIs: EVC and LOC and for photos, drawings and sketches separately


%% category decoding - ROI

if cfg.do.decoding_ROI
    
    cfg.analysis = 'roi';
    cfg.noisenorm = 1; % specify if multivariate noise normalization should be applied
    cfg.hrf_fitting = 1; % specify if hrf fitting betas should be used for loading the residuals
    cfg.parallel = 0; % specify if decoding should be parallelized
    conds = {'Photo'; 'Drawing'; 'Sketch'};
    for this_cond_ind = 1:length(conds) % specifiy condition names "Photo/Drawing/Skech_[1-48]"
        condition_names = cell(1,48);
        cond = conds{this_cond_ind};
        for i=1:48
            condition_names(i) = {[cond, '_', num2str(i)]};
        end
        labels = [1:48]; % set the labels for the decoding
        beta_dir = fullfile(path,'data','fmri','preproc','sub15_betas','fitted'); % path to the single subject betas in the individual subject space
        out_dir = fullfile(path,'data','fmri','decoding','single_sub',cfg.analysis,cond);  % path where results should be saved
        roi_dir = fullfile(path,'data','fmri','preproc','masks'); % folder where ROI and searchlight masks are stored
        cfg.design.function.name = 'make_design_cv'; % function to create the crossvalidation design
        cfg.results.output = {'accuracy_pairwise'}; % desired output of the decoding - here accuracy pairwise
        cfg.files.mask = {fullfile(roi_dir, 'evcmask.nii');fullfile(roi_dir, 'loc_mask.nii')}; % path to the ROI masks
        
        decoding_fmri(condition_names,labels,beta_dir,out_dir,cfg);
    end
end

%% category decoding - searchlight


if cfg.do.decoding_searchlight
    
    cfg.analysis = 'searchlight';
    cfg.noisenorm = 1; % specify if multivariate noise normalization should be applied
    cfg.hrf_fitting = 1; % specify if hrf fitting betas should be used for loading the residuals
    cfg.parallel = 1; % specify if decoding should be parallelized
    conds = {'Photo'; 'Drawing'; 'Sketch'};
    for this_cond_ind = 1:length(conds) % specifiy condition names "Photo/Drawing/Skech_[1-48]"
        condition_names = cell(1,48);
        cond = conds{this_cond_ind};
        for i=1:48
            condition_names(i) = {[cond, '_', num2str(i)]};
        end
        labels = [1:48]; % set the labels for the decoding
        beta_dir = fullfile(path,'data','fmri','preproc','sub15_betas','normalized'); % path to the single subject betas in the MNI space
        out_dir = fullfile(path,'data','fmri','decoding','single_sub',cfg.analysis,cond); % path where results should be saved
        roi_dir = fullfile(path,'data','fmri','preproc','masks'); % folder where ROI and searchlight masks are stored
        cfg.design.function.name = 'make_design_cv'; % function to create the crossvalidation design
        cfg.results.output = {'accuracy_pairwise'}; % desired output of the decoding - here accuracy pairwise
        cfg.files.mask = {fullfile(roi_dir, 'searchlight_mask.nii')}; % path to the searchlight mask
        
        decoding_fmri(condition_names,labels,beta_dir,out_dir,cfg);
    end
end

%% category cross-decoding - ROI

if cfg.do.cross_decoding_ROI
    
    cfg.analysis = 'roi';
    cfg.noisenorm = 1; % specify if multivariate noise normalization should be applied
    cfg.hrf_fitting = 1; % specify if hrf fitting betas should be used for loading the residuals
    cfg.parallel = 0; % specify if decoding should be parallelized
    conds = {'Photo'; 'Drawing'; 'Sketch'};
    for this_cond_ind = 1:length(conds)-1 % specify condition names for training and testing e.g. "Photo_[1-48]” and "Drawing_[1-48]"
        for that_cond_ind = this_cond_ind+1:length(conds)
            condition_names = cell(1,48*2);
            this_cond = conds{this_cond_ind};
            for i=1:48
                condition_names(i) = {[this_cond, '_', num2str(i)]};
            end
            that_cond = conds{that_cond_ind};
            for i=1+48:48*2
                condition_names(i) = {[that_cond, '_', num2str(i-48)]};
            end
            labels = [1:48 1:48]; % set the labels for the decoding
            beta_dir = fullfile(path,'data','fmri','preproc','sub15_betas','fitted'); % path to the single subject betas in the individual subject space
            out_dir = fullfile(path,'data','fmri','crossdecoding','single_sub',cfg.analysis,[this_cond,'_',that_cond]);  % path where results should be saved
            roi_dir = fullfile(path,'data','fmri','preproc','masks'); % folder where ROI and searchlight masks are stored
            cfg.results.output = {'accuracy_pairwise'}; % desired output of the decoding - here accuracy pairwise
            cfg.files.mask = {fullfile(roi_dir, 'evcmask.nii');fullfile(roi_dir, 'loc_mask.nii')}; % path to the ROI masks
            
            crossdecoding_fmri(condition_names,labels,beta_dir,out_dir,cfg);
        end
    end
end

%% category cross-decoding - searchlight


if cfg.do.cross_decoding_searchlight
    
    cfg.analysis = 'searchlight';
    cfg.noisenorm = 1; % specify if multivariate noise normalization should be applied
    cfg.hrf_fitting = 1; % specify if hrf fitting betas should be used for loading the residuals
    cfg.parallel = 1; % specify if decoding should be parallelized
    conds = {'Photo'; 'Drawing'; 'Sketch'};
    for this_cond_ind = 1:length(conds)-1 % specify condition names for training and testing e.g. "Photo_[1-48]” and "Drawing_[1-48]"
        for that_cond_ind = this_cond_ind+1:length(conds)
            condition_names = cell(1,48*2);
            this_cond = conds{this_cond_ind};
            for i=1:48
                condition_names(i) = {[this_cond, '_', num2str(i)]};
            end
            that_cond = conds{that_cond_ind};
            for i=1+48:48*2
                condition_names(i) = {[that_cond, '_', num2str(i-48)]};
            end
            labels = [1:48 1:48]; % set the labels for the decoding
            beta_dir = fullfile(path,'data','fmri','preproc','sub15_betas','normalized'); % path to the single subject betas in the individual subject space
            out_dir = fullfile(path,'data','fmri','crossdecoding','single_sub',cfg.analysis,[this_cond,'_',that_cond]);  % path where results should be saved
            roi_dir = fullfile(path,'data','fmri','preproc','masks'); % folder where ROI and searchlight masks are stored
            cfg.results.output = {'accuracy_pairwise'}; % desired output of the decoding - here accuracy pairwise
            cfg.files.mask = {fullfile(roi_dir, 'searchlight_mask.nii')}; % path to the searchlight mask
            
            crossdecoding_fmri(condition_names,labels,beta_dir,out_dir,cfg);
        end
    end
end

%% RSA - ROI

if cfg.do.RSA
    
    clear cfg
    cfg.analysis = 'roi';
    cfg.hrf_fitting = 1;
    cfg.noisenorm = 1; 
    conds = {'Photo'; 'Drawing'; 'Sketch'};
    for this_cond_ind = 1:length(conds)
        condition_names = cell(1,48);
        cond = conds{this_cond_ind};
        for i=1:48
            condition_names(i) = {[cond, '_', num2str(i)]};
        end
        labels = [1:48];
        beta_dir = fullfile(path,'data','fmri','preproc','sub15_betas','fitted'); % path to the single subject betas in the individual subject space
        out_dir = fullfile(path,'data','fmri','rsa','single_sub',cfg.analysis,cond);  % path where results should be saved
        roi_dir = fullfile(path,'data','fmri','preproc','masks'); % folder where ROI and searchlight masks are stored
        cfg.files.mask = {fullfile(roi_dir, 'evcmask.nii');fullfile(roi_dir, 'loc_mask.nii')}; % path to the ROI masks
        decoding_similarity_pearson(condition_names,labels,beta_dir,out_dir,cfg);
    end
end