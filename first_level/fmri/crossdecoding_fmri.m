% function decoding_sl(condition_names,labels,beta_dir,out_dir,cfg)
%
% Wrapper script for decoding using a searchlight.
%
% Input variables:
%   condition_names: Names of all regressors to be used for classification
%   labels: Labels that should be paired with condition_names (e.g. [-1 1])
%   beta_dir: name where results are which are used for classification
%   out_dir: name of folder where results are saved
%   cfg (optional): config variable used for the decoding

function results = crossdecoding_fmri(condition_names,labels,beta_dir,out_dir,cfg)


if ~exist('cfg','var')
    cfg = decoding_defaults;
else
    cfg = decoding_defaults(cfg);
end

% Get subject number
try
    cfg.sn = str2double(regexp(beta_dir,'(?<=sub).[0-9]','match'));
catch
    warning('Could not extract subject number.')
end

% Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
cfg.results.dir = out_dir;
cfg.results.overwrite = 1;

% Set the filepath where your SPM.mat and all related betas are, e.g. 'c:\exp\glm\model_button'
% done already

% Set the filename of your brain mask (or your ROI masks as cell matrix) 
% for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR 
% for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
try cfg.files.mask;
catch
    cfg.files.mask = fullfile(beta_dir,'mask.img');
    if ~exist(cfg.files.mask,'file')
        cfg.files.mask = fullfile(beta_dir,'mask.nii');
        if ~exist(cfg.files.mask,'file')
            error('Mask not found in %s',cfg.files.mask)
        end
    end
end


if cfg.noisenorm == 1 
    
% These parameters enable the multivariate noise normalization using the
% residuals
cfg.scale.method = 'cov'; % we scale by noise covariance
cfg.scale.estimation = 'separate'; % we scale all data for each run separately while iterating across searchlight spheres
cfg.scale.shrinkage = 'lw2'; % Ledoit-Wolf shrinkage retaining variances

if ~cfg.hrf_fitting
    
[misc.residuals,cfg.files.residuals.chunk] = residuals_from_spm(fullfile(beta_dir,'SPM.mat'),cfg.files.mask); % this only needs to be run once and can be saved and loaded 

elseif cfg.hrf_fitting 
        
    res_names = dir(fullfile(beta_dir, '*Res_*')); 
    
    res_names = {res_names.name}';
    
    misc.residuals = residuals_without_spm(fullfile(beta_dir,res_names),cfg.files.mask); 
    
    cfg.files.residuals.chunk = repelem([1:12], 251); 
    
end 
end

cfg.searchlight.unit = 'voxels'; %mm
cfg.searchlight.radius = 4; %4 this will yield a searchlight radius of 4 voxels
cfg.searchlight.spherical = 1;
cfg.plot_design = 0;
cfg.verbose = 2; % you want all information to be printed on screen

%% Nothing needs to be changed below for a standard leave-one-run out cross
%% validation analysis.

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
regressor_names = design_from_spm(beta_dir);

try
    cfg = decoding_describe_data(cfg,condition_names,labels,regressor_names,beta_dir);
catch
    betas = dir(fullfile(beta_dir,'*wbeta*'));
    betas = {betas.name}';

    % Extract all information for the cfg.files structure (labels will be [1 -1] )
    cfg = decoding_describe_data(cfg,condition_names,labels,regressor_names,fullfile(beta_dir,betas));
end 

cfg.files.chunk = repmat([1 2 3 4],1 ,384/4)'; % little hack to pretend that all betas come from the same chunks (creates the right design) 

cfg.files.twoway = 1; % set this to 1 so that both train-test directions are evaluated 

cfg.files.xclass = [ones(1,length(cfg.files.chunk)/2) ones(1,length(cfg.files.chunk)/2)*2]; 
cfg.design = make_design_xclass_cv(cfg); % This creates the leave-one-run-out cross validation design.

if ~cfg.parallel % check if decoding should be performed in a parallelized way or not
    % Run decoding
    if cfg.noisenorm == 1
        results = decoding(cfg,[],misc);
    else
        results = decoding(cfg);
    end
    
elseif cfg.parallel % if decoding should be parallelized 
    
    % get the number of searchlights for parallelization
    mask = spm_read_vols(spm_vol(cfg.files.mask{1}));
    num_searchlights = size(mask(mask==1),1);
    % start the parallel pool
    pool =  parcluster('local');
    % get number of workers (cpus)
    num_workers = pool.NumWorkers/2;
    searchlights_per_worker = ceil(num_searchlights/num_workers); % Divide the task up into the number of workers
    fprintf('Starting searchlight decoding on %i parallel cores with %i searchlights for each worker', num_workers, searchlights_per_worker)
    parfor crun = 1:num_workers
        results{crun} = decoding_parallel_wrapper(cfg,searchlights_per_worker,crun,misc)
    end
    all_results = results{1};
    for crun = 2:num_workers
        all_results.decoding_subindex = [all_results.decoding_subindex; results{crun}.decoding_subindex];
        all_results.accuracy_pairwise.output(results{crun}.decoding_subindex) = results{crun}.accuracy_pairwise.output(results{crun}.decoding_subindex);
    end
    results = all_results;
    disp('Decoding on the whole brain complete, saving results')
    combine_write_decoding_results(beta_dir,results, fullfile(cfg.results.dir,'res_accuracy_minus_chance.nii'))
    save(fullfile(cfg.results.dir,'res_accuracy_minus_chance.mat'),'results')
end
end       