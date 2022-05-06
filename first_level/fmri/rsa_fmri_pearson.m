function results = rsa_fmri_pearson(condition_names,labels,beta_dir,out_dir,cfg)

if ~exist('cfg','var')
    cfg = decoding_defaults;
else
    cfg = decoding_defaults(cfg);
end

% Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
cfg.results.dir = out_dir;
cfg.results.overwrite = 1;

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

% set everything to calculate (dis)similarity estimates
cfg.decoding.software = 'distance'; % the difference to 'similarity' is that this averages across data with the same label
cfg.decoding.method = 'classification'; % this is more a placeholder
cfg.decoding.train.classification.model_parameters = 'pearson'; % cross-validated Euclidean after noise normalization

% This option below averages across (dis)similarity matrices of each
% cross-validation iteration and across all cells of the lower diagonal
% (i.e. all distance comparisons). If you want the entire matrix, consider
% using 'other_average' which only averages across cross-validation
% iterations. Alternatively, you could use the output 'RSA_beta' which is
% more general purpose, but a little more complex.
cfg.results.output = 'other_average_RDV';

% These parameters carry out the multivariate noise normalization using the
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

%% Nothing needs to be changed below for standard dissimilarity estimates using all data

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
regressor_names = design_from_spm(beta_dir);

% Extract all information for the cfg.files structure (labels will be [1 -1] )
cfg = decoding_describe_data(cfg,condition_names,labels,regressor_names,beta_dir);

% This creates a design without cross-validation 
cfg.design = make_design_similarity(cfg);

% Run decoding
results = decoding(cfg,[],misc);
end 