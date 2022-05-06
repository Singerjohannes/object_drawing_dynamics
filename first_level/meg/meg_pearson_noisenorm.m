function [photo_RDM, drawing_RDM, sketch_RDM, big_RDM] = meg_pearson_noisenorm(photo_data,drawing_data,sketch_data, n_cat)
%
% This code can be used to compute the pearson distance for the MEG data after noise normalization for the object drawing MEG-fMRI experiment.
% We do a number of things:
% 1. We read the data
% 2. We apply noise normalization on the data
% 3. We average the trial data and compute the correlation distance for every timepoint and every
%    category combination - importantly we do this for every type of
%    depiction (photo,drawing,sketch) separately
%
% INPUT:
%   photo/drawing/sketch_data: preprocessed data files in the format
%   obtained from fieldtrip (ft_preprocessing)
%   n_cat : number of categories in the experimental design -> should
%   correspond to the number of unique trigger values in the data files
%   ranging from 1 to n_cat
%
% OUTPUT:
%   photo/drawing/sketch_RDMs: three matrices with the cross-validated distances over all timepoints

% extract data and trial labels for all three conditions separately -
% without condition triggers (all trials with triggers <100)
photo_class_data = photo_data.trial(find(photo_data.trialinfo<100));
photo_trialinfo_exp = photo_data.trialinfo(find(photo_data.trialinfo<100));
drawing_class_data = drawing_data.trial(find(drawing_data.trialinfo<100));
drawing_trialinfo_exp = drawing_data.trialinfo(find(drawing_data.trialinfo<100));
sketch_class_data = sketch_data.trial(find(sketch_data.trialinfo<100));
sketch_trialinfo_exp = sketch_data.trialinfo(find(sketch_data.trialinfo<100));

n_trials = 384; %number of trials per block
n_block = 3; % number of blocks for one condition
n_rep = 24; % number of repititions per stimulus
n_suptrial = n_rep/n_block; %number of suptrials
%exception for case where trials are lost
if any([length(photo_trialinfo_exp) length(drawing_trialinfo_exp) length(sketch_trialinfo_exp)] < (n_trials*n_block)),
    n_suptrial=7;
end

% intialize and fill data matrices
n_channels = size(photo_class_data{1},1);
n_timepoints = size(photo_class_data{1},2);
photo_class_data_mat = NaN(n_cat,n_rep,n_channels,n_timepoints);
drawing_class_data_mat = NaN(n_cat,n_rep,n_channels,n_timepoints);
sketch_class_data_mat = NaN(n_cat,n_rep,n_channels,n_timepoints);

%exception for case where trials are lost
if any([length(photo_trialinfo_exp) length(drawing_trialinfo_exp) length(sketch_trialinfo_exp)] < (n_trials*n_block)),
    photo_class_data_mat = NaN(n_cat,n_suptrial*n_block,n_channels,n_timepoints);
    drawing_class_data_mat = NaN(n_cat,n_suptrial*n_block,n_channels,n_timepoints);
    sketch_class_data_mat = NaN(n_cat,n_suptrial*n_block,n_channels,n_timepoints);
end

% bring data in category x trials x channels x timepoints format

for this_cat = 1:n_cat
    
    %exception for subject where a few trials where lost
    if any([length(photo_trialinfo_exp) length(drawing_trialinfo_exp) length(sketch_trialinfo_exp)] < (n_trials*n_block))
        
        sel_photo_trials = permute(cat(3,photo_class_data{[find(photo_trialinfo_exp == this_cat)]}),[3,1,2]);
        sel_drawing_trials = permute(cat(3,drawing_class_data{[find(drawing_trialinfo_exp == this_cat)]}),[3,1,2]);
        sel_sketch_trials = permute(cat(3,sketch_class_data{[find(sketch_trialinfo_exp == this_cat)]}),[3,1,2]);
        
        photo_class_data_mat(this_cat,:,:,:) = sel_photo_trials(1:n_suptrial*n_block,:,:);
        drawing_class_data_mat(this_cat,:,:,:) = sel_drawing_trials(1:n_suptrial*n_block,:,:);
        sketch_class_data_mat(this_cat,:,:,:) = sel_sketch_trials(1:n_suptrial*n_block,:,:);
    else
        photo_class_data_mat(this_cat,:,:,:) = permute(cat(3,photo_class_data{[find(photo_trialinfo_exp == this_cat)]}),[3,1,2]);
        drawing_class_data_mat(this_cat,:,:,:) = permute(cat(3,drawing_class_data{[find(drawing_trialinfo_exp == this_cat)]}),[3,1,2]);
        sketch_class_data_mat(this_cat,:,:,:) = permute(cat(3,sketch_class_data{[find(sketch_trialinfo_exp == this_cat)]}),[3,1,2]);
    end
end

% concatenate the data for BIG RDM
all_class_data_mat = cat(1,photo_class_data_mat, drawing_class_data_mat, sketch_class_data_mat);

% start outer loop with computing supertrials

n_permutations = 1; % how many permutations for the selection of supertrials

%initialize results

photo_RDM = nan(n_permutations,n_cat,n_cat,n_timepoints);
drawing_RDM = nan(n_permutations,n_cat,n_cat,n_timepoints);
sketch_RDM = nan(n_permutations,n_cat,n_cat,n_timepoints);
big_RDM = nan(n_permutations,n_cat*3,n_cat*3,n_timepoints);

for perm =1:n_permutations
    fprintf('Starting permutation number %i', perm)
    % noise normalize the data
    
    whitened_photo_data = mvnn_whitening(photo_class_data_mat);
    whitened_drawing_data = mvnn_whitening(drawing_class_data_mat);
    whitened_sketch_data = mvnn_whitening(sketch_class_data_mat);
    whitened_all_data = mvnn_whitening(all_class_data_mat);
    
    % loop over time points
    
    for time=size(photo_class_data_mat,4):-1:1
        
        
        % compute average over folds
        photo_RDM(perm,:,:,time) = 1-corr(squeeze(mean(whitened_photo_data(:,:,:,time),2))');
        drawing_RDM(perm,:,:,time ) = 1-corr(squeeze(mean(whitened_drawing_data(:,:,:,time),2))');
        sketch_RDM(perm,:,:,time) = 1-corr(squeeze(mean(whitened_sketch_data(:,:,:,time),2))');
        big_RDM(perm,:,:,time) = 1-corr(squeeze(mean(whitened_all_data(:,:,:,time),2))');
        
        if mod(time,10)==0
            fprintf('.');
        end
    end
    fprintf('\nEnd permutation number: %i',perm);
end
photo_RDM = squeeze(nanmean(photo_RDM,1));
drawing_RDM = squeeze(nanmean(drawing_RDM,1));
sketch_RDM = squeeze(nanmean(sketch_RDM,1));
big_RDM = squeeze(nanmean(big_RDM,1));
fprintf('\ndone!')
end