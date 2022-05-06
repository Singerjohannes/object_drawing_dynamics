% function [photo_accs, drawing_accs, sketch_accs) = meg_cat_temporal_generalization(photo_data,drawing_data,sketch_data, n_cat,n_perm,n_average)
%
% This code can be used to compute the temporal generalization of category decoding on the MEG data for the object drawing MEG-fMRI experiment.
% We do a number of things:
% 1. We bring the data in the right format for the classification
% 2. We compute pseudotrials based on the data
% 3. We noise normalize the data
% 4. We train a SVM classifier for every timepoint and test it on all the other timepoints and repeat that for every timepoint
%    - we do this for every type of depiction separately
%    (photo,drawing,sketch)
%
% INPUT:
%   photo/drawing/sketch_data: preprocessed data files in the format
%   obtained from fieldtrip (ft_preprocessing)
%   n_cat : number of categories in the experimental design -> should
%   correspond to the number of unique trigger values in the data files
%   ranging from 1 to n_cat
%   n_perm: number of times the trials should be randomized and
%   pseudotrials should be computed on these shuffled trials
%   n_average: how many trials should be averaged into one pseudotrials
%   (number of repitions divided by n_average must yield a natural number)
%
% OUTPUT:
%   final_photo/drawing/sketch_accs: three matrices with the classification accuracies over all timepoints
%                               averaged across all combinations of categories

function [final_photo_acc, final_drawing_acc, final_sketch_acc] = meg_cat_temporal_generalization(photo_data,drawing_data,sketch_data, n_cat, n_perm, n_average)

% extract data and trial labels for all three conditions separately -
% without condition triggers (>100)
photo_class_data = photo_data.trial(find(photo_data.trialinfo<100));
photo_trialinfo_exp = photo_data.trialinfo(find(photo_data.trialinfo<100));
drawing_class_data = drawing_data.trial(find(drawing_data.trialinfo<100));
drawing_trialinfo_exp = drawing_data.trialinfo(find(drawing_data.trialinfo<100));
sketch_class_data = sketch_data.trial(find(sketch_data.trialinfo<100));
sketch_trialinfo_exp = sketch_data.trialinfo(find(sketch_data.trialinfo<100));

n_trials = 384; %number of trials per block
n_block = 3; % number of blocks for one condition
n_rep = 24; % number of repititions per stimulus
n_suptrial = n_rep/n_average; %number of suptrials
%exception for case where trials are lost
if any([length(photo_trialinfo_exp) length(drawing_trialinfo_exp) length(sketch_trialinfo_exp)] < (n_trials*n_block)),
    n_suptrial=n_suptrial-1;
end


% intialize and fill data matrices
n_channels = size(photo_class_data{1},1);
n_timepoints = size(photo_class_data{1},2);
photo_class_data_mat = NaN(n_cat,n_rep,n_channels,n_timepoints);
drawing_class_data_mat = NaN(n_cat,n_rep,n_channels,n_timepoints);
sketch_class_data_mat = NaN(n_cat,n_rep,n_channels,n_timepoints);

%exception for case where trials are lost
if any([length(photo_trialinfo_exp) length(drawing_trialinfo_exp) length(sketch_trialinfo_exp)] < (n_trials*n_block)),
    photo_class_data_mat = NaN(n_cat,n_suptrial*n_average,n_channels,n_timepoints);
    drawing_class_data_mat = NaN(n_cat,n_suptrial*n_average,n_channels,n_timepoints);
    sketch_class_data_mat = NaN(n_cat,n_suptrial*n_average,n_channels,n_timepoints);
end

% bring data in category x trials x channels x timepoints format

for this_cat = 1:n_cat
    
    %exception for subject where a few trials where lost
    if any([length(photo_trialinfo_exp) length(drawing_trialinfo_exp) length(sketch_trialinfo_exp)] < (n_trials*n_block))
        
        sel_photo_trials = permute(cat(3,photo_class_data{[find(photo_trialinfo_exp == this_cat)]}),[3,1,2]);
        sel_drawing_trials = permute(cat(3,drawing_class_data{[find(drawing_trialinfo_exp == this_cat)]}),[3,1,2]);
        sel_sketch_trials = permute(cat(3,sketch_class_data{[find(sketch_trialinfo_exp == this_cat)]}),[3,1,2]);
        
        photo_class_data_mat(this_cat,:,:,:) = sel_photo_trials(1:n_suptrial*n_average,:,:);
        drawing_class_data_mat(this_cat,:,:,:) = sel_drawing_trials(1:n_suptrial*n_average,:,:);
        sketch_class_data_mat(this_cat,:,:,:) = sel_sketch_trials(1:n_suptrial*n_average,:,:);
    else
        photo_class_data_mat(this_cat,:,:,:) = permute(cat(3,photo_class_data{[find(photo_trialinfo_exp == this_cat)]}),[3,1,2]);
        drawing_class_data_mat(this_cat,:,:,:) = permute(cat(3,drawing_class_data{[find(drawing_trialinfo_exp == this_cat)]}),[3,1,2]);
        sketch_class_data_mat(this_cat,:,:,:) = permute(cat(3,sketch_class_data{[find(sketch_trialinfo_exp == this_cat)]}),[3,1,2]);
    end
end

% clear some memory
clear photo_class_data drawing_class_data sketch_class_data

% start outer loop with computing supertrials

%initialize results

final_photo_acc = nan(n_perm,n_timepoints,n_timepoints);
final_drawing_acc = nan(n_perm,n_timepoints,n_timepoints);
final_sketch_acc = nan(n_perm,n_timepoints,n_timepoints);

for perm =1:n_perm
    
    fprintf('Starting permutation number %i', perm)
    % randomly re-order trials and average trials
    photo_class_data_pseudotrl = create_pseudotrials(photo_class_data_mat,n_average);
    drawing_class_data_pseudotrl = create_pseudotrials(drawing_class_data_mat,n_average);
    sketch_class_data_pseudotrl = create_pseudotrials(sketch_class_data_mat,n_average);
    
    % noise normalize the data
    whitened_photo_data = mvnn_whitening(photo_class_data_pseudotrl);
    whitened_drawing_data = mvnn_whitening(drawing_class_data_pseudotrl);
    whitened_sketch_data = mvnn_whitening(sketch_class_data_pseudotrl);
    
    
    % reshape data for decoding
    
    whitened_photo_data = reshape(whitened_photo_data, [], n_channels,n_timepoints);
    whitened_drawing_data = reshape(whitened_drawing_data, [], n_channels,n_timepoints);
    whitened_sketch_data = reshape(whitened_sketch_data, [], n_channels,n_timepoints);
    
    % update the trialinfo for the supertrials
    
    trialinfo_exp = repmat([1:n_cat],1,n_suptrial);
    
    % do classification using SVM
    
    %labels for training and testing only need to be specified once
    label_train = repmat([1:48],1,n_suptrial-1)';
    label_test = [1:48];
    
    % create crossvalidation design
    
    design = ones(n_suptrial,n_suptrial);
    
    for set = 1:n_suptrial
        
        design(set,set) = 0;
    end
    
    % classification loop
    
    for time_i=n_timepoints:-1:1
        
        for iter = 1:n_suptrial
            % get train and test indices for the current categories
            trainind = [];
            testind = [];
            for i=1:n_suptrial
                if design(iter,i) == 1
                    trainind = [1+n_cat*(i-1):n_cat*i trainind];
                elseif design(iter,i) ==0
                    testind = [1+n_cat*(i-1):n_cat*i testind];
                end
            end
            
            %train model for every depiction
            model_photo = svmtrain(label_train, whitened_photo_data(trainind,:,time_i),'-s 0 -t 0 -c 1 -b 0 -q'); %#ok<SVMTRAIN>
            model_drawing = svmtrain(label_train, whitened_drawing_data(trainind,:,time_i),'-s 0 -t 0 -c 1 -b 0 -q');
            model_sketch = svmtrain(label_train, whitened_sketch_data(trainind,:,time_i) ,'-s 0 -t 0 -c 1 -b 0 -q');
            
            for time_j = n_timepoints:-1:1
                
                % get prediction for all depictions and all timepoints
                [~,~,photo_decision_values] = svmpredict(label_test',whitened_photo_data(testind,:,time_j),model_photo,'-q');
                [~,~,drawing_decision_values] = svmpredict(label_test',whitened_drawing_data(testind,:,time_j),model_drawing,'-q');
                [~,~,sketch_decision_values] = svmpredict(label_test',whitened_sketch_data(testind,:,time_j),model_sketch,'-q');
                
                %assign ddecoding output to condition specific structures
                photo_decoding_out(iter,time_j).true_labels = label_test;
                photo_decoding_out(iter,time_j).model = model_photo;
                photo_decoding_out(iter,time_j).decision_values = photo_decision_values;
                
                drawing_decoding_out(iter,time_j).true_labels = label_test;
                drawing_decoding_out(iter,time_j).model = model_drawing;
                drawing_decoding_out(iter,time_j).decision_values = drawing_decision_values;
                
                sketch_decoding_out(iter,time_j).true_labels = label_test;
                sketch_decoding_out(iter,time_j).model = model_sketch;
                sketch_decoding_out(iter,time_j).decision_values = sketch_decision_values;
            end
        end
        for time_j =n_timepoints:-1:1
            final_photo_acc(perm,time_i,time_j) = transres_accuracy_pairwise(photo_decoding_out(:,time_j), 0.5);
            final_drawing_acc(perm,time_i,time_j) = transres_accuracy_pairwise(drawing_decoding_out(:,time_j), 0.5);
            final_sketch_acc(perm,time_i,time_j) = transres_accuracy_pairwise(sketch_decoding_out(:,time_j), 0.5);
        end
    end
    if mod(time_i,10)==0
        fprintf('.');
    end
end
fprintf('\nEnd permutation number: %i',perm);
%average over permutations
final_photo_acc = squeeze(mean(final_photo_acc,1));
final_drawing_acc = squeeze(mean(final_drawing_acc,1));
final_sketch_acc = squeeze(mean(final_sketch_acc,1));

% add the matrix and the transpose of the matrix and take the mean ->
% because symmetry is assumed
%final_photo_acc = (final_photo_acc+final_photo_acc')/2;
%final_drawing_acc = (final_drawing_acc+final_drawing_acc')/2;
%final_sketch_acc = (final_sketch_acc+final_sketch_acc')/2;
fprintf('\ndone!')
end
