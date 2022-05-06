% function [photo_accs, drawing_accs, sketch_accs] = meg_cat_crossdecoding_suptrials_multiclass(photo_data,drawing_data,sketch_data, n_cat, n_perm,n_average)
%
% This code can be used to compute the category crossdecoding on the MEG data for the object drawing MEG-fMRI experiment.
% We do a number of things:
% 1. We bring the data in the right format for the classification
% 2. We compute pseudotrials based on the data
% 3. We noise normalize the data
% 4. We train and test a SVM classifier for every timepoint finally take the average of all cv-folds and permutations for
%    every timepoint - importantly we do this for possible pair of types of depiction (photo-drawing,drawing-photo,photo-sketch etc.) separately
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
% OUTPUT:
%   final_photo_drawing/drawing_photo/photo_sketch/
%   sketch_photo/drawing_sketch/sketch_drawing_accs: six vectors (one for each comparison between two types of depiction) with the cross-classification accuracies over all timepoints
%                                                    averaged across all combinations of categories

function [final_photo_drawing_acc, final_drawing_photo_acc, final_photo_sketch_acc, final_sketch_photo_acc, final_drawing_sketch_acc, final_sketch_drawing_acc] = meg_cat_crossdecoding(photo_data,drawing_data,sketch_data,n_cat, n_perm, n_average)

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
n_suptrial = n_rep/n_average; % number of pseudotrials for each category

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


% start outer loop with computing supertrials

%initialize results

final_photo_drawing_acc = nan(n_perm,n_timepoints);
final_drawing_photo_acc = nan(n_perm,n_timepoints);
final_photo_sketch_acc = nan(n_perm,n_timepoints);
final_sketch_photo_acc = nan(n_perm,n_timepoints);
final_drawing_sketch_acc = nan(n_perm,n_timepoints);
final_sketch_drawing_acc = nan(n_perm,n_timepoints);

for perm =1:n_perm
    
    % randomly re-order trials and average trials into supertrials
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
    
    for time=n_timepoints:-1:1
        
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
            
            % select training data for all depictions
            photo_data_train = whitened_photo_data(trainind,:,time);
            drawing_data_train = whitened_drawing_data(trainind,:,time);
            sketch_data_train = whitened_sketch_data(trainind,:,time);
            
            % select test data for all depictions
            photo_data_test = whitened_photo_data(testind,:,time);
            drawing_data_test = whitened_drawing_data(testind,:,time);
            sketch_data_test = whitened_sketch_data(testind,:,time);
            
            % train model for every depiction
            model_photo = svmtrain(label_train, photo_data_train,'-s 0 -t 0 -c 1 -b 0 -q'); %#ok<SVMTRAIN>
            model_drawing = svmtrain(label_train, drawing_data_train,'-s 0 -t 0 -c 1 -b 0 -q');
            model_sketch = svmtrain(label_train, sketch_data_train ,'-s 0 -t 0 -c 1 -b 0 -q');
            
            % get prediction for all depictions
            [~,~,photo_drawing_decision_values] = svmpredict(label_test',drawing_data_test,model_photo,'-q');
            [~,~,drawing_photo_decision_values] = svmpredict(label_test',photo_data_test,model_drawing,'-q');
            [~,~,photo_sketch_decision_values] = svmpredict(label_test',sketch_data_test,model_photo,'-q');
            [~,~,sketch_photo_decision_values] = svmpredict(label_test',photo_data_test,model_sketch,'-q');
            [~,~,drawing_sketch_decision_values] = svmpredict(label_test',sketch_data_test,model_drawing,'-q');
            [~,~,sketch_drawing_decision_values] = svmpredict(label_test',drawing_data_test,model_sketch,'-q');
            
            %assign crossdecoding output to condition specific structures
            photo_drawing_out(iter).true_labels = label_test;
            photo_drawing_out(iter).model = model_photo;
            photo_drawing_out(iter).decision_values = photo_drawing_decision_values;
            
            drawing_photo_out(iter).true_labels = label_test;
            drawing_photo_out(iter).model = model_drawing;
            drawing_photo_out(iter).decision_values = drawing_photo_decision_values;
            
            photo_sketch_out(iter).true_labels = label_test;
            photo_sketch_out(iter).model = model_photo;
            photo_sketch_out(iter).decision_values = photo_sketch_decision_values;
            
            sketch_photo_out(iter).true_labels = label_test;
            sketch_photo_out(iter).model = model_sketch;
            sketch_photo_out(iter).decision_values = sketch_photo_decision_values;
            
            drawing_sketch_out(iter).true_labels = label_test;
            drawing_sketch_out(iter).model = model_drawing;
            drawing_sketch_out(iter).decision_values = drawing_sketch_decision_values;
            
            sketch_drawing_out(iter).true_labels = label_test;
            sketch_drawing_out(iter).model = model_sketch;
            sketch_drawing_out(iter).decision_values = sketch_drawing_decision_values;
            
        end
        
        % get the accruacies 
        final_photo_drawing_acc(perm,time) = transres_accuracy_pairwise(photo_drawing_out, 0.5);
        final_drawing_photo_acc(perm,time) = transres_accuracy_pairwise(drawing_photo_out, 0.5);
        final_photo_sketch_acc(perm,time) = transres_accuracy_pairwise(photo_sketch_out, 0.5);
        final_sketch_photo_acc(perm,time) = transres_accuracy_pairwise(sketch_photo_out, 0.5);
        final_drawing_sketch_acc(perm,time) = transres_accuracy_pairwise(drawing_sketch_out, 0.5);
        final_sketch_drawing_acc(perm,time) = transres_accuracy_pairwise(sketch_drawing_out, 0.5);
        
        if mod(time,10)==0
            fprintf('.');
        end
    end
    fprintf('\nEnd permutation number: %i',perm);
end
% get the mean over all permutationss
final_photo_drawing_acc = mean(final_photo_drawing_acc,1);
final_drawing_photo_acc = mean(final_drawing_photo_acc,1);
final_photo_sketch_acc = mean(final_photo_sketch_acc,1);
final_sketch_photo_acc = mean(final_sketch_photo_acc,1);
final_drawing_sketch_acc = mean(final_drawing_sketch_acc,1);
final_sketch_drawing_acc = mean(final_sketch_drawing_acc,1);

fprintf('\ndone!')
end
