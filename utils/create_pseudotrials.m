% function create_pseudotrials(data, n_trials, n_average)
%
% This function is a helper function which can be used to compute
% pseudotrials from data in the format conditions x trials x channels x timepoints
%
% Input : data = data to compute pseudotrials from in the format
%         described above 
%         n_average = number of trials to average -> e.g if you have 24
%         trials in total per condition and you specify n_average =3 the
%         function will compute 24/3 = 8 pseudotrialsfor each of your
%         conditions
function pseudo_trialC = create_pseudotrials(data,n_average)

    L=size(data,2)/n_average; %must be a real number; number of pseudo trials
    pseudo_trialC=NaN(size(data,1),L,size(data,3),size(data,4)); %pre-allocate memory
    for step=1:L %average by steps
        trial_selector=(1+(step-1)*n_average):(n_average+(step-1)*n_average); %select trials to be averaged
        %pseudo_trialD % permutedD averaged into pseudo trials, i.e. of dimensions M * L * T
        pseudo_trialC(:,step,:,:)= mean(data(:,trial_selector,:,:),2); %assign pseudo trial to pseudo_trial_D
    end
end 