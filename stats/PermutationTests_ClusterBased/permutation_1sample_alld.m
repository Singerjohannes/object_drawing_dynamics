function [pval] = permutation_1sample_alld (data, nperm, cluster_th, significance_th, tail)
% Performs one-sided (>0) or two-sided permutation test on the 'data'.
% The data array can have any number of dimensions (1D, 2D, 3D, 4D etc) in
% the order [observations x variable1 x variable2 x ...]
% Data from each observation are randomly multiplied by +-1 to create
% permutation samples. These samples are converted to pvalues and then the
% cluster_th threshold (in pvalue units) is applied to identify suprathreshold clusters.
%
% INPUT:
%   data: observations x variable1 x variable2 x variable3 x ... (supports all dimemsions)
%   nperm: number of permutations
%   cluster_th: cluster defining threshold (in pvalue units)
%   significance_th: significance threshold (alpha value)
%   tail: string 'right' or 'both'. Default is 'right'.
%
% % for example:
%   load data;
%   nperm = 1000;
%   cluster_th = 0.05;
%   significance_th = 0.05;
%   tail = 'right';
%
% OUTPUT:
%   significantVarMax:  significant clusters (maximum cluster size)
%   significantVarWei:  significant clusters (maximum cluster weights)
%   clusters:           clusters above cluster defining threshold
%   pValMax:            p-value of maximum cluster size found in original data
%   pValWei:            p-value of maximum cluster weight found in original data
%
% Author: Dimitrios Pantazis, December 2015
% (modified): Siying Xie, August 2020

%decide one-sided (right) or two-sided (both) test
if ~exist('tail','var') || strcmp(tail,'right')
    func = '';
else %if two-sided test
    func = 'abs';
end

%initialize
N = ndims(data); %data is observations x variable1 x variable2 x ... (supports all dimensions)
nobservations = size(data,1);
for n = 2:N
    nvariable(n-1) = size(data,n);
end
cln = repmat({':'},1,N-1); %to select the N-1 dimensions

%create permutation samples and convert them to pvalues (perms x variable1 x variable2 x ...)
StatMapPerm = single(zeros([nperm nvariable]));

%first permutation sample is original data
StatMapPerm(1,cln{:}) = mean(data,1) ./ std(data);

%perform permutations
parfor i = 2:nperm
    if ~rem(i,100)
        disp(['Create permutation samples: ' num2str(i) ' out of ' num2str(nperm)]);
    end
    perm = single(sign(rand(nobservations,1)-0.5));
    %multiply by filling extra dimesions
    data_perm = repmat(perm,[1, nvariable]) .* data; % or data_perm = gmultiply(perm,data);
    permsample = mean(data_perm,1) ./ std(data_perm);
    StatMapPerm(i,:) = permsample(:); %does this for multiple dimensions: StatMapPerm(i,:,:) = mean(data_perm,1) ./ std(data_perm);
end

%convert to pvalues
eval([ 'StatMapPermPV = (nperm+1 - tiedrank(' func '(StatMapPerm)))/nperm;' ]);
StatMapPermPV = StatMapPermPV;

pval = StatMapPermPV(1); 
end