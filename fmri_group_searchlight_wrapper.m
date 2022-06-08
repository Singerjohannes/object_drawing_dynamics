%% script for aggregateting the searchlight decoding results 
clear all 
clc 

%setup paths 

path = pwd;
figure_path = fullfile(path,'figures');

% add utils 

addpath(fullfile(path,'utils'));

% add stats functions 

addpath(genpath(fullfile(path,'stats')));

% set plot defaults 

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3, 'DefaultAxesFontName', 'Helvetica','DefaultTextFontName', 'Helvetica') 

% get colormap 
cmap = colormap('redblueTecplot');
close all

%% load decoding results volumes 

data_dir = fullfile(path,'data/fmri/decoding/searchlight');

% load results volumes 

load(fullfile(data_dir,'photo_group_searchlight_results.mat'))
load(fullfile(data_dir,'drawing_group_searchlight_results.mat'))
load(fullfile(data_dir,'sketch_group_searchlight_results.mat'))

%% compute stats for decoding results --> this needs a lot of memory (~45GB RAM), skip this and load the precomputed results if you do not have enough RAM available 

nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'right';

% set rng to a fixed number 
rng(96);

sig_searchlight_photo = permutation_cluster_1sample_weight_alld_less_mem (photo_resultsvol-50, nperm, cluster_th, significance_th, tail);
sig_searchlight_drawing =  permutation_cluster_1sample_weight_alld_less_mem (drawing_resultsvol-50, nperm, cluster_th, significance_th, tail);
sig_searchlight_sketch =  permutation_cluster_1sample_weight_alld_less_mem (sketch_resultsvol-50, nperm, cluster_th, significance_th, tail);

% compute stats differences 
nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'both';

sig_searchlight_photo_drawing = permutation_cluster_1sample_weight_alld_less_mem (photo_resultsvol-drawing_resultsvol, nperm, cluster_th, significance_th, tail);
sig_searchlight_drawing_sketch =  permutation_cluster_1sample_weight_alld_less_mem (drawing_resultsvol-sketch_resultsvol, nperm, cluster_th, significance_th, tail);
sig_searchlight_photo_sketch =  permutation_cluster_1sample_weight_alld_less_mem (photo_resultsvol-sketch_resultsvol, nperm, cluster_th, significance_th, tail);

%% load precomputed significance maps 

load(fullfile(data_dir,'sig_searchlight_photo.mat'))
load(fullfile(data_dir,'sig_searchlight_drawing.mat'))
load(fullfile(data_dir,'sig_searchlight_sketch.mat'))

%% compute conjunction for all depictions

conjunction_decoding = sig_searchlight_photo&sig_searchlight_drawing&sig_searchlight_sketch;

%% plot decoding results 

plot_conjunction_additive_coloring(conjunction_decoding, sig_searchlight_photo, sig_searchlight_drawing, sig_searchlight_sketch, '',cmap)

%print(fullfile(results_path,'conjunction_map_with_overlaps_final_diff_coloring.svg'), ...
%              '-dsvg', '-r600')
          
%% load crossdecoding results volumes 

data_dir = fullfile(path,'data/fmri/crossdecoding/searchlight');

% load results volumes 

load(fullfile(data_dir,'photo_drawing_group_searchlight_results.mat'))
load(fullfile(data_dir,'drawing_sketch_group_searchlight_results.mat'))
load(fullfile(data_dir,'photo_sketch_group_searchlight_results.mat'))

%% compute statistics for crossdecoding with correction for multiple tests --> skip this step if you do not have enough RAM (~45GB) available 

nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'right';

% set rng to a fixed number 
rng(96);

[~,~,~,~,crossdecoding_photo_drawing_clusters,~,~,crossdecoding_photo_drawing_cluster_size, crossdecoding_photo_drawing_cluster_th] = permutation_cluster_1sample_weight_alld (photo_drawing_resultsvol-50, nperm, cluster_th, significance_th, tail);

[~,~,~,~,crossdecoding_drawing_sketch_clusters,~,~,crossdecoding_drawing_sketch_cluster_size, crossdecoding_drawing_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (drawing_sketch_resultsvol-50, nperm, cluster_th, significance_th, tail);

[~,~,~,~,crossdecoding_photo_sketch_clusters,~,~,crossdecoding_photo_sketch_cluster_size, crossdecoding_photo_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (photo_sketch_resultsvol-50, nperm, cluster_th, significance_th, tail);

%get maximum threshold for correcting for multiple comparisons
max_thr = max([crossdecoding_photo_drawing_cluster_th, crossdecoding_photo_sketch_cluster_th, crossdecoding_drawing_sketch_cluster_th]); 

% intialize significance arrays 
sz = size(photo_drawing_resultsvol);
sig_searchlight_photo_drawing = zeros(sz(2:end));
sig_searchlight_photo_sketch = zeros(sz(2:end));
sig_searchlight_drawing_sketch = zeros(sz(2:end));

%now threshold the clusters again based on the maximum threshold 
sig_searchlight_photo_drawing([crossdecoding_photo_drawing_clusters{crossdecoding_photo_drawing_cluster_size>max_thr}]) = 1;
sig_searchlight_photo_sketch([crossdecoding_photo_sketch_clusters{crossdecoding_photo_sketch_cluster_size>max_thr}]) = 1;
sig_searchlight_drawing_sketch([crossdecoding_drawing_sketch_clusters{crossdecoding_drawing_sketch_cluster_size>max_thr}]) = 1;

%% load precomputed significance maps 

load(fullfile(data_dir,'photo_drawing_group_searchlight_stats.mat'))
load(fullfile(data_dir,'drawing_sketch_group_searchlight_stats.mat'))
load(fullfile(data_dir,'photo_sketch_group_searchlight_stats.mat'))

%% compute conjunction for all comparisons

conjunction_crossdecoding = sig_searchlight_photo_drawing&sig_searchlight_drawing_sketch&sig_searchlight_photo_sketch;

%% plot crossdecoding results 

plot_conjunction_additive_coloring(conjunction_crossdecoding, sig_searchlight_photo_drawing, sig_searchlight_drawing_sketch, sig_searchlight_photo_sketch, '',cmap)

%print(fullfile(results_path,'conjunction_map_crossdecoding.jpeg'), ...
%              '-djpeg', '-r600')