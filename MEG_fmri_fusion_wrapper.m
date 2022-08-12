%% meg - fmri fusion wrapper

clear all
clc

%setup paths 

path = pwd;
figure_path = fullfile(path,'figures');

% add utils 

addpath(fullfile(path,'utils'));

% add stats functions 

addpath(genpath(fullfile(path,'stats')));

% setup the decoding toolbox 
try 
    decoding_defaults;
catch 
    tdt_path = input('The Decoding Toolbox seems to be not on your path. Please enter the path to your TDT version:\n','s');
    addpath(tdt_path);
    decoding_defaults;
end 

% set plot defaults 

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3, 'DefaultAxesFontName', 'Helvetica','DefaultTextFontName', 'Helvetica') 

% get colormap 
cmap = colormap('redblueTecplot');
close all

%set correlation type 
corr_type = 'pearson';

% set stats defaults 
nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'both';

%% load results 

fmri_data_dir = fullfile(path,'data/fmri/rsa/');

meg_data_dir =fullfile(path,'data/meg/rsa/');

load(fullfile(fmri_data_dir,'EVC_photo_RDM.mat'))
load(fullfile(fmri_data_dir,'EVC_drawing_RDM.mat'))
load(fullfile(fmri_data_dir,'EVC_sketch_RDM.mat'))

load(fullfile(fmri_data_dir,'LOC_photo_RDM.mat'))
load(fullfile(fmri_data_dir,'LOC_drawing_RDM.mat'))
load(fullfile(fmri_data_dir,'LOC_sketch_RDM.mat'))

load(fullfile(meg_data_dir,'meg_photo_RDM.mat'))
load(fullfile(meg_data_dir,'meg_drawing_RDM.mat'))
load(fullfile(meg_data_dir,'meg_sketch_RDM.mat'))

%% compute fusion - with mean fmri RDM 

for sub_no = 1:size(meg_photo_RDM,4)
    for time = 1:size(meg_photo_RDM,3)

    EVC_fusion(1,sub_no,time) = corr(squareformq(meg_photo_RDM(:,:,time,sub_no)'), mean(EVC_fmri_photo_RDM,2), 'Type', corr_type); 
    LOC_fusion(1,sub_no,time) = corr(squareformq(meg_photo_RDM(:,:,time,sub_no)'), mean(LOC_fmri_photo_RDM,2), 'Type', corr_type); 

    EVC_fusion(2,sub_no,time) = corr(squareformq(meg_drawing_RDM(:,:,time,sub_no)'), mean(EVC_fmri_drawing_RDM,2), 'Type', corr_type); 
    LOC_fusion(2,sub_no,time) = corr(squareformq(meg_drawing_RDM(:,:,time,sub_no)'), mean(LOC_fmri_drawing_RDM,2), 'Type', corr_type); 
        
    EVC_fusion(3,sub_no,time) = corr(squareformq(meg_sketch_RDM(:,:,time,sub_no)'), mean(EVC_fmri_sketch_RDM,2), 'Type', corr_type); 
    LOC_fusion(3,sub_no,time) = corr(squareformq(meg_sketch_RDM(:,:,time,sub_no)'), mean(LOC_fmri_sketch_RDM,2), 'Type', corr_type); 

    end 
end 

%% compute stats on fusion results 

% set rng to a fixed number 
rng(96);

sig_fusion_photo_EVC = permutation_cluster_1sample_weight_alld (squeeze(EVC_fusion(1,:,:)), nperm, cluster_th, significance_th, tail);

sig_fusion_drawing_EVC = permutation_cluster_1sample_weight_alld (squeeze(EVC_fusion(2,:,:)), nperm, cluster_th, significance_th, tail);

sig_fusion_sketch_EVC = permutation_cluster_1sample_weight_alld (squeeze(EVC_fusion(3,:,:)), nperm, cluster_th, significance_th, tail);

sig_fusion_photo_LOC = permutation_cluster_1sample_weight_alld (squeeze(LOC_fusion(1,:,:)), nperm, cluster_th, significance_th, tail);

sig_fusion_drawing_LOC = permutation_cluster_1sample_weight_alld (squeeze(LOC_fusion(2,:,:)), nperm, cluster_th, significance_th, tail);

sig_fusion_sketch_LOC = permutation_cluster_1sample_weight_alld (squeeze(LOC_fusion(3,:,:)), nperm, cluster_th, significance_th, tail);

% compute bootstraps for peak and onset

statsInfo.nperm = 10000;
statsInfo.cluster_th = 0.001;
statsInfo.significance_th = 0.05;
statsInfo.tail = 'both';
statsInfo.stat = [1 0];
nboot = 100000; 

EVC_photo_fusion_boot = bootstrap_fixed_1D(squeeze(EVC_fusion(1,:,:)), [-100:10:1000],nboot,statsInfo); 

EVC_drawing_fusion_boot = bootstrap_fixed_1D(squeeze(EVC_fusion(2,:,:)), [-100:10:1000],nboot,statsInfo); 

EVC_sketch_fusion_boot = bootstrap_fixed_1D(squeeze(EVC_fusion(3,:,:)), [-100:10:1000],nboot,statsInfo); 

LOC_photo_fusion_boot = bootstrap_fixed_1D(squeeze(LOC_fusion(1,:,:)), [-100:10:1000],nboot,statsInfo); 

LOC_drawing_fusion_boot = bootstrap_fixed_1D(squeeze(LOC_fusion(2,:,:)), [-100:10:1000],nboot,statsInfo); 

LOC_sketch_fusion_boot = bootstrap_fixed_1D(squeeze(LOC_fusion(3,:,:)), [-100:10:1000],nboot,statsInfo);

% bootstrap the difference for comparison 

EVC_photo_drawing_fusion_bootdiff = bootstrap_fixed_1D_diff(squeeze(EVC_fusion(1,:,:)),squeeze(EVC_fusion(2,:,:)), [-100:10:1000],nboot,statsInfo); 

EVC_photo_sketch_fusion_bootdiff = bootstrap_fixed_1D_diff(squeeze(EVC_fusion(1,:,:)),squeeze(EVC_fusion(3,:,:)), [-100:10:1000],nboot,statsInfo); 

EVC_drawing_sketch_fusion_bootdiff = bootstrap_fixed_1D_diff(squeeze(EVC_fusion(2,:,:)),squeeze(EVC_fusion(3,:,:)), [-100:10:1000],nboot,statsInfo); 

LOC_photo_drawing_fusion_bootdiff = bootstrap_fixed_1D_diff(squeeze(LOC_fusion(1,:,:)),squeeze(LOC_fusion(2,:,:)), [-100:10:1000],nboot,statsInfo); 

LOC_photo_sketch_fusion_bootdiff = bootstrap_fixed_1D_diff(squeeze(LOC_fusion(1,:,:)),squeeze(LOC_fusion(3,:,:)), [-100:10:1000],nboot,statsInfo); 

LOC_drawing_sketch_fusion_bootdiff = bootstrap_fixed_1D_diff(squeeze(LOC_fusion(2,:,:)),squeeze(LOC_fusion(3,:,:)), [-100:10:1000],nboot,statsInfo); 

% bootstrap difference between EVC and LOC peaks 

EVC_LOC_photo_bootdiff = bootstrap_fixed_1D_diff(squeeze(EVC_fusion(1,:,:)),squeeze(LOC_fusion(1,:,:)), [-100:10:1000],nboot,statsInfo);

EVC_LOC_drawing_bootdiff = bootstrap_fixed_1D_diff(squeeze(EVC_fusion(2,:,:)),squeeze(LOC_fusion(2,:,:)), [-100:10:1000],nboot,statsInfo);

EVC_LOC_sketch_bootdiff = bootstrap_fixed_1D_diff(squeeze(EVC_fusion(3,:,:)),squeeze(LOC_fusion(3,:,:)), [-100:10:1000],nboot,statsInfo);

%% calculate TOST for peak latencies for EVC

%first get peaks 

for sub = 1:size(EVC_fusion,2) 
    
    [~, EVC_photo_peaks(sub)] = max(EVC_fusion(1,sub,:)); 
    
    [~, EVC_drawing_peaks(sub)] = max(EVC_fusion(2,sub,:)); 

    [~, EVC_sketch_peaks(sub)] = max(EVC_fusion(3,sub,:)); 

end

[photo_drawing_equ_p, stat] = tost('one_sample', [-1 1], EVC_photo_peaks-EVC_drawing_peaks);

[photo_sketch_equ_p, stat] = tost('one_sample', [-1 1], EVC_photo_peaks-EVC_sketch_peaks);

[drawing_sketch_equ_p, stat] = tost('one_sample', [-1 1], EVC_drawing_peaks-EVC_sketch_peaks);

[~,~,~, adj_equ_p_EVC] = fdr_bh([photo_drawing_equ_p photo_sketch_equ_p drawing_sketch_equ_p]);


%% calculate TOST for peak latencies for LOC

%first get peaks 

for sub = 1:size(LOC_fusion,2) 
    
    [~, LOC_photo_peaks(sub)] = max(LOC_fusion(1,sub,:)); 
    
    [~, LOC_drawing_peaks(sub)] = max(LOC_fusion(2,sub,:)); 

    [~, LOC_sketch_peaks(sub)] = max(LOC_fusion(3,sub,:)); 

end

[photo_drawing_equ_p, stat] = tost('one_sample', [-1 1], LOC_photo_peaks-LOC_drawing_peaks);

[photo_sketch_equ_p, stat] = tost('one_sample', [-1 1], LOC_photo_peaks-LOC_sketch_peaks);

[drawing_sketch_equ_p, stat] = tost('one_sample', [-1 1], LOC_drawing_peaks-LOC_sketch_peaks);

[~,~,~, adj_equ_p_LOC] = fdr_bh([photo_drawing_equ_p photo_sketch_equ_p drawing_sketch_equ_p]);


%% plot fusion results

fig = figure;
options = [];
options.handle = fig;
options.x_axis = linspace(-0.1,1,111);
options.error = 'sem';
options.color_area = 'black';%[128 193 219]./255;    % Blue theme
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(squeeze(EVC_fusion(1,:,:)),options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(squeeze(EVC_fusion(2,:,:)),options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(squeeze(EVC_fusion(3,:,:)),options);
% plot stats 
sig_fusion_photo_EVC(sig_fusion_photo_EVC==0)= NaN; 
plot(options.x_axis,sig_fusion_photo_EVC*-0.01, 'black');
sig_fusion_drawing_EVC(sig_fusion_drawing_EVC==0)= NaN; 
plot(options.x_axis,sig_fusion_drawing_EVC*-0.02, 'color',cmap(ceil(256),:));
sig_fusion_sketch_EVC(sig_fusion_sketch_EVC==0)= NaN; 
plot(options.x_axis,sig_fusion_sketch_EVC*-0.03, 'color',cmap(ceil(200),:));
xlim([-0.1 1])
ylim([-0.05 0.6])
yline(0,':');
for i= 0:0.1:1
xline(i,':');
end
legend(this_line, 'Photo', 'Drawing', 'Sketch')
title(['EVC'])
ylabel('Pearson Correlation')
xlabel('Time (s)')

print(fullfile(figure_path, ['EVC_fusion.svg']), ...
             '-dsvg', '-r600')

fig2 = figure;
options = [];
options.handle = fig2;
options.x_axis = linspace(-0.1,1,111);
options.error = 'sem';
options.color_area = 'black';%[128 193 219]./255;    % Blue theme
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(squeeze(LOC_fusion(1,:,:)),options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(squeeze(LOC_fusion(2,:,:)),options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(squeeze(LOC_fusion(3,:,:)),options);
sig_fusion_photo_LOC(sig_fusion_photo_LOC==0)= NaN; 
plot(options.x_axis,sig_fusion_photo_LOC*-0.01, 'black');
sig_fusion_drawing_LOC(sig_fusion_drawing_LOC==0)= NaN; 
plot(options.x_axis,sig_fusion_drawing_LOC*-0.02, 'color',cmap(ceil(256),:));
sig_fusion_sketch_LOC(sig_fusion_sketch_LOC==0)= NaN; 
plot(options.x_axis,sig_fusion_sketch_LOC*-0.03, 'color',cmap(ceil(200),:));
yline(0,':');
for i= 0:0.1:1
xline(i,':');
end
xlim([-0.1 1])
ylim([-0.05 0.6])
legend(this_line, 'Photo', 'Drawing', 'Sketch')
title({'LOC'})
ylabel('Pearson Correlation')
xlabel('Time (s)')

print(fullfile(figure_path, ['LO_fusion.svg']), ...
             '-dsvg', '-r600')

%% compute differences between conditions 

% set rng to a fixed number 
rng(96);

EVC_fusion_photo_drawing = squeeze(EVC_fusion(1,:,:) - EVC_fusion(2,:,:));
EVC_fusion_photo_sketch = squeeze(EVC_fusion(1,:,:) - EVC_fusion(3,:,:));
EVC_fusion_drawing_sketch = squeeze(EVC_fusion(2,:,:) - EVC_fusion(3,:,:));

LOC_fusion_photo_drawing = squeeze(LOC_fusion(1,:,:) - LOC_fusion(2,:,:));
LOC_fusion_photo_sketch = squeeze(LOC_fusion(1,:,:) - LOC_fusion(3,:,:));
LOC_fusion_drawing_sketch = squeeze(LOC_fusion(2,:,:) - LOC_fusion(3,:,:));

% run stats 

[~,~,~,~,EVC_photo_drawing_clusters,EVC_photo_drawing_cluster_size, EVC_photo_drawing_cluster_th] = permutation_cluster_1sample_weight_alld (EVC_fusion_photo_drawing, nperm, cluster_th, significance_th, tail);

[~,~,~,~,EVC_photo_sketch_clusters,EVC_photo_sketch_cluster_size, EVC_photo_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (EVC_fusion_photo_sketch, nperm, cluster_th, significance_th, tail);

[~,~,~,~,EVC_drawing_sketch_clusters,EVC_drawing_sketch_cluster_size, EVC_drawing_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (EVC_fusion_drawing_sketch, nperm, cluster_th, significance_th, tail);

[~,~,~,~,LOC_photo_drawing_clusters,LOC_photo_drawing_cluster_size, LOC_photo_drawing_cluster_th] = permutation_cluster_1sample_weight_alld (LOC_fusion_photo_drawing, nperm, cluster_th, significance_th, tail);

[~,~,~,~,LOC_photo_sketch_clusters,LOC_photo_sketch_cluster_size, LOC_photo_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (LOC_fusion_photo_sketch, nperm, cluster_th, significance_th, tail);

[~,~,~,~,LOC_drawing_sketch_clusters,LOC_drawing_sketch_cluster_size, LOC_drawing_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (LOC_fusion_drawing_sketch, nperm, cluster_th, significance_th, tail);

% correct for multiple comparisons for each ROI separately

%get maximum threshold for correcting for multiple comparisons
max_thr_EVC = max([EVC_photo_drawing_cluster_th, EVC_photo_sketch_cluster_th, EVC_drawing_sketch_cluster_th]); 
max_thr_LOC = max([LOC_photo_drawing_cluster_th, LOC_photo_sketch_cluster_th, LOC_drawing_sketch_cluster_th]); 

% intialize significance arrays 
sig_photo_drawing_EVC = zeros(size(sig_fusion_drawing_EVC));
sig_photo_sketch_EVC =  zeros(size(sig_fusion_drawing_EVC));
sig_drawing_sketch_EVC =  zeros(size(sig_fusion_drawing_EVC));
sig_photo_drawing_LOC =  zeros(size(sig_fusion_drawing_EVC));
sig_photo_sketch_LOC =  zeros(size(sig_fusion_drawing_EVC));
sig_drawing_sketch_LOC =  zeros(size(sig_fusion_drawing_EVC));

%now threshold the clusters again based on the maximum threshold 
sig_photo_drawing_EVC([EVC_photo_drawing_clusters{EVC_photo_drawing_cluster_size>max_thr_EVC}]) = 1;
sig_photo_sketch_EVC([EVC_photo_sketch_clusters{EVC_photo_sketch_cluster_size>max_thr_EVC}]) = 1;
sig_drawing_sketch_EVC([EVC_drawing_sketch_clusters{EVC_drawing_sketch_cluster_size>max_thr_EVC}]) = 1;
sig_photo_drawing_LOC([LOC_photo_drawing_clusters{LOC_photo_drawing_cluster_size>max_thr_LOC}]) = 1;
sig_photo_sketch_LOC([LOC_photo_sketch_clusters{LOC_photo_sketch_cluster_size>max_thr_LOC}]) = 1;
sig_drawing_sketch_LOC([LOC_drawing_sketch_clusters{LOC_drawing_sketch_cluster_size>max_thr_LOC}]) = 1;

%% plot results differences

fig = figure;
options = [];
options.handle = fig;
options.x_axis = linspace(-0.1,1,111);
options.error = 'sem';
options.color_area = 'black';%[128 193 219]./255;    % Blue theme
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(EVC_fusion_photo_drawing,options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(EVC_fusion_photo_sketch,options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(EVC_fusion_drawing_sketch,options);
% plot stats 
sig_photo_drawing_EVC(sig_photo_drawing_EVC==0)= NaN; 
plot(options.x_axis,sig_photo_drawing_EVC*-0.02, 'black');
sig_photo_sketch_EVC(sig_photo_sketch_EVC==0)= NaN; 
plot(options.x_axis,sig_photo_sketch_EVC*-0.03, 'color',cmap(ceil(256),:));
sig_drawing_sketch_EVC(sig_drawing_sketch_EVC==0)= NaN; 
plot(options.x_axis,sig_drawing_sketch_EVC*-0.04, 'color',cmap(ceil(200),:));
xlim([-0.1 1])
ylim([-0.08 0.2])
yline(0,':');
for i= 0:0.1:1
xline(i,':');
end
legend(this_line, 'Photo-Drawing', 'Photo-Sketch', 'Drawing-Sketch')
title(['EVC'])
ylabel('Pearson Correlation Diff.')
xlabel('Time (s)')

print(fullfile(figure_path, ['EVC_fusion_diff.svg']), ...
             '-dsvg', '-r600')

fig2 = figure;
options = [];
options.handle = fig2;
options.x_axis = linspace(-0.1,1,111);
options.error = 'sem';
options.color_area = 'black';%[128 193 219]./255;    % Blue theme
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(LOC_fusion_photo_drawing,options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(LOC_fusion_photo_sketch,options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(LOC_fusion_drawing_sketch,options);
sig_photo_drawing_LOC(sig_photo_drawing_LOC==0)= NaN; 
plot(options.x_axis,sig_photo_drawing_LOC*-0.02, 'black');
sig_photo_sketch_LOC(sig_photo_sketch_LOC==0)= NaN; 
plot(options.x_axis,sig_photo_sketch_LOC*-0.03, 'color',cmap(ceil(256),:));
sig_drawing_sketch_LOC(sig_drawing_sketch_LOC==0)= NaN; 
plot(options.x_axis,sig_drawing_sketch_LOC*-0.04, 'color',cmap(ceil(200),:));
yline(0,':');
for i= 0:0.1:1
xline(i,':');
end
xlim([-0.1 1])
ylim([-0.08 0.2])
legend(this_line, 'Photo-Drawing', 'Photo-Sketch', 'Drawing-Sketch')
title('LOC')
ylabel('Pearson Correlation Diff.')
xlabel('Time (s)')

print(fullfile(figure_path, ['LOC_fusion_diff.svg']), ...
             '-dsvg', '-r600')