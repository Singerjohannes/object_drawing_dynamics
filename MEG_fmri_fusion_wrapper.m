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

meg_data_dir = fullfile(path,'data/meg/rsa/');

load(fullfile(fmri_data_dir,'photo_RDM_ROI.mat'))
load(fullfile(fmri_data_dir,'drawing_RDM_ROI.mat'))
load(fullfile(fmri_data_dir,'sketch_RDM_ROI.mat'))

load(fullfile(meg_data_dir,'meg_photo_RDM.mat'))
load(fullfile(meg_data_dir,'meg_drawing_RDM.mat'))
load(fullfile(meg_data_dir,'meg_sketch_RDM.mat'))

%% compute fusion - with mean fmri RDM 

for sub_no = 1:size(meg_photo_RDM,3)
    for time = 1:size(meg_photo_RDM,2)
        for roi = 1:size(fmri_photo_RDM,2)

    fusion(1,sub_no,roi,time) = corr(meg_photo_RDM(:,time,sub_no), mean(fmri_photo_RDM(:,roi,:),3), 'Type', corr_type); 

    fusion(2,sub_no,roi,time) = corr(meg_drawing_RDM(:,time,sub_no), mean(fmri_drawing_RDM(:,roi,:),3), 'Type', corr_type); 
        
    fusion(3,sub_no,roi,time) = corr(meg_sketch_RDM(:,time,sub_no), mean(fmri_sketch_RDM(:,roi,:),3), 'Type', corr_type); 
        end
    end 
end 

%% compute stats on fusion results 

% set rng to a fixed number 
rng(96);


for roi = 1:size(fmri_photo_RDM,2)

sig_fusion_photo(roi,:) = permutation_cluster_1sample_weight_alld (squeeze(fusion(1,:,roi,:)), nperm, cluster_th, significance_th, tail);

sig_fusion_drawing(roi,:) = permutation_cluster_1sample_weight_alld (squeeze(fusion(2,:,roi,:)), nperm, cluster_th, significance_th, tail);

sig_fusion_sketch(roi,:) = permutation_cluster_1sample_weight_alld (squeeze(fusion(3,:,roi,:)), nperm, cluster_th, significance_th, tail);
 
% compute bootstraps for peak and onset

statsInfo.nperm = 10000;
statsInfo.cluster_th = 0.001;
statsInfo.significance_th = 0.05;
statsInfo.tail = 'both';
statsInfo.stat = [1 0];
nboot = 100000; 

if roi <3
    photo_fusion_boot(roi) = bootstrap_fixed_1D(squeeze(fusion(1,:,roi,:)), [-100:10:1000],nboot,statsInfo);
    
    drawing_fusion_boot(roi) = bootstrap_fixed_1D(squeeze(fusion(2,:,roi,:)), [-100:10:1000],nboot,statsInfo);
    
    sketch_fusion_boot(roi) = bootstrap_fixed_1D(squeeze(fusion(3,:,roi,:)), [-100:10:1000],nboot,statsInfo);
    
    
    photo_drawing_fusion_bootdiff(roi) = bootstrap_fixed_1D_diff(squeeze(fusion(1,:,roi,:)),squeeze(fusion(2,:,roi,:)), [-100:10:1000],nboot,statsInfo);
    
    photo_sketch_fusion_bootdiff(roi) = bootstrap_fixed_1D_diff(squeeze(fusion(1,:,roi,:)),squeeze(fusion(3,:,roi,:)), [-100:10:1000],nboot,statsInfo);
    
    drawing_sketch_fusion_bootdiff(roi) = bootstrap_fixed_1D_diff(squeeze(fusion(2,:,roi,:)),squeeze(fusion(3,:,roi,:)), [-100:10:1000],nboot,statsInfo);
    
else
    
    photo_fusion_boot(roi) = bootstrap_fixed_1D(squeeze(fusion(1,:,roi,1:56)), [-100:10:450],nboot,statsInfo); 

    drawing_fusion_boot(roi) = bootstrap_fixed_1D(squeeze(fusion(2,:,roi,1:56)), [-100:10:450],nboot,statsInfo); 

    sketch_fusion_boot(roi) = bootstrap_fixed_1D(squeeze(fusion(3,:,roi,1:56)), [-100:10:450],nboot,statsInfo); 


    photo_drawing_fusion_bootdiff(roi) = bootstrap_fixed_1D_diff(squeeze(fusion(1,:,roi,1:56)),squeeze(fusion(2,:,roi,1:56)), [-100:10:450],nboot,statsInfo); 

    photo_sketch_fusion_bootdiff(roi) = bootstrap_fixed_1D_diff(squeeze(fusion(1,:,roi,1:56)),squeeze(fusion(3,:,roi,1:56)), [-100:10:450],nboot,statsInfo); 

    drawing_sketch_fusion_bootdiff(roi) = bootstrap_fixed_1D_diff(squeeze(fusion(2,:,roi,1:56)),squeeze(fusion(3,:,roi,1:56)), [-100:10:450],nboot,statsInfo); 
end 
end

% bootstrap difference between EVC and LOC peaks 

EVC_LOC_photo_bootdiff = bootstrap_fixed_1D_diff(squeeze(fusion(1,:,1,:)),squeeze(fusion(1,:,2,:)), [-100:10:1000],nboot,statsInfo);

EVC_LOC_drawing_bootdiff = bootstrap_fixed_1D_diff(squeeze(fusion(2,:,1,:)),squeeze(fusion(2,:,2,:)), [-100:10:1000],nboot,statsInfo);

EVC_LOC_sketch_bootdiff = bootstrap_fixed_1D_diff(squeeze(fusion(3,:,1,:)),squeeze(fusion(3,:,2,:)), [-100:10:1000],nboot,statsInfo);

% bootstrap difference between EVC and pIPS peaks 

EVC_pIPS_photo_bootdiff = bootstrap_fixed_1D_diff(squeeze(fusion(1,:,1,1:56)),squeeze(fusion(1,:,3,1:56)), [-100:10:450],nboot,statsInfo);

EVC_pIPS_drawing_bootdiff = bootstrap_fixed_1D_diff(squeeze(fusion(2,:,1,1:56)),squeeze(fusion(2,:,3,1:56)), [-100:10:450],nboot,statsInfo);

EVC_pIPS_sketch_bootdiff = bootstrap_fixed_1D_diff(squeeze(fusion(3,:,1,1:56)),squeeze(fusion(3,:,3,1:51)), [-100:10:450],nboot,statsInfo);

% bootstrap difference between LOC and pIPS peaks 

LOC_pIPS_photo_bootdiff = bootstrap_fixed_1D_diff(squeeze(fusion(1,:,2,1:56)),squeeze(fusion(1,:,3,1:56)), [-100:10:450],nboot,statsInfo);

LOC_pIPS_drawing_bootdiff = bootstrap_fixed_1D_diff(squeeze(fusion(2,:,2,1:56)),squeeze(fusion(2,:,3,1:56)), [-100:10:450],nboot,statsInfo);

LOC_pIPS_sketch_bootdiff = bootstrap_fixed_1D_diff(squeeze(fusion(3,:,2,1:56)),squeeze(fusion(3,:,3,1:56)), [-100:10:450],nboot,statsInfo);


%% calculate TOST for peak latencies in EVC and LOC 

for roi = 1:2
%first get peaks 

for sub = 1:size(fusion,2) 
    
    [~, photo_peaks(sub)] = max(fusion(1,sub,roi,:)); 
    
    [~, drawing_peaks(sub)] = max(fusion(2,sub,roi,:)); 

    [~, sketch_peaks(sub)] = max(fusion(3,sub,roi,:)); 

end

[photo_drawing_equ_p(roi), stat] = tost('one_sample', [-1 1], photo_peaks-drawing_peaks);

[photo_sketch_equ_p(roi), stat] = tost('one_sample', [-1 1], photo_peaks-sketch_peaks);

[drawing_sketch_equ_p(roi), stat] = tost('one_sample', [-1 1], drawing_peaks-sketch_peaks);

[~,~,~, adj_equ_p(roi,:)] = fdr_bh([photo_drawing_equ_p(roi) photo_sketch_equ_p(roi) drawing_sketch_equ_p(roi)]);

end

%% calculate TOST for peak latencies in pIPS 

%first get peaks 

for roi = 1:size(fusion,3)-1
for sub = 1:size(fusion,2)
    
    [~, photo_peaks(sub)] = max(fusion(1,sub,roi,:));
    
    [~, drawing_peaks(sub)] = max(fusion(2,sub,roi,:));
    
    [~, sketch_peaks(sub)] = max(fusion(3,sub,roi,:));

end 
    
    [photo_drawing_equ_p(roi), stat] = tost('one_sample', [-1 1], photo_peaks-drawing_peaks);

    [photo_sketch_equ_p(roi), stat] = tost('one_sample', [-1 1], photo_peaks-sketch_peaks);

    [drawing_sketch_equ_p(roi), stat] = tost('one_sample', [-1 1], drawing_peaks-sketch_peaks);

    [~,~,~, adj_equ_p(roi,:)] = fdr_bh([photo_drawing_equ_p(roi) photo_sketch_equ_p(roi) drawing_sketch_equ_p(roi)]);

end

for sub = 1:size(fusion,2)
    
    [~, photo_peaks(sub)] = max(fusion(1,sub,3,1:51)); 
    
    [~, drawing_peaks(sub)] = max(fusion(2,sub,3,1:51)); 

    [~, sketch_peaks(sub)] = max(fusion(3,sub,3,1:51)); 
end 
    
[photo_drawing_equ_p(3), stat] = tost('one_sample', [-1 1], photo_peaks-drawing_peaks);

[photo_sketch_equ_p(3), stat] = tost('one_sample', [-1 1], photo_peaks-sketch_peaks);

[drawing_sketch_equ_p(3), stat] = tost('one_sample', [-1 1], drawing_peaks-sketch_peaks);

[~,~,~, adj_equ_p(3,:)] = fdr_bh([photo_drawing_equ_p(3) photo_sketch_equ_p(3) drawing_sketch_equ_p(3)]);


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
this_line(1) = plot_areaerrorbar(squeeze(fusion(1,:,1,:)),options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(squeeze(fusion(2,:,1,:)),options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(squeeze(fusion(3,:,1,:)),options);
% plot stats 
sig_fusion_photo(1,sig_fusion_photo(1,:)==0)= NaN; 
plot(options.x_axis,sig_fusion_photo(1,:)*-0.01, 'black');
sig_fusion_drawing(1,sig_fusion_drawing(1,:)==0)= NaN; 
plot(options.x_axis,sig_fusion_drawing(1,:)*-0.02, 'color',cmap(ceil(256),:));
sig_fusion_sketch(1,sig_fusion_sketch(1,:)==0)= NaN; 
plot(options.x_axis,sig_fusion_sketch(1,:)*-0.03, 'color',cmap(ceil(200),:));
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
this_line(1) = plot_areaerrorbar(squeeze(fusion(1,:,2,:)),options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(squeeze(fusion(2,:,2,:)),options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(squeeze(fusion(3,:,2,:)),options);
sig_fusion_photo(2,sig_fusion_photo(2,:)==0)= NaN; 
plot(options.x_axis,sig_fusion_photo(2,:)*-0.01, 'black');
sig_fusion_drawing(2,sig_fusion_drawing(2,:)==0)= NaN; 
plot(options.x_axis,sig_fusion_drawing(2,:)*-0.02, 'color',cmap(ceil(256),:));
sig_fusion_sketch(2,sig_fusion_sketch(2,:)==0)= NaN; 
plot(options.x_axis,sig_fusion_sketch(2,:)*-0.03, 'color',cmap(ceil(200),:));
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

print(fullfile(figure_path, ['LOC_fusion.svg']), ...
             '-dsvg', '-r600')
         
fig3 = figure;
options = [];
options.handle = fig3;
options.x_axis = linspace(-0.1,1,111);
options.error = 'sem';
options.color_area = 'black';%[128 193 219]./255;    % Blue theme
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(squeeze(fusion(1,:,3,:)),options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(squeeze(fusion(2,:,3,:)),options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(squeeze(fusion(3,:,3,:)),options);
sig_fusion_photo(3,sig_fusion_photo(3,:)==0)= NaN; 
plot(options.x_axis,sig_fusion_photo(3,:)*0.11, 'black');
sig_fusion_drawing(3,sig_fusion_drawing(3,:)==0)= NaN; 
plot(options.x_axis,sig_fusion_drawing(3,:)*0.12, 'color',cmap(ceil(256),:));
sig_fusion_sketch(3,sig_fusion_sketch(3,:)==0)= NaN; 
plot(options.x_axis,sig_fusion_sketch(3,:)*0.13, 'color',cmap(ceil(200),:));
yline(0,':');
for i= 0:0.1:1
xline(i,':');
end
xlim([-0.1 1])
ylim([-0.05 0.6])
legend(this_line, 'Photo', 'Drawing', 'Sketch')
title({'pIPS'})
ylabel('Pearson Correlation')
xlabel('Time (s)')

print(fullfile(figure_path, ['pIPS_fusion.svg']), ...
             '-dsvg', '-r600')

%% compute differences between conditions 

% set rng to a fixed number 
rng(96);

for roi = 1:size(fusion,3)
    
fusion_photo_drawing(roi,:,:) = squeeze(fusion(1,:,roi,:) - fusion(2,:,roi,:));
fusion_photo_sketch(roi,:,:) = squeeze(fusion(1,:,roi,:) - fusion(3,:,roi,:));
fusion_drawing_sketch(roi,:,:) = squeeze(fusion(2,:,roi,:) - fusion(3,:,roi,:));

% run stats 

[~,~,~,~,photo_drawing_clusters,photo_drawing_cluster_size, photo_drawing_cluster_th] = permutation_cluster_1sample_weight_alld (squeeze(fusion_photo_drawing(roi,:,:)), nperm, cluster_th, significance_th, tail);

[~,~,~,~,photo_sketch_clusters,photo_sketch_cluster_size, photo_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (squeeze(fusion_photo_sketch(roi,:,:)), nperm, cluster_th, significance_th, tail);

[~,~,~,~,drawing_sketch_clusters,drawing_sketch_cluster_size, drawing_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (squeeze(fusion_drawing_sketch(roi,:,:)), nperm, cluster_th, significance_th, tail);

% correct for multiple comparisons for each ROI separately

%get maximum threshold for correcting for multiple comparisons
max_thr = max([photo_drawing_cluster_th, photo_sketch_cluster_th, drawing_sketch_cluster_th]); 

% intialize significance arrays 
sig_photo_drawing(roi,:) = zeros(size(sig_fusion_drawing(1,:)));
sig_photo_sketch(roi,:) =  zeros(size(sig_fusion_drawing(1,:)));
sig_drawing_sketch(roi,:) =  zeros(size(sig_fusion_drawing(1,:)));

%now threshold the clusters again based on the maximum threshold 
sig_photo_drawing(roi,[photo_drawing_clusters{photo_drawing_cluster_size>max_thr}]) = 1;
sig_photo_sketch(roi,[photo_sketch_clusters{photo_sketch_cluster_size>max_thr}]) = 1;
sig_drawing_sketch(roi,[drawing_sketch_clusters{drawing_sketch_cluster_size>max_thr}]) = 1;
end 
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
this_line(1) = plot_areaerrorbar(squeeze(fusion_photo_drawing(1,:,:)),options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(squeeze(fusion_photo_sketch(1,:,:)),options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(squeeze(fusion_drawing_sketch(1,:,:)),options);
% plot stats 
sig_photo_drawing(1,sig_photo_drawing(1,:)==0)= NaN; 
plot(options.x_axis,sig_photo_drawing(1,:)*-0.02, 'black');
sig_photo_sketch(1,sig_photo_sketch(1,:)==0)= NaN; 
plot(options.x_axis,sig_photo_sketch(1,:)*-0.03, 'color',cmap(ceil(256),:));
sig_drawing_sketch(1,sig_drawing_sketch(1,:)==0)= NaN; 
plot(options.x_axis,sig_drawing_sketch(1,:)*-0.04, 'color',cmap(ceil(200),:));
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
this_line(1) = plot_areaerrorbar(squeeze(fusion_photo_drawing(2,:,:)),options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(squeeze(fusion_photo_sketch(2,:,:)),options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(squeeze(fusion_drawing_sketch(2,:,:)),options);
sig_photo_drawing(2,sig_photo_drawing(2,:)==0)= NaN; 
plot(options.x_axis,sig_photo_drawing(2,:)*-0.02, 'black');
sig_photo_sketch(2,sig_photo_sketch(2,:)==0)= NaN; 
plot(options.x_axis,sig_photo_sketch(2,:)*-0.03, 'color',cmap(ceil(256),:));
sig_drawing_sketch(2,sig_drawing_sketch(2,:)==0)= NaN; 
plot(options.x_axis,sig_drawing_sketch(2,:)*-0.04, 'color',cmap(ceil(200),:));
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
         
fig3 = figure;
options = [];
options.handle = fig3;
options.x_axis = linspace(-0.1,1,111);
options.error = 'sem';
options.color_area = 'black';%[128 193 219]./255;    % Blue theme
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(squeeze(fusion_photo_drawing(3,:,:)),options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(squeeze(fusion_photo_sketch(3,:,:)),options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(squeeze(fusion_drawing_sketch(3,:,:)),options);
sig_photo_drawing(3,sig_photo_drawing(3,:)==0)= NaN; 
plot(options.x_axis,sig_photo_drawing(3,:)*0.1, 'black');
sig_photo_sketch(3,sig_photo_sketch(3,:)==0)= NaN; 
plot(options.x_axis,sig_photo_sketch(3,:)*0.09, 'color',cmap(ceil(256),:));
sig_drawing_sketch(3,sig_drawing_sketch(3,:)==0)= NaN; 
plot(options.x_axis,sig_drawing_sketch(3,:)*0.08, 'color',cmap(ceil(200),:));
yline(0,':');
for i= 0:0.1:1
xline(i,':');
end
xlim([-0.1 1])
ylim([-0.08 0.2])
legend(this_line, 'Photo-Drawing', 'Photo-Sketch', 'Drawing-Sketch')
title('pIPS')
ylabel('Pearson Correlation Diff.')
xlabel('Time (s)')

print(fullfile(figure_path, ['pIPS_fusion_diff.svg']), ...
             '-dsvg', '-r600')