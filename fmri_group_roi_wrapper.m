%% script for aggregating the fmri ROI results 

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


%% load the results 

data_dir = fullfile(path,'data/fmri/decoding/roi'); 

load(fullfile(data_dir,'photo_group_roi_results.mat'))
load(fullfile(data_dir,'drawing_group_roi_results.mat'))
load(fullfile(data_dir,'sketch_group_roi_results.mat'))

%% plot 

roi_names = {'EVC'; 'LOC';'pIPS'};

% bring data in right format for plotting
all_accs = cat(2, mean(photo_group_decoding)', mean(drawing_group_decoding)', mean(sketch_group_decoding)');

for roi = 1:size(photo_group_decoding,2)
photo_se(roi) = [std(photo_group_decoding(:,roi))/sqrt(length(photo_group_decoding))];
drawing_se(roi) = [std(drawing_group_decoding(:,roi))/sqrt(length(drawing_group_decoding))];
sketch_se(roi) = [std(sketch_group_decoding(:,roi))/sqrt(length(sketch_group_decoding))];
end      
all_se = cat(2, photo_se', drawing_se', sketch_se');

% plot
figure
h = bar(all_accs-50, 'grouped','FaceColor', 'flat');
h(1).CData = [0 0 0];
h(2).CData = cmap(ceil(256),:);
h(3).CData = cmap(ceil(200),:);
xticklabels([roi_names])
yticks([0:5:35])
yticklabels([50:5:85])
ylim([0 35])
xlabel('ROI')
ylabel('Decoding Accuracy (%)')
title('Category Decoding - fMRI')

hold on
% Find the number of groups and the number of bars in each group

ngroups = size(all_accs, 1);
nbars = size(all_accs, 2);
% Calculate the width for each bar group

groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, all_accs(:,i)-50, all_se(:,i), 'k', 'linestyle', 'none');
end
legend({'Photos'; 'Drawings'; 'Sketches'} ,'Location','northeast')

print(fullfile(figure_path, ['cat_decoding_ROI.svg']), ...
              '-dsvg', '-r600')

%% compute statistics 

% set stats defaults 
nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'right';

% set rng to a fixed number 
rng(96);

sig_decoding_photo_EVC = permutation_1sample_alld (photo_group_decoding(:,1)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_drawing_EVC = permutation_1sample_alld (drawing_group_decoding(:,1)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_sketch_EVC = permutation_1sample_alld (sketch_group_decoding(:,1)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_photo_LOC = permutation_1sample_alld (photo_group_decoding(:,2)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_drawing_LOC = permutation_1sample_alld (drawing_group_decoding(:,2)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_sketch_LOC = permutation_1sample_alld (sketch_group_decoding(:,2)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_photo_pIPS = permutation_1sample_alld (photo_group_decoding(:,3)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_drawing_pIPS = permutation_1sample_alld (drawing_group_decoding(:,3)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_sketch_pIPS = permutation_1sample_alld (sketch_group_decoding(:,3)-50, nperm, cluster_th, significance_th, tail);


% compute statistics on differences between types of depiction

tail = 'both';

sig_photo_drawing_EVC = permutation_1sample_alld (photo_group_decoding(:,1)-drawing_group_decoding(:,1), nperm, cluster_th, significance_th, tail);

sig_photo_sketch_EVC = permutation_1sample_alld (photo_group_decoding(:,1)-sketch_group_decoding(:,1), nperm, cluster_th, significance_th, tail);

sig_drawing_sketch_EVC = permutation_1sample_alld (drawing_group_decoding(:,1)-sketch_group_decoding(:,1), nperm, cluster_th, significance_th, tail);


sig_photo_drawing_LOC = permutation_1sample_alld (photo_group_decoding(:,2)-drawing_group_decoding(:,2), nperm, cluster_th, significance_th, tail);

sig_photo_sketch_LOC = permutation_1sample_alld (photo_group_decoding(:,2)-sketch_group_decoding(:,2), nperm, cluster_th, significance_th, tail);

sig_drawing_sketch_LOC = permutation_1sample_alld (drawing_group_decoding(:,2)-sketch_group_decoding(:,2), nperm, cluster_th, significance_th, tail);

sig_photo_drawing_pIPS = permutation_1sample_alld (photo_group_decoding(:,3)-drawing_group_decoding(:,3), nperm, cluster_th, significance_th, tail);

sig_photo_sketch_pIPS = permutation_1sample_alld (photo_group_decoding(:,3)-sketch_group_decoding(:,3), nperm, cluster_th, significance_th, tail);

sig_drawing_sketch_pIPS = permutation_1sample_alld (drawing_group_decoding(:,3)-sketch_group_decoding(:,3), nperm, cluster_th, significance_th, tail);

% compute statistics on differences between ROIs

tail = 'both';

sig_photo_EVC_LOC = permutation_1sample_alld (photo_group_decoding(:,1)-photo_group_decoding(:,2), nperm, cluster_th, significance_th, tail);

sig_drawing_EVC_LOC = permutation_1sample_alld (drawing_group_decoding(:,1)-drawing_group_decoding(:,2), nperm, cluster_th, significance_th, tail);

sig_sketch_EVC_LOC = permutation_1sample_alld (sketch_group_decoding(:,1)-sketch_group_decoding(:,2), nperm, cluster_th, significance_th, tail);

sig_photo_EVC_pIPS = permutation_1sample_alld (photo_group_decoding(:,1)-photo_group_decoding(:,3), nperm, cluster_th, significance_th, tail);

sig_drawing_EVC_pIPS = permutation_1sample_alld (drawing_group_decoding(:,1)-drawing_group_decoding(:,3), nperm, cluster_th, significance_th, tail);

sig_sketch_EVC_pIPS = permutation_1sample_alld (sketch_group_decoding(:,1)-sketch_group_decoding(:,3), nperm, cluster_th, significance_th, tail);

sig_photo_LOC_pIPS = permutation_1sample_alld (photo_group_decoding(:,2)-photo_group_decoding(:,3), nperm, cluster_th, significance_th, tail);

sig_drawing_LOC_pIPS = permutation_1sample_alld (drawing_group_decoding(:,2)-drawing_group_decoding(:,3), nperm, cluster_th, significance_th, tail);

sig_sketch_LOC_pIPS = permutation_1sample_alld (sketch_group_decoding(:,2)-sketch_group_decoding(:,3), nperm, cluster_th, significance_th, tail);


% control for multiple comparisons

[~,~,~,adj_p_diff_EVC] = fdr_bh([sig_photo_drawing_EVC sig_photo_sketch_EVC sig_drawing_sketch_EVC]);
[~,~,~,adj_p_diff_LOC] =  fdr_bh([sig_photo_drawing_LOC sig_photo_sketch_LOC sig_drawing_sketch_LOC]);
[~,~,~,adj_p_diff_pIPS] =  fdr_bh([sig_photo_drawing_pIPS sig_photo_sketch_pIPS sig_drawing_sketch_pIPS]);

[~,~,~,adj_p_diff_between_regions] =  fdr_bh([sig_photo_EVC_LOC sig_drawing_EVC_LOC sig_sketch_EVC_LOC,...
                                              sig_photo_EVC_pIPS sig_drawing_EVC_pIPS sig_sketch_EVC_pIPS,...
                                              sig_photo_LOC_pIPS sig_drawing_LOC_pIPS sig_sketch_LOC_pIPS]);

%% load the crossdecoding results 

data_dir = fullfile(path,'data/fmri/crossdecoding/roi'); 

load(fullfile(data_dir,'photo_drawing_group_roi_results.mat'))
load(fullfile(data_dir,'drawing_sketch_group_roi_results.mat'))
load(fullfile(data_dir,'photo_sketch_group_roi_results.mat'))

%% plot 

roi_names = {'EVC'; 'LOC';'pIPS'};

%bring data in right format for plotting
all_accs = cat(2, mean(photo_drawing_group_decoding)', mean(drawing_sketch_group_decoding)', mean(photo_sketch_group_decoding)');

for roi = 1:size(photo_drawing_group_decoding,2)
photo_se(roi) = [std(photo_drawing_group_decoding(:,roi))/sqrt(length(photo_drawing_group_decoding))];
drawing_se(roi) = [std(drawing_sketch_group_decoding(:,roi))/sqrt(length(drawing_sketch_group_decoding))];
sketch_se(roi) = [std(photo_sketch_group_decoding(:,roi))/sqrt(length(photo_sketch_group_decoding))];
end

all_se = cat(2, photo_se', drawing_se', sketch_se');

% plot
figure
h = bar(all_accs-50, 'grouped','FaceColor', 'flat');
h(1).CData = [0 0 0];
h(2).CData = cmap(ceil(256),:);
h(3).CData = cmap(ceil(200),:);
xticklabels([roi_names])
yticks([0:5:35])
yticklabels([50:5:85])
ylim([0 35])
xlabel('ROI')
ylabel('Decoding Accuracy (%)')
title('Category Crossdecoding - fMRI')

hold on
% Find the number of groups and the number of bars in each group

ngroups = size(all_accs, 1);
nbars = size(all_accs, 2);
% Calculate the width for each bar group

groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, all_accs(:,i)-50, all_se(:,i), 'k', 'linestyle', 'none');
end
legend({'Photo-Drawing'; 'Drawing-Sketch'; 'Photo-Sketch'} ,'Location','northeast')

print(fullfile(figure_path, ['cat_crossdecoding_ROI.svg']), ...
              '-dsvg', '-r600')

%% compute statistics 

% set stats defaults 
nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'right';

% set rng to a fixed number 
rng(96);

sig_decoding_photo_drawing_EVC = permutation_1sample_alld (photo_drawing_group_decoding(:,1)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_photo_sketch_EVC = permutation_1sample_alld (photo_sketch_group_decoding(:,1)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_drawing_sketch_EVC = permutation_1sample_alld (drawing_sketch_group_decoding(:,1)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_photo_drawing_LOC = permutation_1sample_alld (photo_drawing_group_decoding(:,2)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_photo_sketch_LOC = permutation_1sample_alld (photo_sketch_group_decoding(:,2)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_drawing_sketch_LOC = permutation_1sample_alld (drawing_sketch_group_decoding(:,2)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_photo_drawing_pIPS = permutation_1sample_alld (photo_drawing_group_decoding(:,3)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_photo_sketch_pIPS = permutation_1sample_alld (photo_sketch_group_decoding(:,3)-50, nperm, cluster_th, significance_th, tail);

sig_decoding_drawing_sketch_pIPS = permutation_1sample_alld (drawing_sketch_group_decoding(:,3)-50, nperm, cluster_th, significance_th, tail);


% control for multiple comparisons

[~,~,~,adj_p_EVC] = fdr_bh([sig_decoding_photo_drawing_EVC sig_decoding_photo_sketch_EVC sig_decoding_drawing_sketch_EVC]);
[~,~,~,adj_p_LOC] = fdr_bh([sig_decoding_photo_drawing_LOC sig_decoding_photo_sketch_LOC sig_decoding_drawing_sketch_LOC]);
[~,~,~,adj_p_pIPS] = fdr_bh([sig_decoding_photo_drawing_pIPS sig_decoding_photo_sketch_pIPS sig_decoding_drawing_sketch_pIPS]);

%% compute the difference between straight decoding and cross-decoding and plot


roi_names = {'EVC'; 'LOC';'pIPS'};

photo_minus_photo_drawing=  photo_group_decoding-photo_drawing_group_decoding;
photo_minus_photo_sketch =  photo_group_decoding-photo_sketch_group_decoding;
drawing_minus_photo_drawing =  drawing_group_decoding-photo_drawing_group_decoding;
drawing_minus_drawing_sketch =  drawing_group_decoding-drawing_sketch_group_decoding;
sketch_minus_drawing_sketch =  sketch_group_decoding-drawing_sketch_group_decoding;
sketch_minus_photo_sketch =  sketch_group_decoding-photo_sketch_group_decoding;

%bring data in right format for plotting
all_accs = cat(2,mean(photo_minus_photo_drawing)',mean(photo_minus_photo_sketch)',mean(drawing_minus_photo_drawing)',mean(drawing_minus_drawing_sketch)',...
                mean(sketch_minus_photo_sketch)',mean(sketch_minus_drawing_sketch)');

all_se = [];          
for roi = 1:length(roi_names)
    
    all_se(roi,:) = [std(photo_minus_photo_drawing(:,roi))/sqrt(length(photo_minus_photo_drawing)),...
              std(photo_minus_photo_sketch(:,roi))/sqrt(length(photo_minus_photo_sketch)),...
              std(drawing_minus_photo_drawing(:,roi))/sqrt(length(drawing_minus_photo_drawing)),...
              std(drawing_minus_drawing_sketch(:,roi))/sqrt(length(drawing_minus_drawing_sketch)),...
              std(sketch_minus_drawing_sketch(:,roi))/sqrt(length(sketch_minus_drawing_sketch)),...
              std(sketch_minus_photo_sketch(:,roi))/sqrt(length(sketch_minus_photo_sketch))]';
end 
% plot
figure
h = bar(all_accs, 'grouped','FaceColor', 'flat');
h(1).CData = [0 0 0];
h(2).CData = cmap(ceil(256),:);
h(3).CData = cmap(ceil(200),:);
h(4).CData = cmap(ceil(170),:);
h(5).CData = cmap(ceil(50),:);
h(6).CData = cmap(ceil(1),:);
xticklabels([roi_names])
yticks([0:5:20])
yticklabels([0:5:20])
ylim([0 15])
xlabel('ROI')
ylabel({'Decoding accuracy','difference (%)'})

hold on
% Find the number of groups and the number of bars in each group

ngroups = size(all_accs, 1);
nbars = size(all_accs, 2);
% Calculate the width for each bar group

groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, all_accs(:,i), all_se(:,i), 'k', 'linestyle', 'none');
end
legend({'Photo minus Photo-Drawing'; 'Photo minus Photo-Sketch'; 'Drawing minus Photo-Drawing'; 'Drawing minus Drawing-Sketch'; 'Sketch minus Photo-Sketch';'Sketch minus Drawing-Sketch'} ,'Location','northeast')

print(fullfile(figure_path, ['straight_vs_cross_decoding_ROI.svg']), ...
              '-dsvg', '-r600')

%% compute statistics for the differences 

% set stats defaults 
nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'right';

% set rng to a fixed number 
rng(96);

for roi = 1:size(photo_group_decoding,2)
    
sig_decoding_photo_minus_photo_drawing(roi) = permutation_1sample_alld (photo_minus_photo_drawing(:,roi), nperm, cluster_th, significance_th, tail);

sig_decoding_photo_minus_photo_sketch(roi) = permutation_1sample_alld (photo_minus_photo_sketch(:,roi), nperm, cluster_th, significance_th, tail);

sig_decoding_drawing_minus_photo_drawing(roi) = permutation_1sample_alld (drawing_minus_photo_drawing(:,roi), nperm, cluster_th, significance_th, tail);

sig_decoding_drawing_minus_drawing_sketch(roi) = permutation_1sample_alld (drawing_minus_drawing_sketch(:,roi), nperm, cluster_th, significance_th, tail);

sig_decoding_sketch_minus_photo_sketch(roi) = permutation_1sample_alld (sketch_minus_photo_sketch(:,roi), nperm, cluster_th, significance_th, tail);

sig_decoding_sketch_minus_drawing_sketch(roi) = permutation_1sample_alld (sketch_minus_drawing_sketch(:,roi), nperm, cluster_th, significance_th, tail);

end 
% control for multiple comparisons

[~,~,~,adj_p_EVC] = fdr_bh([sig_decoding_photo_minus_photo_drawing(1) sig_decoding_photo_minus_photo_sketch(1) sig_decoding_drawing_minus_photo_drawing(1) sig_decoding_drawing_minus_drawing_sketch(1) sig_decoding_sketch_minus_photo_sketch(1) sig_decoding_sketch_minus_drawing_sketch(1)]);
[~,~,~,adj_p_LOC] = fdr_bh([sig_decoding_photo_minus_photo_drawing(2) sig_decoding_photo_minus_photo_sketch(2) sig_decoding_drawing_minus_photo_drawing(2) sig_decoding_drawing_minus_drawing_sketch(2) sig_decoding_sketch_minus_photo_sketch(2) sig_decoding_sketch_minus_drawing_sketch(2)]);
[~,~,~,adj_p_pIPS] = fdr_bh([sig_decoding_photo_minus_photo_drawing(3) sig_decoding_photo_minus_photo_sketch(3) sig_decoding_drawing_minus_photo_drawing(3) sig_decoding_drawing_minus_drawing_sketch(3) sig_decoding_sketch_minus_photo_sketch(3) sig_decoding_sketch_minus_drawing_sketch(3)]);


% now test for differences between EVC and LOC
photo_minus_photo_drawing_diff_EVC_LOC = permutation_1sample_alld (photo_minus_photo_drawing(:,1) - photo_minus_photo_drawing(:,2), nperm, cluster_th, significance_th, tail);
photo_minus_photo_sketch_diff_EVC_LOC = permutation_1sample_alld (photo_minus_photo_sketch(:,1) - photo_minus_photo_sketch(:,2), nperm, cluster_th, significance_th, tail);
drawing_minus_photo_drawing_diff_EVC_LOC = permutation_1sample_alld (drawing_minus_photo_drawing(:,1) - drawing_minus_photo_drawing(:,2), nperm, cluster_th, significance_th, tail);
drawing_minus_drawing_sketch_diff_EVC_LOC = permutation_1sample_alld (drawing_minus_drawing_sketch(:,1) - drawing_minus_drawing_sketch(:,2), nperm, cluster_th, significance_th, tail);
sketch_minus_photo_sketch_diff_EVC_LOC = permutation_1sample_alld (sketch_minus_photo_sketch(:,1) - sketch_minus_photo_sketch(:,2), nperm, cluster_th, significance_th, tail);
sketch_minus_drawing_sketch_diff_EVC_LOC = permutation_1sample_alld (sketch_minus_drawing_sketch(:,1) - sketch_minus_drawing_sketch(:,2), nperm, cluster_th, significance_th, tail);

[~,~,~,adj_p_diff_EVC_LOC] = fdr_bh([photo_minus_photo_drawing_diff_EVC_LOC, photo_minus_photo_sketch_diff_EVC_LOC, drawing_minus_photo_drawing_diff_EVC_LOC,...
                                     drawing_minus_drawing_sketch_diff_EVC_LOC, sketch_minus_photo_sketch_diff_EVC_LOC, sketch_minus_drawing_sketch_diff_EVC_LOC]);

% now test for differences between EVC and pIPS
photo_minus_photo_drawing_diff_EVC_pIPS = permutation_1sample_alld (photo_minus_photo_drawing(:,1) - photo_minus_photo_drawing(:,3), nperm, cluster_th, significance_th, tail);
photo_minus_photo_sketch_diff_EVC_pIPS = permutation_1sample_alld (photo_minus_photo_sketch(:,1) - photo_minus_photo_sketch(:,3), nperm, cluster_th, significance_th, tail);
drawing_minus_photo_drawing_diff_EVC_pIPS = permutation_1sample_alld (drawing_minus_photo_drawing(:,1) - drawing_minus_photo_drawing(:,3), nperm, cluster_th, significance_th, tail);
drawing_minus_drawing_sketch_diff_EVC_pIPS = permutation_1sample_alld (drawing_minus_drawing_sketch(:,1) - drawing_minus_drawing_sketch(:,3), nperm, cluster_th, significance_th, tail);
sketch_minus_photo_sketch_diff_EVC_pIPS = permutation_1sample_alld (sketch_minus_photo_sketch(:,1) - sketch_minus_photo_sketch(:,3), nperm, cluster_th, significance_th, tail);
sketch_minus_drawing_sketch_diff_EVC_pIPS = permutation_1sample_alld (sketch_minus_drawing_sketch(:,1) - sketch_minus_drawing_sketch(:,3), nperm, cluster_th, significance_th, tail);

[~,~,~,adj_p_diff_EVC_pIPS] = fdr_bh([photo_minus_photo_drawing_diff_EVC_pIPS, photo_minus_photo_sketch_diff_EVC_pIPS, drawing_minus_photo_drawing_diff_EVC_pIPS,...
                                     drawing_minus_drawing_sketch_diff_EVC_pIPS, sketch_minus_photo_sketch_diff_EVC_pIPS, sketch_minus_drawing_sketch_diff_EVC_pIPS]);
 
% now test for differences between LOC and pIPS  
photo_minus_photo_drawing_diff_LOC_pIPS = permutation_1sample_alld (photo_minus_photo_drawing(:,1) - photo_minus_photo_drawing(:,3), nperm, cluster_th, significance_th, tail);
photo_minus_photo_sketch_diff_LOC_pIPS = permutation_1sample_alld (photo_minus_photo_sketch(:,1) - photo_minus_photo_sketch(:,3), nperm, cluster_th, significance_th, tail);
drawing_minus_photo_drawing_diff_LOC_pIPS = permutation_1sample_alld (drawing_minus_photo_drawing(:,1) - drawing_minus_photo_drawing(:,3), nperm, cluster_th, significance_th, tail);
drawing_minus_drawing_sketch_diff_LOC_pIPS = permutation_1sample_alld (drawing_minus_drawing_sketch(:,1) - drawing_minus_drawing_sketch(:,3), nperm, cluster_th, significance_th, tail);
sketch_minus_photo_sketch_diff_LOC_pIPS = permutation_1sample_alld (sketch_minus_photo_sketch(:,1) - sketch_minus_photo_sketch(:,3), nperm, cluster_th, significance_th, tail);
sketch_minus_drawing_sketch_diff_LOC_pIPS = permutation_1sample_alld (sketch_minus_drawing_sketch(:,1) - sketch_minus_drawing_sketch(:,3), nperm, cluster_th, significance_th, tail);

[~,~,~,adj_p_diff_LOC_pIPS] = fdr_bh([photo_minus_photo_drawing_diff_LOC_pIPS, photo_minus_photo_sketch_diff_LOC_pIPS, drawing_minus_photo_drawing_diff_LOC_pIPS,...
                                     drawing_minus_drawing_sketch_diff_LOC_pIPS, sketch_minus_photo_sketch_diff_LOC_pIPS, sketch_minus_drawing_sketch_diff_LOC_pIPS]);

