%% meg group analysis wrapper 
% this script contains all the group analyses for the decoding analyses of
% the meg data:
% 1. Category Decoding + Statistics
% 2. Category Cross Decoding + Statistics 
% 3. Temporal Generalization + Statistics 

clear all 
clc

%setup paths 

path = pwd;
figure_path = fullfile(path,'figures');
% create figure path
if ~isdir(figure_path); mkdir(figure_path); end  

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

% setup fieldtrip for later matrix plotting 
addpath(fullfile(path,'fieldtrip')); % change this to your fieldtrip path 
ft_defaults;

% set stats defaults 
nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'right';

%% load group decoding results

data_dir = fullfile(path,'data/meg/decoding/');

load(fullfile(data_dir,'photo_group_results.mat'))
load(fullfile(data_dir,'drawing_group_results.mat'))
load(fullfile(data_dir,'sketch_group_results.mat'))

%% compute stats on decoding results 

% set rng to a fixed number 
rng(96);

[sig_decoding_photo] = permutation_cluster_1sample_weight_alld (photo_group_acc-50, nperm, cluster_th, significance_th, tail);

[sig_decoding_drawing] = permutation_cluster_1sample_weight_alld (drawing_group_acc-50, nperm, cluster_th, significance_th, tail);

[sig_decoding_sketch] = permutation_cluster_1sample_weight_alld (sketch_group_acc-50, nperm, cluster_th, significance_th, tail);

%% bootstrap to get CIs for peak latency

statsInfo.nperm = 10000;
statsInfo.cluster_th = 0.001;
statsInfo.significance_th = 0.05;
statsInfo.tail = 'right';
statsInfo.stat = [1 0]; 
nboot = 100000; 

% set rng to a fixed number 
rng(96);

photo_decoding_boot = bootstrap_fixed_1D(photo_group_acc-50, [-100:10:1000],nboot,statsInfo); 

drawing_decoding_boot = bootstrap_fixed_1D(drawing_group_acc-50, [-100:10:1000],nboot,statsInfo); 

sketch_decoding_boot = bootstrap_fixed_1D(sketch_group_acc-50, [-100:10:1000],nboot,statsInfo); 

% bootstrap the difference for comparison of peak latencies 

photo_drawing_decoding_bootdiff = bootstrap_fixed_1D_diff(photo_group_acc-50,drawing_group_acc-50, [-100:1:1000],nboot,statsInfo); 

photo_sketch_decoding_bootdiff = bootstrap_fixed_1D_diff(photo_group_acc-50,sketch_group_acc-50, [-100:1:1000],nboot,statsInfo); 

drawing_sketch_decoding_bootdiff = bootstrap_fixed_1D_diff(drawing_group_acc-50,sketch_group_acc-50, [-100:1:1000],nboot,statsInfo); 

%% compute TOST for comparing peak latencies 

%first get peaks 

for sub = 1:size(photo_group_acc,1)
    
    [~, photo_peaks(sub)] = max(photo_group_acc(sub,:)); 
    
    [~, drawing_peaks(sub)] = max(drawing_group_acc(sub,:)); 

    [~, sketch_peaks(sub)] = max(sketch_group_acc(sub,:)); 

end

mean_diff_photo_sketch = mean(photo_peaks-sketch_peaks); 

mean_diff_photo_drawing = mean(photo_peaks-drawing_peaks); 

mean_diff_drawing_sketch = mean(drawing_peaks-sketch_peaks); 

[photo_drawing_equ_p, stat] = tost('one_sample', [-1 1], photo_peaks-drawing_peaks);

[photo_sketch_equ_p, stat] = tost('one_sample', [-1 1], photo_peaks-sketch_peaks);

[drawing_sketch_equ_p, stat] = tost('one_sample', [-1 1], drawing_peaks-sketch_peaks);

[~,~,~, adj_equ_p] = fdr_bh([photo_drawing_equ_p photo_sketch_equ_p drawing_sketch_equ_p]);

%% plot results
fig = figure;
options = [];
options.handle = fig;
options.x_axis = linspace(-0.1,1,111);
options.error = 'sem';
options.color_area = 'black';%[128 193 219]./255;    % Blue theme
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(photo_group_acc,options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(drawing_group_acc,options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(sketch_group_acc,options);
% plot stats 
sig_decoding_photo(sig_decoding_photo==0)= NaN; 
plot(options.x_axis,sig_decoding_photo*49, 'black');
sig_decoding_drawing(sig_decoding_drawing==0)= NaN; 
plot(options.x_axis,sig_decoding_drawing*48, 'color',cmap(ceil(256),:));
sig_decoding_sketch(sig_decoding_sketch==0)= NaN; 
plot(options.x_axis,sig_decoding_sketch*47, 'color',cmap(ceil(200),:));
yline(50,':');
for i= 0:0.1:1
xline(i,':');
end
xlim([-0.1 1])
ylim([45 95])
legend(this_line, 'Photos', 'Drawings', 'Sketches')
title('Category Decoding - MEG')
ylabel('Decoding Accuracy (%)')
xlabel('Time (s)')

print(fullfile(figure_path, ['group_cat_decoding.svg']), ...
              '-dsvg', '-r600')

%% single plots 

for i = 1:size(photo_group_acc,1)
    
    sub_id  = num2str(i);

    
    subplot(6,4,i)
    plot(photo_group_acc(i,:))
    hold on
    plot(drawing_group_acc(i,:))
    plot(sketch_group_acc(i,:))
    title(['Decoding for subject: ', sub_id]);
    ylim([45 100])
    yline(50,':');
    xticks(linspace(1,size(photo_group_acc,2),12));
    xticklabels(linspace(-0.1,1,12));
    hold off
end 
suptitle('Single subject decoding')
legend('Photos', 'Drawings', 'Sketches');

%% compute differences between decoding time courses

significance_th = 0.05;
tail = 'both';
% set rng to a fixed number 
rng(96);

photo_minus_drawing = photo_group_acc - drawing_group_acc;
photo_minus_sketch = photo_group_acc - sketch_group_acc;
drawing_minus_sketch = drawing_group_acc - sketch_group_acc; 

[~,~,~,~,photo_minus_drawing_clusters,~,~,photo_minus_drawing_cluster_size, photo_minus_drawing_cluster_th] = permutation_cluster_1sample_weight_alld (photo_minus_drawing, nperm, cluster_th, significance_th, tail);
[~,~,~,~,photo_minus_sketch_clusters,~,~,photo_minus_sketch_cluster_size, photo_minus_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (photo_minus_sketch, nperm, cluster_th, significance_th, tail);
[~,~,~,~, drawing_minus_sketch_clusters,~,~,drawing_minus_sketch_cluster_size, drawing_minus_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (drawing_minus_sketch, nperm, cluster_th, significance_th, tail);

%get maximum threshold for correcting for multiple comparisons
max_thr = max([photo_minus_drawing_cluster_th, photo_minus_sketch_cluster_th, drawing_minus_sketch_cluster_th]); 

% intialize significance arrays 
sig_photo_minus_drawing = zeros(1,size(photo_minus_drawing,2));
sig_photo_minus_sketch = zeros(1,size(photo_minus_drawing,2));
sig_drawing_minus_sketch = zeros(1,size(photo_minus_drawing,2));

%now threshold the clusters again based on the maximum threshold 
sig_photo_minus_drawing([photo_minus_drawing_clusters{photo_minus_drawing_cluster_size>max_thr}]) = 1;
sig_photo_minus_sketch([photo_minus_sketch_clusters{photo_minus_sketch_cluster_size>max_thr}]) = 1;
sig_drawing_minus_sketch([drawing_minus_sketch_clusters{drawing_minus_sketch_cluster_size>max_thr}]) = 1;

%% plot differences 

fig = figure;
options = [];
options.handle = fig;
options.x_axis = linspace(-0.1,1,111);
options.error = 'sem';
options.color_area = 'black';%[128 193 219]./255;    % Blue theme
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(photo_minus_drawing,options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(photo_minus_sketch,options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(drawing_minus_sketch,options);
% plot stats 
sig_photo_minus_drawing(sig_photo_minus_drawing==0)= NaN; 
plot(options.x_axis,sig_photo_minus_drawing*-0.4, 'black');
sig_photo_minus_sketch(sig_photo_minus_sketch==0)= NaN; 
plot(options.x_axis,sig_photo_minus_sketch*-0.6, 'color',cmap(ceil(256),:));
sig_drawing_minus_sketch(sig_drawing_minus_sketch==0)= NaN; 
plot(options.x_axis,sig_drawing_minus_sketch*-0.8, 'color',cmap(ceil(200),:));
yline(0,':')
for i= 0:0.1:1
xline(i,':');
end
xlim([-0.1 1])
ylim([-2 10])
legend(this_line, 'Photo-Drawing', 'Photo-Sketch', 'Drawing-Sketch')
title('Category Decoding Differences - MEG')
ylabel({'Decoding Accuracy', 'Difference (%)'})
xlabel('Time (s)')

print(fullfile(figure_path, ['group_cat_decoding_differences.svg']), ...
              '-dsvg', '-r600')


%% load crossdecoding results

data_dir = fullfile(path,'data/meg/crossdecoding/');

load(fullfile(data_dir,'photo_drawing_group_results.mat'))
load(fullfile(data_dir,'drawing_sketch_group_results.mat'))
load(fullfile(data_dir,'photo_sketch_group_results.mat'))

%% compute stats for crossdecoding 

significance_th = 0.05;
nboot = 100000; 
tail = 'right'; 

% set rng to a fixed number 
rng(96);

[~,~,~,~,crossdecoding_photo_drawing_clusters,~,~,crossdecoding_photo_drawing_cluster_size, crossdecoding_photo_drawing_cluster_th] = permutation_cluster_1sample_weight_alld (photo_drawing_group_acc-50, nperm, cluster_th, significance_th, tail);

[~,~,~,~,crossdecoding_drawing_sketch_clusters,~,~,crossdecoding_drawing_sketch_cluster_size, crossdecoding_drawing_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (drawing_sketch_group_acc-50, nperm, cluster_th, significance_th, tail);

[~,~,~,~,crossdecoding_photo_sketch_clusters,~,~,crossdecoding_photo_sketch_cluster_size, crossdecoding_photo_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (photo_sketch_group_acc-50, nperm, cluster_th, significance_th, tail);

%get maximum threshold for correcting for multiple comparisons
max_thr = max([crossdecoding_photo_drawing_cluster_th, crossdecoding_photo_sketch_cluster_th, crossdecoding_drawing_sketch_cluster_th]); 

% intialize significance arrays 
sig_crossdecoding_photo_drawing = zeros(1,size(photo_drawing_group_acc,2));
sig_crossdecoding_photo_sketch = zeros(1,size(photo_sketch_group_acc,2));
sig_crossdecoding_drawing_sketch = zeros(1,size(drawing_sketch_group_acc,2));

%now threshold the clusters again based on the maximum threshold 
sig_crossdecoding_photo_drawing([crossdecoding_photo_drawing_clusters{crossdecoding_photo_drawing_cluster_size>max_thr}]) = 1;
sig_crossdecoding_photo_sketch([crossdecoding_photo_sketch_clusters{crossdecoding_photo_sketch_cluster_size>max_thr}]) = 1;
sig_crossdecoding_drawing_sketch([crossdecoding_drawing_sketch_clusters{crossdecoding_drawing_sketch_cluster_size>max_thr}]) = 1;

% bootstrap to get CIs for onset and peak 

photo_drawing_decoding_boot = bootstrap_fixed_1D(photo_drawing_group_acc-50, [-100:10:1000],nboot); 

drawing_sketch_decoding_boot = bootstrap_fixed_1D(photo_sketch_group_acc-50, [-100:10:1000],nboot); 

photo_sketch_decoding_boot = bootstrap_fixed_1D(drawing_sketch_group_acc-50, [-100:10:1000],nboot); 

% bootstrap the difference for comparison 

photo_drawing_vs_photo_sketch_bootdiff = bootstrap_fixed_1D_diff(photo_drawing_group_acc-50, photo_sketch_group_acc-50, [-100:10:1000],nboot); 

photo_drawing_vs_drawing_sketch_bootdiff = bootstrap_fixed_1D_diff(photo_drawing_group_acc-50,drawing_sketch_group_acc-50, [-100:10:1000],nboot); 

drawing_sketch_vs_photo_sketch_bootdiff = bootstrap_fixed_1D_diff(drawing_sketch_group_acc-50,photo_sketch_group_acc-50, [-100:10:1000],nboot); 

%% TOST for crossdecoding peaks  

%first get peaks 
for sub = 1:size(photo_drawing_group_acc,1) 
    
    [~, photo_drawing_peaks(sub)] = max(photo_drawing_group_acc(sub,:)); 
    
    [~, drawing_sketch_peaks(sub)] = max(drawing_sketch_group_acc(sub,:)); 

    [~, photo_sketch_peaks(sub)] = max(photo_sketch_group_acc(sub,:)); 

end

[photo_drawing_drawing_sketch_equ_p, stat] = tost('one_sample', [-1 1], photo_drawing_peaks-drawing_sketch_peaks);

[photo_drawing_photo_sketch_equ_p, stat] = tost('one_sample', [-1 1], photo_drawing_peaks-photo_sketch_peaks);

[drawing_sketch_photo_sketch_equ_p, stat] = tost('one_sample', [-1 1], drawing_sketch_peaks-photo_sketch_peaks);

[~,~,~, adj_equ_p] = fdr_bh([photo_drawing_drawing_sketch_equ_p photo_drawing_photo_sketch_equ_p drawing_sketch_photo_sketch_equ_p]);

%% plot crossdecoding with error bars 

clear this_line
fig = figure;
options.handle = fig;
options.x_axis = linspace(-.1,1,size(photo_drawing_group_acc,2));
options.error = 'sem';
options.color_area = 'black';%[128 193 219]./255;    % Blue theme
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(photo_drawing_group_acc,options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(photo_sketch_group_acc,options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(drawing_sketch_group_acc,options);
% plot stats 
sig_crossdecoding_photo_drawing(sig_crossdecoding_photo_drawing==0)= NaN; 
plot(options.x_axis,sig_crossdecoding_photo_drawing*49.3, 'black');
sig_crossdecoding_photo_sketch(sig_crossdecoding_photo_sketch==0)= NaN; 
plot(options.x_axis,sig_crossdecoding_photo_sketch*48.6, 'color',cmap(ceil(256),:));
sig_crossdecoding_drawing_sketch(sig_crossdecoding_drawing_sketch==0)= NaN; 
plot(options.x_axis,sig_crossdecoding_drawing_sketch*47.9, 'color',cmap(ceil(200),:));
yline(50,':');
for i= 0:0.1:1
xline(i,':');
end
xlim([-0.1 1]);
ylim([47 80])
legend(this_line, 'Photo-Drawing', 'Photo-Sketch','Drawing-Sketch')
title('Category Cross-Decoding')
ylabel('Decoding Accuracy (%)')
xlabel('Time (s)')

print(fullfile(figure_path, ['crossdecoding_group_level.svg']),...
             '-dsvg', '-r600')


%% load temporal generalization results

data_dir = fullfile(path,'data/meg/temporal_generalization/');

load(fullfile(data_dir,'photo_temp_gen_group_results.mat'))
load(fullfile(data_dir,'drawing_temp_gen_group_results.mat'))
load(fullfile(data_dir,'sketch_temp_gen_group_results.mat'))

%% compute stats for temporal generalization 

nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'right';

% set rng to a fixed number 
rng(96);

[~,sig_temp_gen_photo] = permutation_cluster_1sample_weight_alld (permute(temp_gen_photo_group_acc-50,[3 1 2]), nperm, cluster_th, significance_th, tail);
[~,sig_temp_gen_drawing] = permutation_cluster_1sample_weight_alld (permute(temp_gen_drawing_group_acc-50,[3 1 2]), nperm, cluster_th, significance_th, tail);
[~,sig_temp_gen_sketch] = permutation_cluster_1sample_weight_alld (permute(temp_gen_sketch_group_acc-50,[3 1 2]), nperm, cluster_th, significance_th, tail);

%% plot temporal generalization with fieldtrip and statistical masks

cmap = colormap('redblueTecplot');
close all
cmap = cmap.^2;
figure('units','normalized','outerposition',[0 0 1 1])
clim = [30 70];
subplot(1,3,1)
ft_plot_matrix(flipud(mean(temp_gen_photo_group_acc,3)), 'clim', clim, 'highlight', flipud(1-sig_temp_gen_photo), 'highlightstyle', 'outline');
xticks(linspace(1,111,12))
xticklabels(linspace(-0.1,1,12))
%xtickformat('%.1f')
yticks(linspace(1,111,12))
yticklabels(linspace(1,-0.1,12))
ylabel('Train time (s)', 'FontSize', 20)
xlabel('Test time (s)', 'FontSize', 20)
title('Photo', 'FontSize',22)
colormap(cmap)
cl = colorbar;
cl.Label.String = 'Decoding Accuracy (%)';
cl.Label.FontSize = 16; 
axis square
hold on 
x = linspace(1,111);
y = linspace(111,1);
plot(x,y, 'black', 'LineWidth', 1);
hold off
subplot(1,3,2)
ft_plot_matrix(flipud(mean(temp_gen_drawing_group_acc,3)), 'clim', clim, 'highlight', flipud(1-sig_temp_gen_drawing), 'highlightstyle', 'outline');
xticks(linspace(1,111,12))
xticklabels(linspace(-0.1,1,12))
%xtickformat('%.1f')
yticks(linspace(1,111,12))
yticklabels(linspace(1,-0.1,12))
title('Drawing', 'FontSize',22)
colorbar
axis square
hold on 
x = linspace(1,111);
y = linspace(111,1);
plot(x,y, 'black', 'LineWidth', 1);
hold off
subplot(1,3,3)
ft_plot_matrix(flipud(mean(temp_gen_sketch_group_acc,3)), 'clim', clim, 'highlight', flipud(1-sig_temp_gen_sketch), 'highlightstyle', 'outline');
xticks(linspace(1,111,12))
xticklabels(linspace(-0.1,1,12))
%xtickformat('%.1f')
yticks(linspace(1,111,12))
yticklabels(linspace(1,-0.1,12))
title('Sketch', 'FontSize',22)
colorbar
axis square
hold on 
x = linspace(1,111);
y = linspace(111,1);
plot(x,y, 'black', 'LineWidth', 1);
hold off

print(fullfile(figure_path, ['group_temp_gen.svg']), ...
             '-dsvg', '-r600')

         
%% run statistics on differences between temporal generalization results  

tail= 'both';
significance_th = 0.05; 

% set rng to a fixed number 
rng(96);

[~,sigWei_temp_gen_photo_minus_drawing,~,~,temp_gen_photo_minus_drawing_clusters,temp_gen_photo_minus_drawing_cluster_size,temp_gen_photo_minus_drawing_cluster_th] = permutation_cluster_1sample_weight_alld (permute(temp_gen_photo_group_acc-temp_gen_drawing_group_acc,[3 1 2]), nperm, cluster_th, significance_th, tail);
[~,sigWei_temp_gen_drawing_minus_sketch,~,~,temp_gen_drawing_minus_sketch_clusters,temp_gen_drawing_minus_sketch_cluster_size,temp_gen_drawing_minus_sketch_cluster_th]= permutation_cluster_1sample_weight_alld (permute(temp_gen_drawing_group_acc-temp_gen_sketch_group_acc,[3 1 2]), nperm, cluster_th, significance_th, tail);
[~,sigWei_temp_gen_photo_minus_sketch,~,~,temp_gen_photo_minus_sketch_clusters,temp_gen_photo_minus_sketch_cluster_size,temp_gen_photo_minus_sketch_cluster_th] = permutation_cluster_1sample_weight_alld (permute(temp_gen_photo_group_acc-temp_gen_sketch_group_acc,[3 1 2]), nperm, cluster_th, significance_th, tail);

%get maximum threshold for correcting for multiple comparisons
max_thr = max([temp_gen_photo_minus_drawing_cluster_th, temp_gen_photo_minus_sketch_cluster_th, temp_gen_drawing_minus_sketch_cluster_th]); 

% intialize significance arrays 
sig_temp_gen_photo_minus_drawing = zeros(size(temp_gen_photo_group_acc(:,:,1)));
sig_temp_gen_photo_minus_sketch = zeros(size(temp_gen_photo_group_acc(:,:,1)));
sig_temp_gen_drawing_minus_sketch = zeros(size(temp_gen_photo_group_acc(:,:,1)));

%now threshold the clusters again based on the maximum threshold 
sig_temp_gen_photo_minus_drawing([temp_gen_photo_minus_drawing_clusters{temp_gen_photo_minus_drawing_cluster_size>max_thr}]) = 1;
sig_temp_gen_photo_minus_sketch([temp_gen_photo_minus_sketch_clusters{temp_gen_photo_minus_sketch_cluster_size>max_thr}]) = 1;
sig_temp_gen_drawing_minus_sketch([temp_gen_drawing_minus_sketch_clusters{temp_gen_drawing_minus_sketch_cluster_size>max_thr}]) = 1;

%% plot the difference between temporal generalization results

clim = [-5 5];
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
ft_plot_matrix(flipud(mean(temp_gen_photo_group_acc,3)-mean(temp_gen_drawing_group_acc,3)), 'clim', clim, 'highlight', flipud(1-sig_temp_gen_photo_minus_drawing), 'highlightstyle', 'outline');
xticks(linspace(1,111,12))
xticklabels(linspace(-0.1,1,12))
%xtickformat('%.1f')
yticks(linspace(1,111,12))
yticklabels(linspace(1,-0.1,12))
ylabel('Train time (s)', 'FontSize', 20)
xlabel('Test time (s)', 'FontSize', 20)
colormap(cmap)
cl = colorbar;
cl.Label.String = 'Decoding Accuracy Difference (%)';
cl.Label.FontSize = 16; 
title('Photo-Drawing','FontSize',22)
x = linspace(1,111);
y = linspace(111,1);
hold on 
plot(x,y, 'black', 'LineWidth', 1);
hold off
axis square
subplot(1,3,2)
ft_plot_matrix(flipud(mean(temp_gen_photo_group_acc,3)-mean(temp_gen_sketch_group_acc,3)), 'clim', clim, 'highlight', flipud(1-sig_temp_gen_photo_minus_sketch), 'highlightstyle', 'outline');
xticks(linspace(1,111,12))
xticklabels(linspace(-0.1,1,12))
%xtickformat('%.1f')
yticks(linspace(1,111,12))
yticklabels(linspace(1,-0.1,12))
title('Photo-Sketch','FontSize',22)
colorbar
axis square
hold on 
x = linspace(1,111);
y = linspace(111,1);
plot(x,y, 'black', 'LineWidth', 1);
subplot(1,3,3)
ft_plot_matrix(flipud(mean(temp_gen_drawing_group_acc,3)-mean(temp_gen_sketch_group_acc,3)), 'clim', clim, 'highlight', flipud(1-sig_temp_gen_drawing_minus_sketch), 'highlightstyle', 'outline');
xticks(linspace(1,111,12))
xticklabels(linspace(-0.1,1,12))
%xtickformat('%.1f')
yticks(linspace(1,111,12))
yticklabels(linspace(1,-0.1,12))
title('Drawing-Sketch','FontSize',22)
colorbar
hold on 
x = linspace(1,111);
y = linspace(111,1);
plot(x,y, 'black', 'LineWidth', 1);
axis square

print(fullfile(figure_path, ['group_temp_gen_diff.svg']), ...
             '-dsvg', '-r600')
