%% Script for running all the first level analyses of the MEG data
% possible analysis steps are:
% Category Decoding 
% Category Crossdecoding 
% Temporal Generalization 
% Temporal Cross-Generalization (Training on one condition testing on the
% other) 
% RSA 

clear all
clc

%setup paths 

path = pwd;
figure_path = fullfile(path,'figures');

% add utils 

addpath(fullfile(path,'utils'));

% add first level functions

addpath(fullfile(path,'first_level','meg'));

% setup the decoding toolbox 
try 
    decoding_defaults;
catch 
    tdt_path = input('The Decoding Toolbox seems to be not on your path. Please enter the path to your TDT version:\n','s');
    addpath(tdt_path);
    decoding_defaults;
end 

% setup fieldtrip 
try 
    ft_defaults;
catch
    ft_path = input('fieldtrip seems to be not on your path. Please enter the path to your fieldtrip version:\n','s');
    addpath(ft_path);
    ft_defaults;
end

% set plot defaults 

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3, 'DefaultAxesFontName', 'Helvetica','DefaultTextFontName', 'Helvetica') 

% get colormap 
cmap = colormap('redblueTecplot');
close all

% specify decoding parameters 

n_cat = 48;
n_perm = 5; 
n_average = 2; 

% specify subs to exclude 

excluded_subs = {};

% specify which steps to compute 

cfg.do.decoding = 1; % runs the category decoding for a single subject for photos, drawings and sketches separately 
cfg.do.cross_decoding = 1; % runs the category crossdecoding for a single subject for three comparisons: photo-drawing, drawing-sketch and photo-sketch separately 
cfg.do.temp_gen = 1; % runs the temporal generalization analysis for a single subject for photos, drawings and sketches separately 
cfg.do.RSA = 1; % computes representational dissimilarity matrices (RDMs) for a single subject for each time point and for photos, drawings and sketches separately 

%% Category Decoding

if cfg.do.decoding 

sub_id = 'od12a'; % sample data is provided for the subject "od12a" and can be retrieved from OSF

%load data

load(fullfile(path,'data','meg','preproc', [sub_id,'_preproc_photo_data.mat']));
load(fullfile(path,'data','meg','preproc', [sub_id,'_preproc_drawing_data.mat']));
load(fullfile(path,'data','meg','preproc', [sub_id,'_preproc_sketch_data.mat']));

% do decoding
[photo_accs, drawing_accs, sketch_accs] = meg_cat_decoding(photo_data_all,drawing_data_all,sketch_data_all, n_cat, n_perm, n_average);

% create path for decoding results if not already there
if ~isdir(fullfile(path,'data','meg','decoding')), mkdir(fullfile(path,'data','meg','decoding')), end

%plot

figure 
plot(photo_accs)
hold on 
plot(drawing_accs)
plot(sketch_accs)
yline(50,':');
xline(11,':');
xticks(linspace(1,111,12))
xlim([1 111])
xticklabels(linspace(-0.1,1,12))
legend('Photos', 'Drawings', 'Sketches');
title(['Decoding for subject: ', sub_id]);
print(fullfile(path,'data','meg','decoding', [sub_id,'_decoding_final.jpg']),'-djpeg')

% save results  
save(fullfile(path,'data','meg','decoding', [sub_id,'_decoding_accs.mat']), 'photo_accs', 'drawing_accs', 'sketch_accs')

end 
%% Category Crossdecoding 

if cfg.do.cross_decoding

sub_id = 'od12a'; % sample data is provided for the subject "od12a" and can be retrieved from OSF 

%load data

load(fullfile(path,'data','meg','preproc', [sub_id,'_preproc_photo_data.mat']));
load(fullfile(path,'data','meg','preproc', [sub_id,'_preproc_drawing_data.mat']));
load(fullfile(path,'data','meg','preproc', [sub_id,'_preproc_sketch_data.mat']));

% do decoding
[photo_drawing_accs, drawing_photo_accs, photo_sketch_accs, sketch_photo_accs, drawing_sketch_accs, sketch_drawing_accs] = meg_cat_crossdecoding(photo_data_all,drawing_data_all,sketch_data_all, n_cat, n_perm, n_average);

% create path for decoding results if not already there
if ~isdir(fullfile(path,'data','meg','crossdecoding')), mkdir(fullfile(path,'data','meg','crossdecoding')), end

% plot
figure 
plot(photo_drawing_accs)
hold on 
plot(photo_sketch_accs)
plot(drawing_photo_accs)
plot(sketch_photo_accs)
plot(drawing_sketch_accs)
plot(sketch_drawing_accs)
yline(50,':');
xticks(linspace(1,111,12));
xlim([1 111]);
xticklabels(linspace(-0.1,1,12));
legend('Photo-Drawing', 'Photo-Sketch', 'Drawing-Photo','Sketch-Photo', 'Drawing-Sketch','Sketch-Drawing');
title(['Cross-Decoding for subject: ', sub_id]);
print(fullfile(path,'data','meg','crossdecoding', [sub_id,'_crossdecoding.jpg']),'-djpeg')
close all

% save results  
save(fullfile(path,'data','meg','crossdecoding', [sub_id,'_crossdecoding.mat']), 'photo_drawing_accs', 'drawing_photo_accs', 'photo_sketch_accs', 'sketch_photo_accs','drawing_sketch_accs','sketch_drawing_accs')

end 

%% Temporal Generalization 

if cfg.do.temp_gen 

sub_id = 'od12a'; % sample data is provided for the subject "od12a" and can be retrieved from OSF 

%load data

load(fullfile(path,'data','meg','preproc', [sub_id,'_preproc_photo_data.mat']));
load(fullfile(path,'data','meg','preproc', [sub_id,'_preproc_drawing_data.mat']));
load(fullfile(path,'data','meg','preproc', [sub_id,'_preproc_sketch_data.mat']));

% do decoding
[photo_accs, drawing_accs, sketch_accs] = meg_cat_temporal_generalization(photo_data_all,drawing_data_all,sketch_data_all, n_cat,n_perm, n_average);

% create path for decoding results if not already there
if ~isdir(fullfile(path,'data','meg','temporal_generalization')), mkdir(fullfile(path,'data','meg','temporal_generalization')), end

% plot

figure 
clim = [20 80];
imagesc(flipud(photo_accs), clim)
xticks(linspace(1,111,12))
xticklabels(linspace(-0.1,1,12))
yticks(linspace(1,111,12))
yticklabels(linspace(1,-0.1,12))
title('Temporal Generalization Photo')
print(fullfile(path,'data','meg','temporal_generalization', [sub_id,'_temp_gen_photo.jpg']),'-djpeg')
figure
imagesc(flipud(drawing_accs), clim)
xticks(linspace(1,111,12))
xticklabels(linspace(-0.1,1,12))
yticks(linspace(1,111,12))
yticklabels(linspace(1,-0.1,12))
title('Temporal Generalization Drawing')
print(fullfile(path,'data','meg','temporal_generalization', [sub_id,'_temp_gen_drawing.jpg']),'-djpeg')
figure
imagesc(flipud(sketch_accs), clim)
xticks(linspace(1,111,12))
xticklabels(linspace(-0.1,1,12))
yticks(linspace(1,111,12))
yticklabels(linspace(1,-0.1,12))
title('Temporal Generalization Sketch')
print(fullfile(path,'data','meg','temporal_generalization', [sub_id,'_temp_gen_sketch.jpg']),'-djpeg')

% save results  

save(fullfile(path,'data','meg','temporal_generalization', [sub_id,'_temp_gen_accs.mat']), 'photo_accs', 'drawing_accs', 'sketch_accs')


end 

%% RSA - pearson correlation as a distance measure, with noise normalization 

if cfg.do.RSA

sub_id = 'od12a'; % sample data is provided for the subject "od12a" and can be retrieved from OSF 

%load data

load(fullfile(path,'data','meg','preproc', [sub_id,'_preproc_photo_data.mat']));
load(fullfile(path,'data','meg','preproc', [sub_id,'_preproc_drawing_data.mat']));
load(fullfile(path,'data','meg','preproc', [sub_id,'_preproc_sketch_data.mat']));

% create path for decoding results if not already there
if ~isdir(fullfile(path,'data','meg','rsa')), mkdir(fullfile(path,'data','meg','rsa')), end

% do decoding
[photo_RDM, drawing_RDM, sketch_RDM] = meg_pearson_noisenorm(photo_data_all,drawing_data_all,sketch_data_all, n_cat);


%plot the RDMs for a single timepoint
for time = 20
 figure 
 subplot(1,3,1)
 imagesc(photo_RDM(:,:,time))
 subplot(1,3,2)
 imagesc(drawing_RDM(:,:,time))
 subplot(1,3,3)
 imagesc(sketch_RDM(:,:,time))
end 

% save results  
save(fullfile(path,'data','meg','rsa', [sub_id,'_pearson_noisenorm_RDMs.mat']), 'photo_RDM','drawing_RDM','sketch_RDM')


end 
