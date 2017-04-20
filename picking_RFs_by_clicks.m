% Rishabh Jain
% Picking RFs from a large set 


% clc;
% clear all;
close all;


% set parameters
filters= 'mixed';

neurons= 2500;
MV_counter= 1000000;
r= 1;
max_radius_val= 35;
sLR= 0.015;
percent_val= 0.07;

grid_OUTX= round(sqrt(neurons));
grid_OUTY= round(sqrt(neurons));

sensory_types= 1;
feature_dimX=  7+ r*2;
feature_dimY=  7+ r*2;
ip_dimen= [feature_dimX feature_dimY];

GRID= 1:grid_OUTX*grid_OUTY;
GRID= reshape(GRID, [grid_OUTX grid_OUTY]);

% file_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/development_with_NI_images/data/backup-data/');
% file_pre= strcat('NI-wts-frame=',num2str(MV_counter,'%06d'), ...;
%     '-neurons=', num2str(neurons, '%d'), ...
%     '-nbd_radius_val=', num2str(r, '%d'), ...
%     '-max_radius_val=', num2str(max_radius_val, '%d'), ...
%     '-sLR=',num2str(sLR, '%1.3f'), ...
%     '-percent=', num2str(percent_val, '%1.2f'), ...
%     '-filters_', filters);
% fOut= strcat(file_path, file_pre,'.mat');
% load(fOut);

% flatten the WT vector
sz_local_WTs= size(WT);

% normalize the range
weight_vectors_min= min(WT, [], 2);
weight_vectors_max= max(WT, [], 2);

wts_range= weight_vectors_max- weight_vectors_min;

WT_min_subtracted= bsxfun(@minus,  WT, weight_vectors_min);
WTs_normed= bsxfun(@rdivide,  WT_min_subtracted, wts_range);
WTs_normed= WTs_normed.*255;

seed= 1;
WTs_sorted= pairwise_image_sorting(neurons, seed, WT, WTs_normed, feature_dimX);
figure(1); WTs_visualize(neurons, WTs_sorted, feature_dimX);



%% PICK YOUR FAVORITES ...
how_many= 4;

for v=1:how_many
    zoom on;
    pause();
    zoom off;
    [x(v), y(v)]= ginput(1);
    display(v);
    %zoom out;
end;


xc= ceil(x/feature_dimX);
yc= ceil(y/feature_dimY);
[seed_n]= sub2ind([size(GRID)], xc, yc);

WTs_picked=[];
WTs_picked(1:how_many, :)= WTs_sorted(seed_n, :);

figure(2); WTs_visualize(how_many, WTs_picked, feature_dimX);
%figure(2); WTs_visualize_linear(how_many, 1, WTs_picked, feature_dimX);
print('picked_shape_RFs_categories_5.eps', '-depsc');
%export_fig(strcat('picked_RFs.eps'), '-transparent', '-m2');

