% Rishabh Jain
% Sorting RFs


% clc;
% clear all;

filters= 'mixed';

neurons= 2500;
MV_counter= 250000;
r= 1;
max_radius_val= 35;
sLR= 0.01;
percent_val= 0.12;

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

% sort by intensity first
weight_vectors_min= min(WT, [], 2);
weight_vectors_max= max(WT, [], 2);

wts_range= weight_vectors_max- weight_vectors_min;

WT_min_subtracted= bsxfun(@minus,  WT, weight_vectors_min);
WTs_normed= bsxfun(@rdivide,  WT_min_subtracted, wts_range);
WTs_normed= WTs_normed.*255;

seed= 1;
WTs_sorted= pairwise_image_sorting(neurons, seed, WT, WTs_normed, feature_dimX);
figure(2); WTs_visualize(neurons, WTs_sorted, feature_dimX);






% INTERACTIVE SORTING LOOP....
WTs_normed= WTs_sorted;

clicks= 1;
while (clicks<100)
clicks= clicks +1;

[x y]= ginput(2);
xc= ceil(x/feature_dimX);
yc= ceil(y/feature_dimY);
[seed_n]= sub2ind([size(GRID)], xc, yc);



% pick the best you think and swap
cell_n= seed_n(1);
replacement= seed_n(2);

temp= WTs_sorted(cell_n, :);
WTs_sorted(cell_n, :)= WTs_sorted(replacement, :);
WTs_sorted(replacement, :)= temp;


figure(2); WTs_visualize(neurons, WTs_sorted, feature_dimX);

% RUN again after swap ....
WT= WTs_sorted;
WTs_normed= WTs_sorted;

WTs_sorted= pairwise_image_sorting(neurons, seed_n(1), WT, WTs_normed, feature_dimX);
figure(3); WTs_visualize(neurons, WTs_sorted, feature_dimX);


end;
 










