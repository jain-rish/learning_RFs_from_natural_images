% Rishabh Jain
% LNC, USC
% Feb 2015

clc;
clear all;
close all;

%% General parameters initialization
addpath('/amnt/foam/foamd0/rishabh/SOM_simulations/August_2011_Simulations/export_fig/');

max_radius_val= 1;
neurons= 5200;
supplied_parameters_natural_inputs_reduced_v4_size;

filters= 'mixed';
map= 'hybrid';


switch filters
    case 'mixed'
        sensory_types=  3;
    case 'curved'
        sensory_types=  2;
    case 'straight'
        sensory_types = 1;
end


d=  2*max_radius_val+1;
rf_half= floor(gabor_size/2);

% neurons of interest
n_ids= 36-max_radius_val: 36+max_radius_val;
xc= 36;
yc= 36;




%% General parameters and file reading

knee= 40;
offset= 0.6;
s_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/spring_2015_data/hybrid-data-awesome/single-scale/');
%s_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/spring_2015_data/spatial-data-awesome/');
s_pre= strcat('v1-neighborhood-activity', ...
    '-max_radius_val_', num2str(max_radius_val, '%d'), ...
    '-offset_', num2str(offset, '%1.2f'), ...
    '-knee_', num2str(knee, '%d'),'-substrate_', filters, '_norm_', map, '_NNlinear');


fOut= strcat(s_path, s_pre,'.mat');
load(fOut);


data_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/spring_2015_data/hybrid-data-awesome/single-scale/');
%data_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/spring_2015_data/spatial-data-awesome/');
data_pre= strcat('WTs_v1-filters', ...
              '-max_radius_val_', num2str(max_radius_val, '%d'), ...
              '-filters_', filters, ...
              '-map_', map);
fdata= strcat(data_path, data_pre,'.mat');
load(fdata);




%% Natural image read-out from the database

% read the directory information
dirinfo= dir('/amnt/wave/waved1/image_data/corel/');
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
subdirinfo = cell(length(dirinfo));
for K = 1:length(dirinfo)
    thisdir = dirinfo(K).name;
    subdirinfo{K} = dir(fullfile('/amnt/wave/waved1/image_data/corel/', thisdir, '*.ppm'));
end

%len = uint16(round((patch_size-1)/2));
len = single(round((patch_size-1)/2));



% Matrix the top n-stimuli
count= 100;


%D= v1nbrhood_actv(xs, ys, :);
D= v1nbrhood_actv(max_radius_val+1, max_radius_val+1, :); % get the center of the box ...

[sortvals, sortids] = sort(squeeze(D), 'descend');

idx_start= 0;

top_folder=     img_folder(sortids(idx_start+ 1: idx_start+ count));
top_subfolder=  img_subfolder(sortids(idx_start+ 1: idx_start+ count));
top_positions=  img_positions(sortids(idx_start+ 1: idx_start+ count), :);


% scores for display
scores= sortvals(1:100);
scores= reshape(scores, 10, 10);
%scores= scores';





%% FIGURE 1 (Only red boxes)

% normalize the visual coordinates
% +1 instead of +visV1_left since the image is already cut-off

MxI= 1;
MyI= 1;
MxA= V1_grid_OUTX;
MyA= V1_grid_OUTY;
% x = (xc - MxI)/(MxA - MxI)*(visV1right- visV1left)+ 1;
% y = (yc - MyI)/(MyA - MyI)*(visV1right- visV1left)+ 1;
% img_xc= round(x);
% img_yc= round(y);
q= 33;
img_xc= floor(q/2);
img_yc= floor(q/2);

[top_image_patches bigImage]= image_displayer(max_radius_val, img_xc, img_yc, top_folder, top_subfolder, top_positions, dirinfo, subdirinfo, len);


nRows = 10;
nCols = 10;
%q= visV1right- visV1left+1;
q= (d+ gabor_size -1) +2;

figure;
imshow(bigImage, [0 255]); hold on;
for row = 1:nRows
    for colmn= 1:nCols
        rectangle('Position', [(colmn-1)*q+ 2-0.5...
                               (row-1)*q+ 2-0.5...
                                d+2*rf_half d+2*rf_half], 'EdgeColor', 'r', 'LineWidth', 0.1);
    end;
end;

    
% save the file
s_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/spring_2015_data/hybrid-data-awesome/single-scale/');
s_pre= strcat('NI', ...
              '-max_radius_val_', num2str(max_radius_val, '%d'), ...
              '-substrate_', filters);
          

% produce a CENTERED figure
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');

papersize = get(gcf, 'PaperSize');
width =  16;         % Initialize a variable for width.
height = 16;          % Initialize a variable for height
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

export_fig(strcat(s_path, s_pre, '-hybrid-coverage','.eps'), '-transparent', '-m2');






%% FIGURE 2 (red-boxes with context)
q= 33;
img_xc= floor(q/2);
img_yc= floor(q/2);

x = (xc - MxI)/(MxA - MxI)*(33- 1)+ 1;
y = (yc - MyI)/(MyA - MyI)*(33- 1)+ 1;

[top_image_patches bigImage]= image_displayer_larger_context(img_xc, img_yc, top_folder, top_subfolder, top_positions, dirinfo, subdirinfo, len);

figure;
imshow(bigImage, [0 255]); hold on;
for row = 1:nRows
    for colmn= 1:nCols
        rectangle('Position', [(colmn-1)*q+ round(y)-1.5- (max_radius_val)- rf_half...
                               (row-1)*q+   round(x)-1.5- (max_radius_val)- rf_half...
                               d+2*rf_half d+2*rf_half], 'EdgeColor', 'r', 'LineWidth', 0.4);
    end;
end;


% save the file
s_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/spring_2015_data/hybrid-data-awesome/single-scale/');
s_pre= strcat('NI', ...
              '-max_radius_val_', num2str(max_radius_val, '%d'), ...
              '-substrate_', filters);
          

% produce a CENTERED figure
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');

papersize = get(gcf, 'PaperSize');
width =  16;         % Initialize a variable for width.
height = 16;          % Initialize a variable for height
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

export_fig(strcat(s_path, s_pre, '-hybrid-larger-context','.eps'), '-transparent', '-m2');






%% Set general parameters depending on curve/straight and hybrid/spatial

figure;
imshow(bigImage, [0 255]); hold on;


% version 3
viscn= (visV1right-visV1left)/2 +1;
[xo_v yo_v]= meshgrid( linspace(viscn -max_radius_val, viscn +max_radius_val, d), linspace(viscn -max_radius_val, viscn +max_radius_val, d));
xo_layer= zeros(V1_grid_OUTX, V1_grid_OUTY); xo_layer(n_ids, n_ids)= xo_v; 
yo_layer= zeros(V1_grid_OUTX, V1_grid_OUTY); yo_layer(n_ids, n_ids)= yo_v;




% make visualization for each NI 
for row = 1:nRows
    for colmn= 1:nCols
        
        % get individual images from the BIG matrix.
        single_img= bigImage((colmn-1)*q+1:(colmn-1)*q +q, (row-1)*q+1:(row-1)*q +q);
 
        switch filters
            case 'mixed'
                oriens_list= 0: 22.5: 359;
                V1_image= get_v1_activitymap_mix_filters_NORM_only_sines(WTs_V1, sensory_types, max_radius_val, single_img, xo_layer, yo_layer, oriens_list, gabor_size, map);                
            case 'curved'
                oriens_list= 0: 22.5: 359;
                V1_image= get_v1_activitymap_curved_filters_NORM_only_sines(WTs_V1, sensory_types, max_radius_val, single_img, xo_layer, yo_layer, oriens_list, gabor_size, map);
            case 'straight'
                % for computational efficiency
                oriens_list= 0: 22.5: 359;
                V1_image= get_v1_activitymap_straight_filters_NORM_only_sines(WTs_V1, sensory_types, max_radius_val, single_img, xo_layer, yo_layer, oriens_list, gabor_size, map);      
        end
        
        
        %f-I curve
        % V1_image = 2.2*( V1_image- offset)./(1- exp(-knee*(V1_image- offset)));
        
        % further processing in the firing-rate display code
        matrix_info(1)= row;
        matrix_info(2)= colmn;        
        
        firing_rates_on_natural_images_feb_2015(matrix_info, filters, map, V1_image, sensory_types, xo_layer, yo_layer);
        
        % make the red-circle
        % t= 0:.001:2*pi;
        % plot( (colmn-1)*q+   round(y)+rad*sin(t), ...
        %       (row-1)*q+     round(x)+rad*cos(t), 'r', 'LineWidth', 0.25);
        rectangle('Position', [(row-1)*q+     round(y)-max_radius_val-0.5   ...
                               (colmn-1)*q+   round(x)-max_radius_val-0.5 d d], 'EdgeColor', 'r', 'LineWidth', 1);
        
        text( (row)*q-3, (colmn-1)*q+2, num2str(scores(colmn, row), '%1.1f'), 'FontSize', 10) ; % r is the radius
        
    end;
end;

%%

hold off;

% produce a CENTERED figure
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');

papersize = get(gcf, 'PaperSize');
width =  16;         % Initialize a variable for width.
height = 16;          % Initialize a variable for height
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

s_pre= strcat('v1-top100-images', ...
    '-max_radius_val_', num2str(max_radius_val, '%d'), ...
    '-offset_', num2str(offset, '%1.2f'), ...
    '-knee_', num2str(knee, '%d'),'-substrate_', filters, '_norm_', map, '_NNlinear');
%export_fig(strcat(s_path, s_pre,'.eps'), '-transparent', '-m2');

h= gcf;
saveas(h, strcat(s_path, s_pre, '.fig'),'fig');
%openfig( strcat(s_path, s_pre,'.fig'), 'new', 'visible');



% figure;
% plot((idx_start+1: idx_start+count), sortvals(idx_start+1: idx_start+count), 'r.-');
% ylabel('Total neighborhood-Activity','fontname','Arial','fontsize',20);
% xlabel('Image Index','fontname','Arial','fontsize',20);




