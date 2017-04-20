% Rishabh Jain
% LNC, USC
% Feb 2015

function NI_pick_from_V1_filters_options_NORM_sines_hybrid(neurons, max_radius_val, knee, filters)

supplied_parameters_natural_inputs_reduced_v4_size;
warning('off', 'all');


%% Built V1 map

clear saved_data;
map= 'hybrid';

switch filters    
    case 'mixed'
        sensory_types=  3;
        strght_rf_AR=   2; 
        curved_rf_AR_1= 0.12;  
        curved_r1=      0.16;
        offset_curved=  0.42;  % chosen after calibrating with straights
        offset_straights=  0.60;          
        
    case 'curved'
        sensory_types=  2;
        curved_rf_AR_1= 0.12; 
        curved_r1=      0.16;
        offset_curved=  0.42; 
        
    case 'straight'
        sensory_types=  1;
        strght_rf_AR=   2;        
        offset_straights=  0.60;  
end

% Orientation steps...
oriens_list= 0: 22.5: 359;


% neurons of interest
n_ids= 36-max_radius_val: 36+max_radius_val;
d= 2*max_radius_val+1;


% Controlled RF center layout (10th Nov. 2014)
% version 3
viscn= visV1left+ (visV1right-visV1left)/2 +1;
[yo_v xo_v]= meshgrid( linspace(viscn -max_radius_val, viscn +max_radius_val, d), linspace(viscn -max_radius_val, viscn +max_radius_val, d));
xo_layer= zeros(V1_grid_OUTX, V1_grid_OUTY); xo_layer(n_ids, n_ids)= xo_v; 
yo_layer= zeros(V1_grid_OUTX, V1_grid_OUTY); yo_layer(n_ids, n_ids)= yo_v;

% calculate the RF offsets to "accurately" create gabors 
offsetx= ceil(xo_layer)- xo_layer;
offsety= ceil(yo_layer)- yo_layer;



%% Initialization of the V1 layer

WTs_V1=     zeros(gabor_size, gabor_size, d, d, 8, 1, 'single');


% V1 layer
% Initialize with GABOR map
for x=1:1:d
    for y=1:1:d
        
        OUTX= 36-(max_radius_val+1)+ x;
        OUTY= 36-(max_radius_val+1)+ y;
                
        % CURVED
        if (strcmp('curved', filters))
        
        %oriens_list= 0: 22.5: 359;      
        for or= 1:length(oriens_list)            
        sin_gabor=  curved_sin_gabor_dark_out(gabor_size, 2, oriens_list(or), curved_rf_AR_1, curved_r1, offsetx(OUTX, OUTY), offsety(OUTX, OUTY));
        % normalized gabors
        sin_gabor_0_mean = sin_gabor - mean(sin_gabor(:));
        sqrt_sin_gabor   = sqrt((sin_gabor_0_mean(:))'*sin_gabor_0_mean(:));
        sin_gabor_norm   = sin_gabor_0_mean ./sqrt_sin_gabor;
        
        WTs_V1(:, :, x, y, or, 1) = sin_gabor_norm;
        end;
        
        for or= 1:length(oriens_list)  
        sin_gabor=  curved_sin_gabor_light_out(gabor_size, 2, oriens_list(or), curved_rf_AR_1, curved_r1, offsetx(OUTX, OUTY), offsety(OUTX, OUTY));  
        % normalized gabors
        sin_gabor_0_mean = sin_gabor - mean(sin_gabor(:));
        sqrt_sin_gabor   = sqrt((sin_gabor_0_mean(:))'*sin_gabor_0_mean(:));
        sin_gabor_norm   = sin_gabor_0_mean ./sqrt_sin_gabor;
        
        WTs_V1(:, :, x, y, or, 2) = sin_gabor_norm;
        end;        
           
        
        
        % STRAIGHT
        elseif (strcmp('straight', filters))
            
        %oriens_list= 0: 22.5: 179;            
        for or= 1:length(oriens_list)  
        sin_gabor= make_sin_gabor_adjust_rfc(gabor_size, 1, oriens_list(or) -90, strght_rf_AR, offsetx(OUTX, OUTY), offsety(OUTX, OUTY));   
        % normalized gabors
        sin_gabor_0_mean = sin_gabor - mean(sin_gabor(:));
        sqrt_sin_gabor   = sqrt((sin_gabor_0_mean(:))'*sin_gabor_0_mean(:));
        sin_gabor_norm   = sin_gabor_0_mean ./sqrt_sin_gabor;
        
        WTs_V1(:, :, x, y, or, 1) = sin_gabor_norm;
        end;
        
        
        
        
        % MIX- Both straight and curved ...
        
        elseif (strcmp('mixed', filters))
        %oriens_list= 0: 22.5: 359;        
        for or= 1:length(oriens_list)            
        % dark_out
        sin_gabor=  curved_sin_gabor_dark_out(gabor_size, 2, oriens_list(or), curved_rf_AR_1, curved_r1, offsetx(OUTX, OUTY), offsety(OUTX, OUTY));
        % normalized gabors
        sin_gabor_0_mean = sin_gabor - mean(sin_gabor(:));
        sqrt_sin_gabor   = sqrt((sin_gabor_0_mean(:))'*sin_gabor_0_mean(:));
        sin_gabor_norm   = sin_gabor_0_mean ./sqrt_sin_gabor;
        
        WTs_V1(:, :, x, y, or, 1) = sin_gabor_norm;
        end;
        
        for or= 1:length(oriens_list)  
        % light_out 
        sin_gabor=  curved_sin_gabor_light_out(gabor_size, 2, oriens_list(or), curved_rf_AR_1, curved_r1, offsetx(OUTX, OUTY), offsety(OUTX, OUTY));  
        % normalized gabors
        sin_gabor_0_mean = sin_gabor - mean(sin_gabor(:));
        sqrt_sin_gabor   = sqrt((sin_gabor_0_mean(:))'*sin_gabor_0_mean(:));
        sin_gabor_norm   = sin_gabor_0_mean ./sqrt_sin_gabor;
        
        WTs_V1(:, :, x, y, or, 2) = sin_gabor_norm;
        end;        
                    
                                   
        for or= 1:length(oriens_list)  
        sin_gabor= make_sin_gabor_adjust_rfc(gabor_size, 1, oriens_list(or) -90, strght_rf_AR, offsetx(OUTX, OUTY), offsety(OUTX, OUTY));   
        % normalized gabors
        sin_gabor_0_mean = sin_gabor - mean(sin_gabor(:));
        sqrt_sin_gabor   = sqrt((sin_gabor_0_mean(:))'*sin_gabor_0_mean(:));
        sin_gabor_norm   = sin_gabor_0_mean ./sqrt_sin_gabor;
        
        WTs_V1(:, :, x, y, or, 3) = sin_gabor_norm;
        end;
        
        
 
        end;    
                        
    end;
end;




%% MAIN ALGORITHM .... 

% read the directory information
dirinfo= dir('/amnt/wave/waved1/image_data/corel/');
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
subdirinfo = cell(length(dirinfo));
for K = 1:length(dirinfo)
    thisdir = dirinfo(K).name;
    subdirinfo{K} = dir(fullfile('/amnt/wave/waved1/image_data/corel/', thisdir, '*.ppm'));
end


MV_max= 120000; %length(ppmfiles)*60;
% MV_counter= 0; %length(ppmfiles)*60;

repeats_per_image= 120;
len = uint16(round((patch_size-1)/2));


% % initialize some variables
v1nbrhood_actv=  zeros(d, d, 1, MV_max, 'single');

% individual cell firing rates (only for calibration)
% v1cell_actv=  zeros(d, d, 16, sensory_types, MV_max, 'single');


img_folder= zeros(MV_max, 1, 'single');
img_subfolder= zeros(MV_max, 1, 'single');
img_positions= zeros(MV_max, 3, 'single');



MV_counter= 1;

for repeat_image= 1:repeats_per_image
    
    while (MV_counter <=MV_max)
        
        folder= randsample(length(dirinfo) -2, 1, 'true');
        folder= folder +2; % avoid the . and ..
        subfolder= randsample(length(subdirinfo{folder}), 1, 'true');
        
        [Irotated small_d]= image_read(folder, subfolder, dirinfo, subdirinfo);
        
        % step 2: pick a random patch from the image set
        randx= randi([len+1 small_d/4-len-1], repeats_per_image,1);
        randy= randi([len+1 small_d/4-len-1], repeats_per_image,1);
        randz= randi([1 15], repeats_per_image, 1);
        
        image_patch(:,:) = Irotated(randx(repeat_image)-len:randx(repeat_image)+len, ...
                                    randy(repeat_image)-len:randy(repeat_image)+len, ...
                                    randz(repeat_image));
        
        % normalize (7th Oct. 2014)
        image_patch= image_patch./255;
        
        
        %% correlation of V1 RFs with natural image patches
        
        if (strcmp('curved', filters))
        V1_image= get_v1_activitymap_curved_filters_NORM_only_sines(WTs_V1, sensory_types, max_radius_val, image_patch, xo_layer, yo_layer, oriens_list, gabor_size, map);          
        V1_image = (2.2).*( V1_image- offset_curved)./(1- exp(-knee*(V1_image- offset_curved)));
        
        elseif (strcmp('straight', filters))
        V1_image= get_v1_activitymap_straight_filters_NORM_only_sines(WTs_V1, sensory_types, max_radius_val, image_patch, xo_layer, yo_layer, oriens_list, gabor_size, map);  
        V1_image = (2.2).*( V1_image- offset_straights)./(1- exp(-knee*(V1_image- offset_straights)));
        
        elseif (strcmp('mixed', filters))
        V1_image= get_v1_activitymap_mix_filters_NORM_only_sines(WTs_V1, sensory_types, max_radius_val, image_patch, xo_layer, yo_layer, oriens_list, gabor_size, map);         
        % f-I curve  ( different thresholds for curved and straights ) 
        V1_image(:,:,:,1:2) = (2.2).*( V1_image(:,:,:,1:2)- offset_curved)./(1- exp(-knee*(V1_image(:,:,:,1:2)- offset_curved)));
        V1_image(:,:,:,  3) = (2.2).*( V1_image(:,:,:,  3)- offset_straights)./(1- exp(-knee*(V1_image(:,:,:,3)- offset_straights)));
        end;
        

        %% f-I curve        
                
        % V1_image = (2.2).*( V1_image- offset)./(1- exp(-knee*(V1_image- offset)));
        

        %% maximal local energy calculation        
        % maximum responsive channel
        if (strcmp('curved', filters) || strcmp('mixed', filters))
        maxV1_image1 = max(V1_image, [], 4);
        % maximum responsive orientation
        maxV1_image = max(maxV1_image1, [], 3);
        
        
        else        
        maxV1_image = max(V1_image, [], 3);
        end;
        
        % local nbd. SUM
        B = padarray(maxV1_image, [max_radius_val+1 max_radius_val+1]);
        s = cumsum(B,1);
        c = s(1+d:end-1,:)-s(1:end-d-1,:);
        s = cumsum(c,2);
        v1local_nbd_summed = s(:,1+d:end-1)-s(:,1:end-d-1);     

       
        %% store ...
        % individual cell firing rates (only for calibration)
        % v1cell_actv(:, :, :, :, MV_counter)=  V1_image;
        
        % population firing rates
        v1nbrhood_actv(:, :, MV_counter)=  v1local_nbd_summed;
        
        img_folder(MV_counter)=           folder;
        img_subfolder(MV_counter)=        subfolder;
        img_positions(MV_counter, :)=     [randx(repeat_image) randy(repeat_image) randz(repeat_image)];
        
        MV_counter= MV_counter+ 1;
        
    end;
    
end;


        
        
%%
% save the variables
% make the path
s_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/spring_2015_data/');
s_pre= strcat('v1-neighborhood-activity', ...
              '-max_radius_val_', num2str(max_radius_val, '%d'), ...
              '-knee_', num2str(knee, '%d'), ...
              '-substrate_', filters,'_norm_hybrid_NNlinear');
          
fOut= strcat(s_path, s_pre,'.mat');
save(fOut, 'v1nbrhood_actv', 'img_folder', 'img_subfolder', 'img_positions');



