% Rishabh Jain
% LNC
% April 2015

% This version does not run 8 orientations of the image patch per trial as
% last WORKS_v2 version 2 does


function SOM_v69_hybrid_random_inputs_WORKS_v3(neurons, filters, r, max_radius_val, sLR, percent_val)
% clc;
% clear all;
% profile on;

grid_OUTX= round(sqrt(neurons));
grid_OUTY= round(sqrt(neurons));

feature_dimX=  7+ r*2;
feature_dimY=  7+ r*2;
ip_dimen= [feature_dimX feature_dimY];



%% General parameters and file reading

map= 'hybrid';
gabor_size =feature_dimX;

%filters= 'mixed';

% mixed parameters
strght_rf_AR=   2;
curved_rf_AR_1= 0.12;
curved_r1=      0.16;
offset_curved=  0.42;  % chosen after calibrating with straights
offset_straights=  0.60;




% file read ...
knee= 40;
data_read_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/development_with_NI_images/data/');
s_pre= strcat('v1-neighborhood-activity', ...
    '-max_radius_val_', num2str(r, '%d'), ...
    '-knee_', num2str(knee, '%d'),'-substrate_', filters, '_norm_', map, '_NNlinear');


fOut= strcat(data_read_path, s_pre,'.mat');
load(fOut);






%% Update parameters with supplied ones from UNIX script

sensory_types=  1;
MV_max= 150000; 
save_step= round(linspace(1, MV_max, 4));


sampling= round(MV_max/10);

%t_phase1= round(MV_max/100);
%pval(1:t_phase1)= 100; % phase 1
pval(1:sampling)= linspace(percent_val, 0.02, sampling); % Assuming 0.025% is 1 neuron


% (7th May 2015) Reversed the learning rate curve- High LR with more neighbors and low LR towards the end 
% LR_space= linspace(0.33, 1, 100);

LR_space= linspace(1, 1, 100); % 2-fold jump in 1 to 100%
LR_curve= ( linspace(LR_space(round(percent_val*150)), 1, sampling) ); % lower the percentile, higher the jump
LR(1:sampling)= sLR .*LR_curve; % scaling LR

total_nbrs= round((2*max_radius_val+1)*(2*max_radius_val+1) /sensory_types);

% radius= round(linspace(max_radius_val, end_radius_val, sampling));
warning off all;







%% Pre-allocation of arrays

%percent_training=  Inf([MV_max 1]);
%mean_pop_response= Inf([MV_max 1]);

% set the random streams based on your parameter
stream = RandStream('mt19937ar','Seed',sum(100*clock));
RandStream.setDefaultStream(stream);



% For neighbor calculation
%GRID= 1:grid_OUTX*grid_OUTY;
%GRID= reshape(GRID, [grid_OUTX grid_OUTY]);
% pad with -1
%Padded_grid= padarray(GRID,[max_radius_val max_radius_val],-1);


hits= zeros(grid_OUTX*grid_OUTY, 1);




%% Initialization of V4

% jitter the center postions
nlim= 0.25;
a=-nlim; b=nlim;

WT=        zeros(grid_OUTX, grid_OUTY, ip_dimen(1), ip_dimen(2), sensory_types, 'single');
%WT_matrix= zeros(ip_dimen(1)*grid_OUTX, ip_dimen(2)*grid_OUTY, sensory_types, 'single');

% Initialize with Coarse "blobby" map
for OUTX=1:1:grid_OUTX
    for OUTY=1:1:grid_OUTY        
      
        % make a blob
        %sensory_ip= exp(-(((xo_v(OUTX, OUTY)-y).^2 + (yo_v(OUTX, OUTY)-x).^2))./(2*init_stim_siigma*init_stim_siigma));        
        
        % make a sine gabor
        or= 359*rand(1);
        ch = randi(3, 1,1);
        if (ch==1)
        sensory_ip=  make_sin_gabor_adjust_rfc(gabor_size, 1, or -90, strght_rf_AR, 0, 0);  
        elseif (ch==2)
        sensory_ip=  curved_sin_gabor_light_out(gabor_size, 2, or, curved_rf_AR_1, curved_r1, 0, 0); 
        elseif (ch==3)
        sensory_ip=  curved_sin_gabor_dark_out(gabor_size, 2, or, curved_rf_AR_1, curved_r1, 0, 0);
        end;
        
        % add some noise
        random_jitter= a+ (b-a).*rand(size(sensory_ip));
        sensory_ip= sensory_ip+ random_jitter;

        
        WT(OUTX, OUTY, :,:,:) = random_jitter;

       
    end;
end;






%% Natural images database filled

%gabor_size = 7;
patch_size = 33;


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
im_count= 1200000;
count= .03* im_count;


%D= v1nbrhood_actv(r+1, r+1, :); % get the center of the box ...
v1nbrhood_actv=[]; % save MEMORY


% NO SORTING ....
sortids= randperm(im_count);
idx_start= 0;

top_folder=     img_folder(sortids(idx_start+ 1: idx_start+ count));
top_subfolder=  img_subfolder(sortids(idx_start+ 1: idx_start+ count));
top_positions=  img_positions(sortids(idx_start+ 1: idx_start+ count), :);


% random indices
r_ids = randi(count, [MV_max 1]);






%% MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic


MV_counter= 1;
%repeats_per_image= 8;


% reshape WT matrix for learning (in vectorized way)
WT= reshape(WT, grid_OUTX*grid_OUTY, ip_dimen(1)*ip_dimen(2)*sensory_types);

oriens = linspace(0, 360, 16);
L= length(oriens);



while (MV_counter <=MV_max)
    
    %INPUT reading from NI database
    Irotated= image_read(top_folder(r_ids(MV_counter)), top_subfolder(r_ids(MV_counter)), dirinfo, subdirinfo);
    
    % pick random orientations
    % r_orien = rand(repeats_per_image, 1)*359;
    r_orien= randi(L-1, 1, 1);
    
    X = Irotated(top_positions(r_ids(MV_counter), 1)-len: top_positions(r_ids(MV_counter), 1)+len, ...
                 top_positions(r_ids(MV_counter), 2)-len: top_positions(r_ids(MV_counter), 2)+len, ...
                 top_positions(r_ids(MV_counter), 3) );
    
    
    im= imrotate(X, oriens(r_orien), 'crop');
    
    img_xc= floor(patch_size/2);
    img_yc= floor(patch_size/2);
    sensory_ip= im(img_xc- floor(feature_dimX/2): img_xc+ floor(feature_dimX/2), img_yc- floor(feature_dimX/2): img_yc+ floor(feature_dimX/2));
    
    
    
    
    %% Do the CORRELATION match
    
    
    % c. normalized sensory ip
    sensory_ip_0_mean = sensory_ip - mean(sensory_ip(:));
    sqrt_sensory_ip   = sqrt((sensory_ip_0_mean(:))'*sensory_ip_0_mean(:));
    sensory_ip_norm   = sensory_ip_0_mean ./sqrt_sensory_ip;
    
    
    % d. do the correlation calculation
    weight_vectors_0_mean= bsxfun(@minus,  WT,   mean(WT, 2));
    %sqrt_weight_vectors=   sqrt( dot(weight_vectors_0_mean', weight_vectors_0_mean') );
    sqrt_weight_vectors=   sqrt( sum(weight_vectors_0_mean'.* weight_vectors_0_mean') );
    weight_vectors_norm=   bsxfun(@rdivide, weight_vectors_0_mean', sqrt_weight_vectors);
    
    correlation_matrix=  weight_vectors_norm'* sensory_ip_norm(:);
    
    
    
    
    
    %% LEARNING ALGORITHM
    
    
    % COMBINED PHASE 1 and 2
    
    [w_values sorted_w_indices]= sort(correlation_matrix(:), 'descend');
    
    % Compute the ID
    pidx= ceil((MV_counter)/(MV_max) *length(LR));
    limit_ids= ceil(pval(pidx)/100*total_nbrs);
    
    
    % get the coordinantes
    coords= sorted_w_indices(1:limit_ids);
    
    juices= w_values(1:limit_ids);
    
    
    % Compute the ID
    idx= ceil(MV_counter/MV_max *length(LR));
    
    % Adjust the LR of the winner ids
    R_learning= LR(idx).*juices;
    
    % data logging
    % percent_training(MV_counter)= limit_ids/total_nbrs;
    % mean_pop_response(MV_counter)= mean(juices);     
    hits(sorted_w_indices(1:limit_ids))= hits(sorted_w_indices(1:limit_ids)) +1;
    
    
    learning_WTs= WT(coords, :);
    
    
    %% save the results before the weight matrix gets changed
    if (nnz(save_step==MV_counter)>0)
        
        % file save ...
        file_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/development_with_NI_images/data/');
        file_pre= strcat('NI-wts-frame=',num2str(MV_counter,'%06d'), ...
            '-neurons=', num2str(neurons, '%d'), ...
            '-nbd_radius_val=', num2str(r, '%d'), ...
            '-max_radius_val=', num2str(max_radius_val, '%d'), ...
            '-sLR=',num2str(sLR, '%1.3f'), ...
            '-percent=', num2str(percent_val, '%1.2f'), ...
            '-random_inputs');
        
        fOut= strcat(file_path, file_pre,'.mat');
        
        % save the variables
        save(fOut,'WT','correlation_matrix', 'hits');
    end;
    
    
    
    %% Do the learning (vectorized)
    sz_LR_wts= size(learning_WTs);
    sensory_matrix_big= repmat(sensory_ip(:)', [sz_LR_wts(1) 1]);
    
    weight_vectors=   bsxfun(@times, learning_WTs, (1-R_learning) );
    target=           bsxfun(@times, sensory_matrix_big, R_learning );
    WT_v=   weight_vectors+  target;
    
    WT(coords, :) =WT_v;
    
    
    
    MV_counter= MV_counter+ 1;
    
end;


%toc
profile report
