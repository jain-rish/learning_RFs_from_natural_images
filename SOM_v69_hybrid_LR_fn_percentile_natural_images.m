% Rishabh Jain
% LNC
% April 2015



function SOM_v69_hybrid_LR_fn_percentile_natural_images(neurons, max_radius_val, sLR, spikenoise, percent_val)
% clc;
% clear all;
% profile on;

grid_OUTX= round(sqrt(neurons));
grid_OUTY= round(sqrt(neurons));

sensory_types= 1;
feature_dimX=  9;
feature_dimY=  9;
ip_dimen= [feature_dimX feature_dimY];



%% General parameters and file reading

filters= 'mixed';
map= 'hybrid';
r= 1;

% file read ...
knee= 40;
data_read_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/spring_2015_data/hybrid-data-awesome/single-scale/');
s_pre= strcat('v1-neighborhood-activity', ...
    '-max_radius_val_', num2str(r, '%d'), ...
    '-knee_', num2str(knee, '%d'),'-substrate_', filters, '_norm_', map, '_NNlinear');


fOut= strcat(data_read_path, s_pre,'.mat');
load(fOut);







%% Update parameters with supplied ones from UNIX script

% simulation running time ...

% fr= 75000; 
% MV_max= fr*sensory_types; 
% 
% % save more early on ...
% save_step1= round(linspace(1, MV_max/5, 20));
% save_step2= round(linspace(MV_max/5+save_step1(end)-save_step1(end-1), MV_max, 5));
% save_step= [save_step1 save_step2];

MV_max= 300000; 
save_step= round(linspace(1, MV_max, 3));


sampling= round(MV_max/10);

%t_phase1= round(MV_max/100);
%pval(1:t_phase1)= 100; % phase 1
pval(1:sampling)= linspace(percent_val, 0.025, sampling); % Assuming 0.025% is 1 neuron


%LR(1:t_phase1)= sLR; % phase 1
LR_space= linspace(1, 0.2, 100); % 5-fold jump in 1 to 100%
LR_curve= ( linspace(LR_space(round(percent_val*150)), 1, sampling) ); % lower the percentile, higher the jump
LR(1:sampling)= sLR .*LR_curve; % scaling LR


total_nbrs= round((2*max_radius_val+1)*(2*max_radius_val+1) /sensory_types);

% radius= round(linspace(max_radius_val, end_radius_val, sampling));
warning off all;







%% Pre-allocation of arrays

percent_training=  Inf([MV_max 1]);
mean_pop_response= Inf([MV_max 1]);

% set the random streams based on your parameter
stream = RandStream('mt19937ar','Seed',sum(100*clock));
RandStream.setDefaultStream(stream);



% For neighbor calculation
GRID= 1:grid_OUTX*grid_OUTY;
GRID= reshape(GRID, [grid_OUTX grid_OUTY]);
% pad with -1
Padded_grid= padarray(GRID,[max_radius_val max_radius_val],-1);


hits= zeros(grid_OUTX, grid_OUTY);
hits_winners= zeros(grid_OUTX, grid_OUTY);

% tracking related
track_radius= 1;
track_ids_Xlist= 25-track_radius: 25+track_radius;
track_ids_Ylist= 25-track_radius: 25+track_radius;

for xid=1:length(track_ids_Xlist)
    for yid=1:length(track_ids_Ylist)
        tracking_ids(xid, yid).time=1;
    end
end





%% Initialization of V4

vleft= round(feature_dimX/2); vright= round(feature_dimX/2);
[xo_v yo_v]= meshgrid( linspace(vleft, vright, grid_OUTX), linspace(vleft, vright, grid_OUTY));

% jitter the center postions
nlim= 0.1;
a=-nlim; b=nlim;

random_jitter_x= a+ (b-a).*rand(size(xo_v));
xo_v= xo_v+ random_jitter_x;

random_jitter_y= a+ (b-a).*rand(size(yo_v));
yo_v= yo_v+ random_jitter_y;

[x y]= meshgrid( linspace(1, feature_dimX, ip_dimen(1)), linspace(1, feature_dimY, ip_dimen(2)));

WT=      zeros(grid_OUTX, grid_OUTY, ip_dimen(1), ip_dimen(2), sensory_types, 'single');
WT_matrix= zeros(ip_dimen(1)*grid_OUTX, ip_dimen(2)*grid_OUTY, sensory_types, 'single');

% Initialize with Coarse "blobby" map
for OUTX=1:1:grid_OUTX
    for OUTY=1:1:grid_OUTY        
      
        % make a blob
        %sensory_ip= exp(-(((xo_v(OUTX, OUTY)-y).^2 + (yo_v(OUTX, OUTY)-x).^2))./(2*init_stim_siigma*init_stim_siigma));        
        
        % make a sine gabor
        or= 359*rand(1);
        sensory_ip= make_sin_gabor_adjust_rfc(9, 1, or, 2, 0, 0); 
        WT(OUTX, OUTY, :,:,:) = sensory_ip;

       
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
count= 750;


%D= v1nbrhood_actv(xs, ys, :);
D= v1nbrhood_actv(r+1, r+1, :); % get the center of the box ...

[sortvals, sortids] = sort(squeeze(D), 'descend');

idx_start= 0;

top_folder=     img_folder(sortids(idx_start+ 1: idx_start+ count));
top_subfolder=  img_subfolder(sortids(idx_start+ 1: idx_start+ count));
top_positions=  img_positions(sortids(idx_start+ 1: idx_start+ count), :);

% random indices
r_ids = randi(count, [MV_max 1]);





%% MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic



for MV_counter=1:MV_max    
    
    % reshape WT matrix for easy access to (likely) winner neighborhood
    WT= reshape(WT, [grid_OUTX, grid_OUTY, ip_dimen(1), ip_dimen(2), sensory_types]);
    
    
    %% INPUT reading from NI database
    
        
    Irotated= image_read(top_folder(r_ids(MV_counter)), top_subfolder(r_ids(MV_counter)), dirinfo, subdirinfo);
    
    X = Irotated(top_positions(r_ids(MV_counter), 1)-len: top_positions(r_ids(MV_counter), 1)+len, ...
                                                          top_positions(r_ids(MV_counter), 2)-len: top_positions(r_ids(MV_counter), 2)+len, ...
                                                          top_positions(r_ids(MV_counter), 3) );
    
    % pick a random orientation
    r_orien = rand(1)*359;                                               
    X= imrotate(X, r_orien, 'crop');
    im= squeeze(X);

    img_xc= floor(patch_size/2);
    img_yc= floor(patch_size/2);
    sensory_ip= im(img_xc- floor(feature_dimX/2): img_xc+ floor(feature_dimX/2), img_yc- floor(feature_dimX/2): img_yc+ floor(feature_dimX/2));
    
    
    
    
    %% Do the CORRELATION match                
    
    % a. Based on input, calculate approx. neural coordinates of the winner
    %MxI= vleft;  MyI= vleft;
    %MxA= vright; MyA= vright;
       
    % normalize the visual coordinates
    %likely_winy = round((xy(1) - MxI)/(MxA - MxI)*(grid_OUTX-1)+1);
    %likely_winx = round((xy(2) - MyI)/(MyA - MyI)*(grid_OUTY-1)+1);
    
    likely_winy = round(grid_OUTX/2);
    likely_winx = round(grid_OUTY/2);
    
    % b. calculate neural coordinates of the local region
    localXs= likely_winx-max_radius_val:likely_winx+ max_radius_val;
    localYs= likely_winy-max_radius_val:likely_winy+ max_radius_val;
    
    localXs= localXs(localXs>0 & localXs <=grid_OUTX);
    localYs= localYs(localYs>0 & localYs <=grid_OUTY);
    
    %local_correlation_matrix= -1*ones(length(localXs), length(localYs));    
    
              
  
    % c. normalized sensory ip
    sensory_ip_0_mean = sensory_ip - mean(sensory_ip(:));
    sqrt_sensory_ip   = sqrt((sensory_ip_0_mean(:))'*sensory_ip_0_mean(:));
    sensory_ip_norm   = sensory_ip_0_mean ./sqrt_sensory_ip;
    
    
    % d. do the correlation calculation
    local_WTs= WT(localXs, localYs,:,:,:);
    sz_local_WTs= size(local_WTs);
    reshaped_local_WTs= reshape(local_WTs, [sz_local_WTs(1)*sz_local_WTs(2) sz_local_WTs(3)*sz_local_WTs(4)*sensory_types]);    
    weight_vectors_0_mean= bsxfun(@minus,  reshaped_local_WTs,   mean(reshaped_local_WTs, 2));
    %sqrt_weight_vectors=   sqrt( dot(weight_vectors_0_mean', weight_vectors_0_mean') );    
    sqrt_weight_vectors=   sqrt( sum(weight_vectors_0_mean'.* weight_vectors_0_mean') );        
    weight_vectors_norm=   bsxfun(@rdivide, weight_vectors_0_mean', sqrt_weight_vectors);
    
    local_correlation_matrix=  weight_vectors_norm'* sensory_ip_norm(:);
    local_correlation_matrix=  reshape(local_correlation_matrix, [length(localXs) length(localYs) ]);
    
    
    
    
    
    %% LEARNING ALGORITHM
    
    % April 2013- correlation matrix needs zeros not -1 
    correlation_matrix= zeros(grid_OUTX, grid_OUTY);
    % copy the local calculation to the BIG "fake" correlation matrix
    correlation_matrix(localXs, localYs)= local_correlation_matrix;
    
    % choose a winner in the map based on "correlation score"
    [spatial_winner spatial_w_index]= max(correlation_matrix(:));
    [SP_WINX SP_WINY]= ind2sub(size(correlation_matrix), spatial_w_index);
    
    % document the winner    
    hits_winners(SP_WINX, SP_WINY)= hits_winners(SP_WINX, SP_WINY) +1;
    
    % Add the offset
    WINX_off= SP_WINX+ max_radius_val;
    WINY_off= SP_WINY+ max_radius_val;
    
      
    % List of spatial neighbors
    neighbors= Padded_grid((WINX_off)-max_radius_val:(WINX_off)+max_radius_val, (WINY_off)-max_radius_val:(WINY_off)+max_radius_val);
    neighbors= neighbors(neighbors~=-1);
 
    nbd_correlation_matrix= correlation_matrix(neighbors);  
    


    % COMBINED PHASE 1 and 2
    
    %add NOISE to correlation values
    random_noise= normrnd(0, spikenoise, size(nbd_correlation_matrix)); % GAUSSIAN NOISE
    nbd_correlation_matrix= nbd_correlation_matrix + random_noise;
    
    
    [w_values sorted_w_indices]= sort(nbd_correlation_matrix(:), 'descend');
    
    % Compute the ID
    pidx= ceil((MV_counter)/(MV_max) *length(LR));
    limit_ids= round(pval(pidx)/100*total_nbrs);
    
    % to avoid index exceeding dimensions "error"
    if (limit_ids> length(w_values))
        limit_ids= length(w_values);
    end
    
    % get the coordinantes
    coords= neighbors(sorted_w_indices(1:limit_ids));
    [L1X L1Y]= ind2sub(size(correlation_matrix), coords);
    
    juices= w_values(1:limit_ids);
    
    
    % Compute the ID
    idx= ceil(MV_counter/MV_max *length(LR));
    
    % Adjust the LR of the winner ids
    R_learning= LR(idx).*juices;
    
    % data logging
    percent_training(MV_counter)= limit_ids/total_nbrs;
    mean_pop_response(MV_counter)= mean(juices);
    
     
    
    % reshape WT matrix for learning (in vectorized way)
    WT= reshape(WT, grid_OUTX*grid_OUTY, ip_dimen(1)*ip_dimen(2)*sensory_types);
    learning_WTs= WT(sub2ind(size(correlation_matrix), L1X, L1Y),:); 
        
    

    
    %% Visualization
%     % save the results before the weight matrix gets changed
if (nnz(save_step==MV_counter)>0)
    %
    %         close all;
    %
    %         % reshape WT for visualization
    %         for OUTX=1:1:grid_OUTX
    %             for OUTY=1:1:grid_OUTY
    %
    %                 [cell_n]= sub2ind([size(GRID)], OUTX, OUTY);
    %
    %                 rf= WT(cell_n, :);
    %                 rf= reshape(rf, [ip_dimen(1), ip_dimen(2), sensory_types]);
    %
    %                 wtX_strt=  (OUTX-1)*ip_dimen(1)+1;
    %                 wtX_end=   (OUTX-1)*ip_dimen(1)+ip_dimen(1);
    %                 wtY_strt=  (OUTY-1)*ip_dimen(2)+1;
    %                 wtY_end=   (OUTY-1)*ip_dimen(2)+ip_dimen(2);
    %
    %                 WT_matrix(wtX_strt: wtX_end, wtY_strt:wtY_end, :)= rf;
    %
    %             end
    %         end
    %
    %
        

    % file save ...
    file_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/development_with_NI_images/');
    file_pre= strcat('NI-wts-frame=',num2str(MV_counter,'%06d'), ...;
        '-neurons=', num2str(neurons, '%d'), ...
        '-max_radius_val=', num2str(max_radius_val, '%d'), ...
        '-sLR=',num2str(sLR, '%1.3f'), ...
        '-percent=', num2str(percent_val, '%1.2f'));
    fOut= strcat(file_path, file_pre,'.mat');
    
    % save the variables
    save(fOut,'WT','correlation_matrix', 'hits', 'hits_winners', 'percent_training', 'mean_pop_response', 'tracking_ids');           
    
    
    %         center_of_mass_fit_ellipse_v10('hybrid',s_path, WT_matrix, MV_counter, hits, hits_winners, percent_training, mean_pop_response, tracking_ids);
    %
    %         % reshape WT again for simulations
    %         for OUTX=1:1:grid_OUTX
    %             for OUTY=1:1:grid_OUTY
    %
    %                 wtX_strt=  (OUTX-1)*ip_dimen(1)+1;
    %                 wtX_end=   (OUTX-1)*ip_dimen(1)+ip_dimen(1);
    %                 wtY_strt=  (OUTY-1)*ip_dimen(2)+1;
    %                 wtY_end=   (OUTY-1)*ip_dimen(2)+ip_dimen(2);
    %
    %                 rf= WT_matrix(wtX_strt: wtX_end, wtY_strt:wtY_end, :);
    %
    %                 [cell_n]= sub2ind([size(GRID)], OUTX, OUTY);
    %
    %                 WT(cell_n, :)= rf(:);
    %
    %             end
    %         end
    %
    %     end;
    %
    %
    %     % tracking a neighborhood
    %     for xid= 1:length(track_ids_Xlist)
    %         for yid= 1:length(track_ids_Ylist)
    %
    %             if (nnz(SP_WINX== track_ids_Xlist(xid))>0 && nnz(SP_WINY== track_ids_Ylist(yid))>0 )
    %
    %                 T= tracking_ids(xid, yid).time;
    %                 %tracking_ids(xid, yid).snap(T,:)= snapshotter(correlation_matrix, neighbors, neighbors(winner_ids), local_threshold_counter(xid, yid));
    %                 tracking_ids(xid, yid).corr_val(T)= correlation_matrix(SP_WINX, SP_WINY);
    %                 tracking_ids(xid, yid).MV_counter(T)= MV_counter;
    %
    %                 if(isempty(max(R_learning(:))))
    %                     tracking_ids(xid, yid).LR(T)= 0;
    %                 else
    %                     tracking_ids(xid, yid).LR(T)= max(R_learning(:));
    %                 end
    %
    %                 T= T+1;
    %                 tracking_ids(xid, yid).time= T;
    %
    %             end;
    %         end;
    %     end
end;
    
    
            
    
    
    %% Do the learning (vectorized) 
   sz_LR_wts= size(learning_WTs);
   sensory_matrix_big= repmat(sensory_ip(:)', [sz_LR_wts(1) 1]);
    
   weight_vectors=   bsxfun(@times, learning_WTs, (1-R_learning) );
   target=           bsxfun(@times, sensory_matrix_big, R_learning );
   WT_v=   weight_vectors+  target;
   
   WT(sub2ind(size(correlation_matrix), L1X, L1Y),:) =WT_v;

            

end;


%toc
profile report
