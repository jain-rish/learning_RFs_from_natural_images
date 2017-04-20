% Rishabh Jain
% Sorting RFs

neurons= 1225;
MV_counter= 300000;
max_radius_val= 35;
sLR= 0.1;
percent_val= 0.12;

grid_OUTX= round(sqrt(neurons));
grid_OUTY= round(sqrt(neurons));

sensory_types= 1;
feature_dimX=  9;
feature_dimY=  9;
ip_dimen= [feature_dimX feature_dimY];

GRID= 1:grid_OUTX*grid_OUTY;
GRID= reshape(GRID, [grid_OUTX grid_OUTY]);

%file_path= strcat('/amnt/foam/foamd0/rishabh/Phase_2_multimap/NI_search_many_filters_radius_expts/development_with_NI_images/');
file_pre= strcat('NI-wts-frame=',num2str(MV_counter,'%06d'), ...;
    '-neurons=', num2str(neurons, '%d'), ...
    '-max_radius_val=', num2str(max_radius_val, '%d'), ...
    '-sLR=',num2str(sLR, '%1.3f'), ...
    '-percent=', num2str(percent_val, '%1.2f'));
fOut= strcat(file_pre,'.mat');



load(fOut);
% flatten the WT vector
sz_local_WTs= size(WT);

% 
% WTs_sorted= WT;
% for OUTX=1:1:grid_OUTX
%     for OUTY=1:1:grid_OUTY
%         
%         [cell_n]= sub2ind([size(GRID)], OUTX, OUTY);
%         img_candidate=  WTs_sorted(cell_n, :);
%          
%         weight_vectors_0_mean= bsxfun(@minus,  WTs_sorted,   mean(WT, 2));
%         sqrt_weight_vectors=   sqrt( sum(weight_vectors_0_mean'.* weight_vectors_0_mean') );
%         weight_vectors_norm=   bsxfun(@rdivide, weight_vectors_0_mean', sqrt_weight_vectors);
%         
%         img_0_mean = img_candidate - mean(img_candidate(:));
%         sqrt_img   = sqrt((img_0_mean(:))'*img_0_mean(:));
%         img_norm   = img_0_mean ./sqrt_img;
%                 
%         img_correlation_matrix=  weight_vectors_norm'* img_norm(:);
%         [w_values sorted_w_indices]= sort(img_correlation_matrix(:), 'descend');
%         WTs_sorted= WTs_sorted(sorted_w_indices, :);
%         
%     end;
% end;

%%

boundary=4;
img_xc= 5;
img_yc= 5;


% sort by intensity first
weight_vectors_min= min(WT, [], 2);
weight_vectors_max= max(WT, [], 2);

wts_range= weight_vectors_max- weight_vectors_min;

WT_min_subtracted= bsxfun(@minus,  WT, weight_vectors_min);
WTs_sorted= bsxfun(@rdivide,  WT_min_subtracted, wts_range);
WTs_sorted= WTs_sorted.*255;


clusters = kmeans(WTs_sorted, 250, 'distance', 'correlation');
[~, clusters_ind] = sort(clusters, 'ascend');


for i= 1:1:grid_OUTX*grid_OUTY
    clusters(i) = entropy(double(WTs_sorted(i,:)));
end;

[~, clusters_ind] = sort(clusters, 'ascend');

% [X Y]= meshgrid(1:grid_OUTX, 1:grid_OUTY);
% d= sqrt(X.^2+ Y.^2);
% [~,d_ind] = sort(d(:));


% reshape WT for visualization
for OUTX=1:1:grid_OUTX
    for OUTY=1:1:grid_OUTY
        
        [cell_n]= sub2ind([size(GRID)], OUTX, OUTY);
        
        rf= WTs_sorted(clusters_ind(cell_n), :);
        rf= reshape(rf, [ip_dimen(1), ip_dimen(2), sensory_types]);
        
        wtX_strt=  (OUTX-1)*ip_dimen(1)+1;
        wtX_end=   (OUTX-1)*ip_dimen(1)+ip_dimen(1);
        wtY_strt=  (OUTY-1)*ip_dimen(2)+1;
        wtY_end=   (OUTY-1)*ip_dimen(2)+ip_dimen(2);
        
        %im= rf;
        
        %im_mean= mean(im(:));
        %im= im- im_mean;
        %im_block= im(img_xc-boundary: img_xc+boundary, img_yc-boundary: img_yc+boundary);
        %min_box= min(im_block(:));
        %max_box= max(im_block(:));
        %im_new= (im- min_box)./(max_box- min_box)*(255 -0)+ 0;
        
        
        WT_matrix(wtX_strt: wtX_end, wtY_strt:wtY_end, :)= rf;
        
        
    end
end

figure; imshow(WT_matrix', []);




N = feature_dimX *sqrt(neurons);
x = linspace(0.5, N+0.5, N/feature_dimX+1);
y = linspace(0.5, N+0.5, N/feature_dimY+1);

% Horizontal grid 
for k = 1:length(y)
  line([x(1) x(end)], [y(k) y(k)])
end

% Vertical grid
for k = 1:length(y)
  line([x(k) x(k)], [y(1) y(end)])
end

axis square



