function WTs_sorted= pairwise_image_sorting(neurons, seed, WT, WTs_normed, feature_dimX)

grid_OUTX= round(sqrt(neurons));
grid_OUTY= round(sqrt(neurons));

sensory_types= 1;
%feature_dimX=  9;
%feature_dimY=  9;
ip_dimen= [feature_dimX feature_dimX];

GRID= 1:grid_OUTX*grid_OUTY;
GRID= reshape(GRID, [grid_OUTX grid_OUTY]);



% pairwise-sorting
WTs_sorted= nan(size(WT));
img_candidate=  WT(seed, :);
temp= WTs_normed(1:seed, :);
WTs_normed(1:seed, :)= nan;

for cell_n= seed:1:grid_OUTX*grid_OUTY-1   
    weight_vectors_0_mean= bsxfun(@minus,  WTs_normed,   mean(WTs_normed, 2));
    sqrt_weight_vectors=   sqrt( sum(weight_vectors_0_mean'.* weight_vectors_0_mean') );
    weight_vectors_norm=   bsxfun(@rdivide, weight_vectors_0_mean', sqrt_weight_vectors);
    
    img_0_mean = img_candidate - mean(img_candidate(:));
    sqrt_img   = sqrt((img_0_mean(:))'*img_0_mean(:));
    img_norm   = img_0_mean ./sqrt_img;
    
    img_correlation_matrix=  (weight_vectors_norm'* img_norm(:));
    img_correlation_matrix(isnan(img_correlation_matrix)) = 0;
    [w_values sorted_w_indices]= sort(img_correlation_matrix(:), 'descend');
    
    WTs_sorted(cell_n+1, :)= WTs_normed(sorted_w_indices(1), :);
    WTs_normed(sorted_w_indices(1), :)= nan;
    img_candidate=  WT(sorted_w_indices(1), :);
end;

WTs_sorted(1:seed, :)= temp;
% [X Y]= meshgrid(1:grid_OUTX, 1:grid_OUTY);
% d= sqrt(X.^2+ Y.^2);
% [~,d_ind] = sort(d(:));



