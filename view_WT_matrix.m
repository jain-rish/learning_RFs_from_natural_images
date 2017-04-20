%%

neurons= 625;

grid_OUTX= round(sqrt(neurons));
grid_OUTY= round(sqrt(neurons));

% For neighbor calculation
GRID= 1:grid_OUTX*grid_OUTY;
GRID= reshape(GRID, [grid_OUTX grid_OUTY]);


sensory_types= 1;
feature_dimX=  9;
feature_dimY=  9;
ip_dimen= [feature_dimX feature_dimY];

boundary=4;
img_xc= 5;
img_yc= 5;

% reshape WT for visualization
for OUTX=1:1:grid_OUTX
    for OUTY=1:1:grid_OUTY
        
        [cell_n]= sub2ind([size(GRID)], OUTX, OUTY);
        
        rf= WT(cell_n, :);
        rf= reshape(rf, [ip_dimen(1), ip_dimen(2), sensory_types]);
        
        wtX_strt=  (OUTX-1)*ip_dimen(1)+1;
        wtX_end=   (OUTX-1)*ip_dimen(1)+ip_dimen(1);
        wtY_strt=  (OUTY-1)*ip_dimen(2)+1;
        wtY_end=   (OUTY-1)*ip_dimen(2)+ip_dimen(2);
        
        
        im= rf;
        
        im_mean= mean(im(:));
        im= im- im_mean;
        im_block= im(img_xc-boundary: img_xc+boundary, img_yc-boundary: img_yc+boundary);
        min_box= min(im_block(:));
        max_box= max(im_block(:));
        im_new= (im- min_box)./(max_box- min_box)*(255 -0)+ 0;
        
        %     im_new= im_new(img_xc- (max_radius_val+1 +3 ): img_xc+ (max_radius_val+1 +3 ), ...
        %                    img_yc- (max_radius_val+1 +3 ): img_yc+ (max_radius_val+1 +3 ));
        %
        
        
        WT_matrix(wtX_strt: wtX_end, wtY_strt:wtY_end, :)= im_new;
        
    end
end

figure; imshow(WT_matrix, [])