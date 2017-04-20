%% concatenate many images into a large matrix


function [top_image_patches bigImage]= image_displayer(max_radius_val, img_xc, img_yc, top_folder, top_subfolder, top_positions, dirinfo, subdirinfo, len)


boundary= 4;
supplied_parameters_natural_inputs_reduced_v4_size;

% for zooming IN, see the red circle clearly
% visV1left=  11;
% visV1right= 23;


nRows = 10;
nCols = 10;
imgCell = cell(nRows,nCols);

for t = 1:nRows*nCols

    Irotated= image_read(top_folder(t), top_subfolder(t), dirinfo, subdirinfo);
    
    top_image_patches(t, :,:) = Irotated(top_positions(t, 1)-len: top_positions(t, 1)+len, ...
                                         top_positions(t, 2)-len: top_positions(t, 2)+len, ...
                                         top_positions(t, 3) );
    
    im= squeeze(top_image_patches(t,:,:));
    
    im_mean= mean(im(:));
    im= im- im_mean;    
    im_block= im(img_xc-boundary: img_xc+boundary, img_yc-boundary: img_yc+boundary);
    min_box= min(im_block(:));
    max_box= max(im_block(:));
    im_new= (im- min_box)./(max_box- min_box)*(255 -0)+ 0;
    
    im_new= im_new(img_xc- (max_radius_val+1 +3 ): img_xc+ (max_radius_val+1 +3 ), ...
                   img_yc- (max_radius_val+1 +3 ): img_yc+ (max_radius_val+1 +3 ));    
    
    
    %im_block= im(12-3: 12+3, 12-3: 12+3);
    %im_std= std(im_block(:));    
    %im_std= std(im(:));   
    %im_stdzied= 30*(im./im_std);    
    %im_final= im_stdzied+ 200;
    
    %# add the image to imgCell
    %# images will filled first into all rows of column one
    %# then into all rows of column 2, etc
    imgCell{t} = im_new;

end

%# if you want the images to be arranged along rows instead of 
%# columns, you can transpose imgCell here
%# imgCell = imgCell';

%# catenate into big image
bigImage = cell2mat(imgCell);

