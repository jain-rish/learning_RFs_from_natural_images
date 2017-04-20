function WT_matrix= WTs_visualize_linear(grid_OUTX, grid_OUTY, WTs_sorted, feature_dimX)

%grid_OUTX= round(sqrt(neurons));
%grid_OUTY= round(sqrt(neurons));

sensory_types= 1;
%feature_dimX=  9;
%feature_dimY=  9;
ip_dimen= [feature_dimX feature_dimX];

GRID= 1:grid_OUTX*grid_OUTY;
GRID= reshape(GRID, [grid_OUTX grid_OUTY]);


% reshape WT for visualization
for OUTX=1:1:grid_OUTX
    for OUTY=1:1:grid_OUTY
    
    [cell_n]= sub2ind([size(GRID)], OUTX, OUTY);
    
    rf= WTs_sorted(cell_n, :);
    rf= reshape(rf, [ip_dimen(1), ip_dimen(2), sensory_types]);
    
    wtX_strt=  (OUTX-1)*ip_dimen(1)+1;
    wtX_end=   (OUTX-1)*ip_dimen(1)+ip_dimen(1);
    wtY_strt=  (OUTY-1)*ip_dimen(2)+1;
    wtY_end=   (OUTY-1)*ip_dimen(2)+ip_dimen(2);
    
    WT_matrix(wtX_strt: wtX_end, wtY_strt:wtY_end, :)= rf;
    end;
end;

imshow(WT_matrix', []);

N = feature_dimX *grid_OUTX;
M = feature_dimX *grid_OUTY;
x = linspace(0.5, N+0.5, N/feature_dimX+1);
y = linspace(0.5, M+0.5, M/feature_dimX+1);

% Horizontal grid 
for k = 1:length(y)
  line([x(1) x(end)], [y(k) y(k)])
end

% Vertical grid
for k = 1:length(x)
  line([x(k) x(k)], [y(1) y(end)])
end
%axis square
