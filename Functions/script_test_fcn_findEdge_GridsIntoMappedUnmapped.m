% script_test_fcn_findEdge_GridsIntoMappedUnmapped.m
% Exercises the function: fcn_findEdge_GridsIntoMappedUnmapped.m
% Revision history:
% 2024_07_15 - Aneesh Batchu
% -- Wrote the example 
% 2024_07_15 - Jiabao Zhao
% -- Functionlize the code 
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it
%% INPUTS

% figure number
fig_num_LLA = 3001;

% LiDAR data
LiDAR_outer_edge = vertcat(LiDAR_ENU_cell{1400:1410});

% Input points (LiDAR data)
input_points = LiDAR_outer_edge; 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1.25; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries
grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

% As of now. 
point_density = 60;

% This function is used here because some input of this test depends on
% the output of it
[gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, grids_with_zero_points,...
grids_greater_than_zero_points, gridCenters_zero_point_density,...
gridCenters_greater_than_zero_point_density] = fcn_geometry_findGridsWithPoints(input_points,...
grid_size,grid_boundaries,-1);


% Below this the MAIN function 
[original_grids_with_required_point_density, gridCenters_low_point_density,...
gridCenters_required_point_density, current_grids_with_low_point_density,...
current_grids_with_required_point_density] = fcn_geometry_GridsIntoMappedUnmapped...
(point_density, total_N_points_in_each_grid, grids_greater_than_zero_points, gridCenters,fig_num_LLA);

assert(iscell(gridIndices_cell_array))
assert(isequal(length(gridCenters(1,:)),2))
assert(isequal(length(gridIndices(1,:)),1))
assert(isequal(length(grid_AABBs(1,:)),4))
