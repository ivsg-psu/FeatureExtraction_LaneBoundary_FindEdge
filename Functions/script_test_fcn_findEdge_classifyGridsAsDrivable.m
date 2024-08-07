% script_test_fcn_findEdge_classifyGridsAsDrivable.m
% Exercises the function: fcn_findEdge_classifyGridsAsDrivable.m
% REVISION HISTORY: 
% This function was written on 2024_07_15 by Aneesh Batchu
% -- Seperated this code from fcn_findEdge_classifyGridsAsDrivable
% Funclionalize this code on 7/16/2024 by Jiabao Zhao
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it

% figure number
fig_num_LLA = 3001;

% LiDAR data
LiDAR_outer_edge = vertcat(LiDAR_ENU_cell{1400:1410});

% Input points (LiDAR data)
input_points = LiDAR_outer_edge(:,1:2); 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1.25; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_findEdge_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries
% grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

grid_boundaries = [Min_x Max_x Min_y Max_y]; 

% As of now. 
point_density = 20;

input_points = LiDAR_outer_edge; 

std_threshold = []; 

theta_threshold = 7*pi/180;


% These two function are used here because some input of this test depends on
% the output of them
[gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, grids_with_zero_points,...
grids_greater_than_zero_points, gridCenters_zero_point_density,...
gridCenters_greater_than_zero_point_density] = fcn_findEdge_findGridsWithPoints(input_points,...
grid_size,grid_boundaries,-1);

[original_grids_with_required_point_density, gridCenters_low_point_density,...
gridCenters_required_point_density, current_grids_with_low_point_density,...
current_grids_with_required_point_density] = fcn_findEdge_GridsIntoMappedUnmapped...
(point_density, total_N_points_in_each_grid, grids_greater_than_zero_points, gridCenters,-1);


% Below is the MAIN function 
[standard_deviation_in_z, angle_btw_unit_normals_and_vertical, original_drivable_grids,...
original_non_drivable_grids, current_drivable_grid_numbers_in_mapped_grids,...
current_non_drivable_grid_numbers_in_mapped_grids, gridCenters_drivable_grids,...
gridCenters_non_drivable_grids] = fcn_findEdge_classifyGridsAsDrivable...
(original_grids_with_required_point_density,input_points,std_threshold,theta_threshold,gridCenters,fig_num_LLA);
