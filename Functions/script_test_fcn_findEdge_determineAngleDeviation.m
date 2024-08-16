%script_test_fcn_findEdge_determineAngleDeviation
% Exercises the function: fcn_findEdge_determineAngleDeviation
%
% Revision History:
% 2024_08_07 -Jiabao Zhao
% -- wrote the script orginally 
%% Test 1 
% load the data and give figure number
fig_num = 1;
fig_num2 = 11;
fig_num3 = 111;
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), []);


% Extract scan lines
station_1 = 1400;
station_2 = 1450;
range_of_LiDAR = 10;
[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);
scanLineRange = [station1_minus_range_index station2_plus_range_index];
ringsRange = [];

% selected LIDAR range
[~, ~, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]),([]));
[~, ~, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2,LIDAR_intensity,[]);


% prepare for grids
LiDAR_allPoints = [LIDAR_ENU(in_domain,:), LIDAR_scanLineAndRingID(in_domain,:)];
LiDAR_allPoints = LiDAR_allPoints(~isnan(LiDAR_allPoints(:,1)),:);
LIDAR_scanLines = LiDAR_allPoints(:,4:5);
input_points = LiDAR_allPoints(:,1:2);
grid_size = 1; %0.8;%1;%1.25; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,~,~] = fcn_findEdge_findMaxMinOfXYZ(input_points,-1);

% The 2D grid boundaries required for the current analysis
grid_boundaries = [Min_x Max_x Min_y Max_y];

% Seperate the data into grids
[gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, ~,...
    grids_greater_than_zero_points, ~,...
    gridCenters_greater_than_zero_point_density, gridIndices, grid_AABBs] = fcn_findEdge_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,-1);

% Grids with required point density and low point density
point_density = floor(20*((grid_size^2)/(0.3^2)));
[grid_indices_with_required_point_density, ~] = ...
    fcn_findEdge_classifyGridsBasedOnDensity(grids_greater_than_zero_points,total_N_points_in_each_grid,point_density,gridCenters,[],[],[]);

[total_scan_lines_in_each_grid_with_more_than_zero_points] = ....
    fcn_findEdge_calcNumberOfGridScanLines(gridIndices,LIDAR_scanLines,grids_greater_than_zero_points);
[grid_indices_with_more_than_one_scan_line, ~] = ...
    fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters,[],[],[]);

% Finding the orthogonal distances of points in the remaining rings by projecting orthogonally from one ring 
[transverse_span_threshold,transverse_span_each_grid] = fcn_findEdge_determineTransverseSpanThreshold...
    (grids_greater_than_zero_points, grid_AABBs, grid_size, gridIndices, input_points, LIDAR_scanLines,...
    gridCenters_greater_than_zero_point_density,[],[],[]);

% Grids with required point density and low point density
[grid_indices_with_more_than_transverse_span_threshold, ~] =...
    fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid,transverse_span_threshold,grids_greater_than_zero_points,gridCenters,[],[],[]);

% classify Qualified Grids
[current_qualified_grids,~,original_qualified_grids,gridCenters_qualified_grids,~] = ...
    fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,...
    grid_indices_with_more_than_transverse_span_threshold,...
    grids_greater_than_zero_points,gridCenters,[],[],[]);

% "inpolygon" is used to find the grids within the boundary points
shift = 5;
Nscans = length(VehiclePose(:,1));
boundary_points_driven_path = fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose, scanLineRange, Nscans, shift, ([]), ([]));
[in_qg,~] = inpolygon(gridCenters_qualified_grids(:,1),gridCenters_qualified_grids(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));

% Original grid numbers of driven path
% original_grid_numbers_of_driven_path = original_qualified_grids(in_qg);

% Current grid numbers in driven path
current_grid_numbers_of_driven_path = current_qualified_grids(in_qg); %find(in);
gridCenters_driven_path = [gridCenters_qualified_grids(in_qg,1),gridCenters_qualified_grids(in_qg,2)];

% Qualified grid conditions - angle deviation
[angle_btw_unit_normals_and_vertical, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path,...
    angle_btw_unit_normals_and_vertical_driven_path] = fcn_findEdge_determineAngleDeviation...
    (LiDAR_allPoints, gridIndices_cell_array, original_qualified_grids,...
    gridCenters_qualified_grids,current_qualified_grids,gridCenters_driven_path, ...
    grid_AABBs, grid_size, gridIndices, current_grid_numbers_of_driven_path, fig_num,fig_num2,fig_num3);

assert(~isempty(angle_btw_unit_normals_and_vertical));
assert(~isempty(mean_angle_btw_unit_normals_and_vertical_driven_path));
assert(~isempty(angle_btw_unit_normals_and_vertical_driven_path));

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
temp_h2 = figure(fig_num2);
assert(~isempty(get(temp_h2,'Children')))
temp_h3 = figure(fig_num3);
assert(~isempty(get(temp_h3,'Children')))

%% Test 2 no plotting

fig_num = 2;
fig_num2 = 22;
fig_num3 = 222;
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), []);


% Extract scan lines
station_1 = 1400;
station_2 = 1450;
range_of_LiDAR = 10;
[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);
scanLineRange = [station1_minus_range_index station2_plus_range_index];
ringsRange = [];

% selected LIDAR range
[~, ~, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]),([]));
[~, ~, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2,LIDAR_intensity,[]);


% prepare for grids
LiDAR_allPoints = [LIDAR_ENU(in_domain,:), LIDAR_scanLineAndRingID(in_domain,:)];
LiDAR_allPoints = LiDAR_allPoints(~isnan(LiDAR_allPoints(:,1)),:);
LIDAR_scanLines = LiDAR_allPoints(:,4:5);
input_points = LiDAR_allPoints(:,1:2);
grid_size = 1; %0.8;%1;%1.25; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,~,~] = fcn_findEdge_findMaxMinOfXYZ(input_points,-1);

% The 2D grid boundaries required for the current analysis
grid_boundaries = [Min_x Max_x Min_y Max_y];

% Seperate the data into grids
[gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, grids_with_zero_points,...
    grids_greater_than_zero_points, gridCenters_zero_point_density,...
    gridCenters_greater_than_zero_point_density, gridIndices, grid_AABBs] = fcn_findEdge_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,-1);

% Grids with required point density and low point density
point_density = floor(20*((grid_size^2)/(0.3^2)));
[grid_indices_with_required_point_density, gridCenters_low_point_density] = ...
    fcn_findEdge_classifyGridsBasedOnDensity(grids_greater_than_zero_points,total_N_points_in_each_grid,point_density,gridCenters,[],[],[]);

[total_scan_lines_in_each_grid_with_more_than_zero_points] = ....
    fcn_findEdge_calcNumberOfGridScanLines(gridIndices,LIDAR_scanLines,grids_greater_than_zero_points);
[grid_indices_with_more_than_one_scan_line, gridCenters_with_one_scan_line] = ...
    fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters,[],[],[]);

% Finding the orthogonal distances of points in the remaining rings by projecting orthogonally from one ring 
[transverse_span_threshold,transverse_span_each_grid] = fcn_findEdge_determineTransverseSpanThreshold...
    (grids_greater_than_zero_points, grid_AABBs, grid_size, gridIndices, input_points, LIDAR_scanLines,...
    gridCenters_greater_than_zero_point_density,[],[],[]);

% Grids with required point density and low point density
[grid_indices_with_more_than_transverse_span_threshold, gridCenters_with_less_than_transverse_span_threshold] =...
    fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid,transverse_span_threshold,grids_greater_than_zero_points,gridCenters,[],[],[]);

% classify Qualified Grids
[current_qualified_grids,current_unqualified_grids,original_qualified_grids,gridCenters_qualified_grids,gridCenters_unqualified_grids] = ...
    fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,...
    grid_indices_with_more_than_transverse_span_threshold,...
    grids_greater_than_zero_points,gridCenters,[],[],[]);

% "inpolygon" is used to find the grids within the boundary points
shift = 5;
Nscans = length(VehiclePose(:,1));
boundary_points_driven_path = fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose, scanLineRange, Nscans, shift, ([]), ([]));
[in_qg,on_qg] = inpolygon(gridCenters_qualified_grids(:,1),gridCenters_qualified_grids(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));

% Original grid numbers of driven path
original_grid_numbers_of_driven_path = original_qualified_grids(in_qg);

% Current grid numbers in driven path
current_grid_numbers_of_driven_path = current_qualified_grids(in_qg); %find(in);
gridCenters_driven_path = [gridCenters_qualified_grids(in_qg,1),gridCenters_qualified_grids(in_qg,2)];

% Qualified grid conditions - angle deviation
[angle_btw_unit_normals_and_vertical, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path,...
    angle_btw_unit_normals_and_vertical_driven_path] = fcn_findEdge_determineAngleDeviation...
    (LiDAR_allPoints, gridIndices_cell_array, original_qualified_grids,...
    gridCenters_qualified_grids,current_qualified_grids,gridCenters_driven_path, ...
    grid_AABBs, grid_size, gridIndices, current_grid_numbers_of_driven_path, -1,-1,-1);

assert(~isempty(angle_btw_unit_normals_and_vertical));
assert(~isempty(mean_angle_btw_unit_normals_and_vertical_driven_path));
assert(~isempty(angle_btw_unit_normals_and_vertical_driven_path));

% Check that the figure plotted
temp_h = figure(fig_num);
assert(isempty(get(temp_h,'Children')))
temp_h2 = figure(fig_num2);
assert(isempty(get(temp_h2,'Children')))
temp_h3 = figure(fig_num3);
assert(isempty(get(temp_h3,'Children')))
close(fig_num)
close(fig_num2)
close(fig_num3)

%% Test 3 Different grid size

fig_num = 3;
fig_num2 = 33;
fig_num3 = 333;
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), []);


% Extract scan lines
station_1 = 1400;
station_2 = 1450;
range_of_LiDAR = 10;
[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);
scanLineRange = [station1_minus_range_index station2_plus_range_index];
ringsRange = [];

% selected LIDAR range
[~, ~, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]),([]));
[~, ~, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2,LIDAR_intensity,[]);


% prepare for grids
LiDAR_allPoints = [LIDAR_ENU(in_domain,:), LIDAR_scanLineAndRingID(in_domain,:)];
LiDAR_allPoints = LiDAR_allPoints(~isnan(LiDAR_allPoints(:,1)),:);
LIDAR_scanLines = LiDAR_allPoints(:,4:5);
input_points = LiDAR_allPoints(:,1:2);
grid_size = 3; %0.8;%1;%1.25; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,~,~] = fcn_findEdge_findMaxMinOfXYZ(input_points,-1);

% The 2D grid boundaries required for the current analysis
grid_boundaries = [Min_x Max_x Min_y Max_y];

% Seperate the data into grids
[gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, ~,...
    grids_greater_than_zero_points, ~,...
    gridCenters_greater_than_zero_point_density, gridIndices, grid_AABBs] = fcn_findEdge_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,-1);

% Grids with required point density and low point density
point_density = floor(20*((grid_size^2)/(0.3^2)));
[grid_indices_with_required_point_density, ~] = ...
    fcn_findEdge_classifyGridsBasedOnDensity(grids_greater_than_zero_points,total_N_points_in_each_grid,point_density,gridCenters,[],[],[]);

[total_scan_lines_in_each_grid_with_more_than_zero_points] = ....
    fcn_findEdge_calcNumberOfGridScanLines(gridIndices,LIDAR_scanLines,grids_greater_than_zero_points);
[grid_indices_with_more_than_one_scan_line, ~] = ...
    fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters,[],[],[]);

% Finding the orthogonal distances of points in the remaining rings by projecting orthogonally from one ring 
[transverse_span_threshold,transverse_span_each_grid] = fcn_findEdge_determineTransverseSpanThreshold...
    (grids_greater_than_zero_points, grid_AABBs, grid_size, gridIndices, input_points, LIDAR_scanLines,...
    gridCenters_greater_than_zero_point_density,[],[],[]);

% Grids with required point density and low point density
[grid_indices_with_more_than_transverse_span_threshold, ~] =...
    fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid,transverse_span_threshold,grids_greater_than_zero_points,gridCenters,[],[],[]);

% classify Qualified Grids
[current_qualified_grids,~,original_qualified_grids,gridCenters_qualified_grids,~] = ...
    fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,...
    grid_indices_with_more_than_transverse_span_threshold,...
    grids_greater_than_zero_points,gridCenters,[],[],[]);

% "inpolygon" is used to find the grids within the boundary points
shift = 5;
Nscans = length(VehiclePose(:,1));
boundary_points_driven_path = fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose, scanLineRange, Nscans, shift, ([]), ([]));
[in_qg,~] = inpolygon(gridCenters_qualified_grids(:,1),gridCenters_qualified_grids(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));

% Original grid numbers of driven path
% original_grid_numbers_of_driven_path = original_qualified_grids(in_qg);

% Current grid numbers in driven path
current_grid_numbers_of_driven_path = current_qualified_grids(in_qg); %find(in);
gridCenters_driven_path = [gridCenters_qualified_grids(in_qg,1),gridCenters_qualified_grids(in_qg,2)];

% Qualified grid conditions - angle deviation
[angle_btw_unit_normals_and_vertical, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path,...
    angle_btw_unit_normals_and_vertical_driven_path] = fcn_findEdge_determineAngleDeviation...
    (LiDAR_allPoints, gridIndices_cell_array, original_qualified_grids,...
    gridCenters_qualified_grids,current_qualified_grids,gridCenters_driven_path, ...
    grid_AABBs, grid_size, gridIndices, current_grid_numbers_of_driven_path, fig_num,fig_num2,fig_num3);

assert(~isempty(angle_btw_unit_normals_and_vertical));
assert(~isempty(mean_angle_btw_unit_normals_and_vertical_driven_path));
assert(~isempty(angle_btw_unit_normals_and_vertical_driven_path));

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
temp_h2 = figure(fig_num2);
assert(~isempty(get(temp_h2,'Children')))
temp_h3 = figure(fig_num3);
assert(~isempty(get(temp_h3,'Children')))