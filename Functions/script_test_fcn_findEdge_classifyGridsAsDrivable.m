% script_test_fcn_findEdge_classifyGridsAsDrivable
% Exercises the function: fcn_findEdge_classifyGridsAsDrivable
% REVISION HISTORY: 
% This function was written on 2024_07_15 by Aneesh Batchu
% -- Seperated this code from fcn_findEdge_classifyGridsAsDrivable
% Funclionalize this code on 7/16/2024 by Jiabao Zhao
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it
% 2024_08_15 - Jiabao Zhao
% -- rewrote the code

%% Test 1
% figure number
fig_num = 1;
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

% Determined std_threshold
std_threshold = 0.1; 
input_points = LiDAR_allPoints(:,1:3); 

% Qualified grid conditions - angle deviation
[angle_btw_unit_normals_and_vertical, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path,...
    angle_btw_unit_normals_and_vertical_driven_path] = fcn_findEdge_determineAngleDeviation...
    (LiDAR_allPoints, gridIndices_cell_array, original_qualified_grids,...
    gridCenters_qualified_grids,current_qualified_grids,gridCenters_driven_path, ...
    grid_AABBs, grid_size, gridIndices, current_grid_numbers_of_driven_path, -1,-1,-1);

% Histogram of angle deviation
theta_threshold = fcn_findEdge_histogramAngleDeviation(angle_btw_unit_normals_and_vertical, ...
    angle_btw_unit_normals_and_vertical_driven_path, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path, -1);

% classify mapped grids into drivable and non-drivable
[standard_deviation_in_z, ... 
    angle_btw_unit_normals_and_vertical, ...
    original_drivable_grids, ...
    original_non_drivable_grids, ...
    current_drivable_grid_numbers_in_mapped_grids, ...
    current_non_drivable_grid_numbers_in_mapped_grids, ...
    current_failed_grid_numbers_in_mapped_grids, ...
    current_uncertain_grid_numbers_in_mapped_grids, ...
    gridCenters_failed_grids, ... 
    ~, ...
    gridCenters_drivable_grids, ...
    gridCenters_non_drivable_grids, ... 
    concatenate_gridCenters_drivable_non_drivable_grids] = ...
    fcn_findEdge_classifyGridsAsDrivable(gridIndices_cell_array, original_qualified_grids, input_points, std_threshold, theta_threshold, gridCenters, fig_num);

assert(~isempty(standard_deviation_in_z))
assert(~isempty(angle_btw_unit_normals_and_vertical))
assert(~isempty(original_drivable_grids))
assert(~isempty(original_non_drivable_grids))
assert(~isempty(current_drivable_grid_numbers_in_mapped_grids))
assert(~isempty(current_non_drivable_grid_numbers_in_mapped_grids))
assert(~isempty(current_failed_grid_numbers_in_mapped_grids))
assert(~isempty(current_uncertain_grid_numbers_in_mapped_grids))
assert(~isempty(gridCenters_failed_grids))
assert(~isempty(gridCenters_drivable_grids))
assert(~isempty(gridCenters_non_drivable_grids))
assert(~isempty(concatenate_gridCenters_drivable_non_drivable_grids))
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% Test 2 no plotting
% figure number
fig_num = 1;
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

% Determined std_threshold
std_threshold = 0.1; 
input_points = LiDAR_allPoints(:,1:3); 

% Qualified grid conditions - angle deviation
[angle_btw_unit_normals_and_vertical, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path,...
    angle_btw_unit_normals_and_vertical_driven_path] = fcn_findEdge_determineAngleDeviation...
    (LiDAR_allPoints, gridIndices_cell_array, original_qualified_grids,...
    gridCenters_qualified_grids,current_qualified_grids,gridCenters_driven_path, ...
    grid_AABBs, grid_size, gridIndices, current_grid_numbers_of_driven_path, -1,-1,-1);

% Histogram of angle deviation
theta_threshold = fcn_findEdge_histogramAngleDeviation(angle_btw_unit_normals_and_vertical, ...
    angle_btw_unit_normals_and_vertical_driven_path, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path, -1);

% classify mapped grids into drivable and non-drivable
[standard_deviation_in_z, ... 
    angle_btw_unit_normals_and_vertical, ...
    original_drivable_grids, ...
    original_non_drivable_grids, ...
    current_drivable_grid_numbers_in_mapped_grids, ...
    current_non_drivable_grid_numbers_in_mapped_grids, ...
    current_failed_grid_numbers_in_mapped_grids, ...
    current_uncertain_grid_numbers_in_mapped_grids, ...
    gridCenters_failed_grids, ... 
    ~, ...
    gridCenters_drivable_grids, ...
    gridCenters_non_drivable_grids, ... 
    concatenate_gridCenters_drivable_non_drivable_grids] = ...
    fcn_findEdge_classifyGridsAsDrivable(gridIndices_cell_array, original_qualified_grids, input_points, std_threshold, theta_threshold, gridCenters, -1);

assert(~isempty(standard_deviation_in_z))
assert(~isempty(angle_btw_unit_normals_and_vertical))
assert(~isempty(original_drivable_grids))
assert(~isempty(original_non_drivable_grids))
assert(~isempty(current_drivable_grid_numbers_in_mapped_grids))
assert(~isempty(current_non_drivable_grid_numbers_in_mapped_grids))
assert(~isempty(current_failed_grid_numbers_in_mapped_grids))
assert(~isempty(current_uncertain_grid_numbers_in_mapped_grids))
assert(~isempty(gridCenters_failed_grids))
assert(~isempty(gridCenters_drivable_grids))
assert(~isempty(gridCenters_non_drivable_grids))
assert(~isempty(concatenate_gridCenters_drivable_non_drivable_grids))
temp_h = figure(fig_num);
assert(isempty(get(temp_h,'Children')))
close(fig_num)
%% Test 3 different std_threshold = 0.2; 
% figure number
fig_num = 1;
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

% Determined std_threshold
std_threshold = 0.2; 
input_points = LiDAR_allPoints(:,1:3); 

% Qualified grid conditions - angle deviation
[angle_btw_unit_normals_and_vertical, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path,...
    angle_btw_unit_normals_and_vertical_driven_path] = fcn_findEdge_determineAngleDeviation...
    (LiDAR_allPoints, gridIndices_cell_array, original_qualified_grids,...
    gridCenters_qualified_grids,current_qualified_grids,gridCenters_driven_path, ...
    grid_AABBs, grid_size, gridIndices, current_grid_numbers_of_driven_path, -1,-1,-1);

% Histogram of angle deviation
theta_threshold = fcn_findEdge_histogramAngleDeviation(angle_btw_unit_normals_and_vertical, ...
    angle_btw_unit_normals_and_vertical_driven_path, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path, -1);

% classify mapped grids into drivable and non-drivable
[standard_deviation_in_z, ... 
    angle_btw_unit_normals_and_vertical, ...
    original_drivable_grids, ...
    original_non_drivable_grids, ...
    current_drivable_grid_numbers_in_mapped_grids, ...
    current_non_drivable_grid_numbers_in_mapped_grids, ...
    current_failed_grid_numbers_in_mapped_grids, ...
    current_uncertain_grid_numbers_in_mapped_grids, ...
    gridCenters_failed_grids, ... 
    gridCenters_uncertain_grids, ...
    gridCenters_drivable_grids, ...
    gridCenters_non_drivable_grids, ... 
    concatenate_gridCenters_drivable_non_drivable_grids] = ...
    fcn_findEdge_classifyGridsAsDrivable(gridIndices_cell_array, original_qualified_grids, input_points, std_threshold, theta_threshold, gridCenters, fig_num);

assert(~isempty(standard_deviation_in_z))
assert(~isempty(angle_btw_unit_normals_and_vertical))
assert(~isempty(original_drivable_grids))
assert(~isempty(original_non_drivable_grids))
assert(~isempty(current_drivable_grid_numbers_in_mapped_grids))
assert(~isempty(current_non_drivable_grid_numbers_in_mapped_grids))
assert(~isempty(current_failed_grid_numbers_in_mapped_grids))
assert(~isempty(current_uncertain_grid_numbers_in_mapped_grids))
assert(~isempty(gridCenters_failed_grids))
assert(~isempty(gridCenters_drivable_grids))
assert(~isempty(gridCenters_non_drivable_grids))
assert(~isempty(concatenate_gridCenters_drivable_non_drivable_grids))
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))