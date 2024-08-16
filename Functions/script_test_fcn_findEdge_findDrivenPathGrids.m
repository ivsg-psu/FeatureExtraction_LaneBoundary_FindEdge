%script_test_fcn_findEdge_findDrivenPathGrids
% Exercises the function:fcn_findEdge_findDrivenPathGrids
%
% Revision History:
% 2024_08_11 - Jiabao Zhao
% -- Wrote the script
%% Test 1 
% load the data
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), -1, -1);

% initialize scanlineRnage
station_1 = 1400; 
station_2 = 1450;
range_of_LiDAR = 10;
[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);  
scanLineRange = [station1_minus_range_index station2_plus_range_index]; 

Nscans = length(VehiclePose(:,1));
shift = 5;
boundary_points_driven_path = fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose, scanLineRange, Nscans, shift, -1, (-1));

% Nscans = length(VehiclePose(:,1));
% Set defaults for which scans to extract
% scanLineRange = [1400 1450];
ringsRange = []; % If leave empty, it loads all rings
% Extract scan lines
[~, ~, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), -1,-1);

[~, ~, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2,LIDAR_intensity, -1);


% These are concatenated LiDAR points of chosen scans and cells in the
% first step. 
LiDAR_allPoints = [LIDAR_ENU(in_domain,:), LIDAR_scanLineAndRingID(in_domain,:)];
% remove NANs
LiDAR_allPoints = LiDAR_allPoints(~isnan(LiDAR_allPoints(:,1)),:);
% scanLine_rings without NaNs
% LIDAR_scanLines = LiDAR_allPoints(:,4:5); 
% Input points for seperating the data into grids. The points are in 2D as
% the analysis is carried out in 2D
input_points = LiDAR_allPoints(:,1:2); 
% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1; %0.8;%1;%1.25; %1.26
% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,~,~] = fcn_findEdge_findMaxMinOfXYZ(input_points,-1);
% The grids are plotted only within the boundaries in #D
% grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 
% The 2D grid boundaries required for the current analysis
grid_boundaries = [Min_x Max_x Min_y Max_y]; 
% Seperate the data into grids
[~, total_N_points_in_each_grid, ~, ~,...
    grids_greater_than_zero_points, ~,...
    gridCenters_greater_than_zero_point_density, ~, ~] = fcn_findEdge_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,-1);

% find the drivavle path grids
fig_num = 110;
ENU_3D_fig_num = 111;
format = sprintf('''.'',''MarkerSize'',30,''Color'',[0.8 0.8 0.8]');
format1 = sprintf('''.'',''MarkerSize'',30,''Color'',[0 1 0]');
current_grids_greater_than_zero_points = 1:length(grids_greater_than_zero_points); 

[~, ~,total_points_in_each_grid_in_the_driven_path, total_points_in_each_grid_with_points_greater_than_zero]...
    = fcn_findEdge_findDrivenPathGrids(gridCenters_greater_than_zero_point_density, boundary_points_driven_path,...
    grids_greater_than_zero_points, current_grids_greater_than_zero_points, total_N_points_in_each_grid, (format), (format1),[],[], (fig_num), (ENU_3D_fig_num));

assert(~isempty(total_points_in_each_grid_in_the_driven_path))
assert(~isempty(total_points_in_each_grid_with_points_greater_than_zero))
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
temp_h2 = figure(ENU_3D_fig_num);
assert(~isempty(get(temp_h2,'Children')))

%% Test 2 no plotting 

% load the data
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), -1, -1);

% initialize scanlineRnage
station_1 = 1400; 
station_2 = 1450;
range_of_LiDAR = 10;
[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);  
scanLineRange = [station1_minus_range_index station2_plus_range_index]; 

Nscans = length(VehiclePose(:,1));
shift = 5;
boundary_points_driven_path = fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose, scanLineRange, Nscans, shift, -1, (-1));

% Nscans = length(VehiclePose(:,1));
% Set defaults for which scans to extract
% scanLineRange = [1400 1450];
ringsRange = []; % If leave empty, it loads all rings
% Extract scan lines
[~, ~, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), -1,-1);

[~, ~, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2,LIDAR_intensity, -1);


% These are concatenated LiDAR points of chosen scans and cells in the
% first step. 
LiDAR_allPoints = [LIDAR_ENU(in_domain,:), LIDAR_scanLineAndRingID(in_domain,:)];
% remove NANs
LiDAR_allPoints = LiDAR_allPoints(~isnan(LiDAR_allPoints(:,1)),:);
% scanLine_rings without NaNs
% LIDAR_scanLines = LiDAR_allPoints(:,4:5); 
% Input points for seperating the data into grids. The points are in 2D as
% the analysis is carried out in 2D
input_points = LiDAR_allPoints(:,1:2); 
% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1; %0.8;%1;%1.25; %1.26
% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,~,~] = fcn_findEdge_findMaxMinOfXYZ(input_points,-1);
% The grids are plotted only within the boundaries in #D
% grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 
% The 2D grid boundaries required for the current analysis
grid_boundaries = [Min_x Max_x Min_y Max_y]; 
% Seperate the data into grids
[gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, grids_with_zero_points,...
    grids_greater_than_zero_points, gridCenters_zero_point_density,...
    gridCenters_greater_than_zero_point_density, gridIndices, grid_AABBs] = fcn_findEdge_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,-1);

% find the drivavle path grids
fig_num = 110;
ENU_3D_fig_num = 111;
format = sprintf('''.'',''MarkerSize'',30,''Color'',[0.8 0.8 0.8]');
format1 = sprintf('''.'',''MarkerSize'',30,''Color'',[0 1 0]');
current_grids_greater_than_zero_points = 1:length(grids_greater_than_zero_points); 

[~, ~,total_points_in_each_grid_in_the_driven_path, total_points_in_each_grid_with_points_greater_than_zero]...
    = fcn_findEdge_findDrivenPathGrids(gridCenters_greater_than_zero_point_density, boundary_points_driven_path,...
    grids_greater_than_zero_points, current_grids_greater_than_zero_points, total_N_points_in_each_grid, (format), (format1),[],[], (-1), (-1));

assert(~isempty(total_points_in_each_grid_in_the_driven_path))
assert(~isempty(total_points_in_each_grid_with_points_greater_than_zero))
temp_h = figure(fig_num);
assert(isempty(get(temp_h,'Children')))
temp_h2 = figure(ENU_3D_fig_num);
assert(isempty(get(temp_h2,'Children')))
close(fig_num)
close(ENU_3D_fig_num)