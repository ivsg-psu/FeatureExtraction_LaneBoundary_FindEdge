% script_test_fcn_findEdge_determineTransverseSpanThreshold
% Exercises the function: fcn_findEdge_determineTransverseSpanThreshold
% Revision history:
% 2024_08_15 - Jiabao Zhao
% -- wrote the code originally 

%% Test 1 
% load the data
fig_num = 1; 
fig_num2 = 11;
fig_num3 = 111;
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), [], []);


% Extract scan lines
station_1 = 1400; 
station_2 = 1450;
range_of_LiDAR = 10;
[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);  
scanLineRange = [station1_minus_range_index station2_plus_range_index]; 
ringsRange = []; 
[~, ~, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]), ([]));
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
[~, ~, ~, ~,...
    grids_greater_than_zero_points, ~,...
    gridCenters_greater_than_zero_point_density, gridIndices, grid_AABBs] = fcn_findEdge_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,[]);

[transverse_span_threshold,transverse_span_each_grid] = fcn_findEdge_determineTransverseSpanThreshold...
    (grids_greater_than_zero_points, grid_AABBs, grid_size, gridIndices, input_points, LIDAR_scanLines,...
    gridCenters_greater_than_zero_point_density,fig_num,fig_num2,fig_num3);

assert(isequal(transverse_span_threshold,0.15));
assert(isequal(length(transverse_span_each_grid(1,:)),1));

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
temp_h2 = figure(fig_num2);
assert(~isempty(get(temp_h2,'Children')))
temp_h3 = figure(fig_num3);
assert(~isempty(get(temp_h3,'Children')))
%%  Test 2 different grid size
 
% load the data
fig_num = 2; 
fig_num2 = 22;
fig_num3 = 222;
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), [], []);


% Extract scan lines
station_1 = 1400; 
station_2 = 1450;
range_of_LiDAR = 10;
[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);  
scanLineRange = [station1_minus_range_index station2_plus_range_index]; 
ringsRange = []; 
[~, ~, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]), ([]));
[~, ~, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2,LIDAR_intensity,[]);


% prepare for grids 
LiDAR_allPoints = [LIDAR_ENU(in_domain,:), LIDAR_scanLineAndRingID(in_domain,:)];
LiDAR_allPoints = LiDAR_allPoints(~isnan(LiDAR_allPoints(:,1)),:);
LIDAR_scanLines = LiDAR_allPoints(:,4:5); 
input_points = LiDAR_allPoints(:,1:2); 
grid_size = 10; %0.8;%1;%1.25; %1.26
% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,~,~] = fcn_findEdge_findMaxMinOfXYZ(input_points,-1);
% The 2D grid boundaries required for the current analysis
grid_boundaries = [Min_x Max_x Min_y Max_y]; 

% Seperate the data into grids
[~, ~, ~, ~,...
    grids_greater_than_zero_points, ~,...
    gridCenters_greater_than_zero_point_density, gridIndices, grid_AABBs] = fcn_findEdge_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,[]);

[transverse_span_threshold,transverse_span_each_grid] = fcn_findEdge_determineTransverseSpanThreshold...
    (grids_greater_than_zero_points, grid_AABBs, grid_size, gridIndices, input_points, LIDAR_scanLines,...
    gridCenters_greater_than_zero_point_density,fig_num,fig_num2,fig_num3);

assert(isequal(transverse_span_threshold,0.15));
assert(isequal(length(transverse_span_each_grid(1,:)),1));

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
temp_h2 = figure(fig_num2);
assert(~isempty(get(temp_h2,'Children')))
temp_h3 = figure(fig_num3);
assert(~isempty(get(temp_h3,'Children')))
%%  Test 3 no plotting
 
% load the data
fig_num = 3; 
fig_num2 = 33;
fig_num3 = 333;
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), [], []);


% Extract scan lines
station_1 = 1400; 
station_2 = 1450;
range_of_LiDAR = 10;
[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);  
scanLineRange = [station1_minus_range_index station2_plus_range_index]; 
ringsRange = []; 
[~, ~, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]), ([]));
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
[~, ~, ~, ~,...
    grids_greater_than_zero_points, ~,...
    gridCenters_greater_than_zero_point_density, gridIndices, grid_AABBs] = fcn_findEdge_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,[]);

[transverse_span_threshold,transverse_span_each_grid] = fcn_findEdge_determineTransverseSpanThreshold...
    (grids_greater_than_zero_points, grid_AABBs, grid_size, gridIndices, input_points, LIDAR_scanLines,...
    gridCenters_greater_than_zero_point_density,[],[],[]);

assert(isequal(transverse_span_threshold,0.15));
assert(isequal(length(transverse_span_each_grid(1,:)),1));

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