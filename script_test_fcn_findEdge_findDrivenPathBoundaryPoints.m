%script_test_fcn_findEdge_findDrivenPathBoundaryPoints
% Exercises the function: fcn_findEdge_findDrivenPathBoundaryPoints
%
% Revision History:
% 2024_08_14 -Jiabao Zhao
% -- created the test script originally 
%% Test 1 

% load the data
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
fig_num = [];
[VehiclePose, ~] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num));

% get the scanlineRange
station_1 = 1400; 
station_2 = 1450;
range_of_LiDAR = 10;
[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);  
scanLineRange = [station1_minus_range_index station2_plus_range_index]; 

% number of scan lines
Nscans = length(VehiclePose(:,1));

% scalar for shifting the boundar points
shift = 5;

% Find the driven path
fig_num = 1;
fig_num2 = 11;
boundary_points_driven_path = fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose, scanLineRange, Nscans, shift, (fig_num), (fig_num2));

assert(isequal(length(boundary_points_driven_path(1,:)),3))
% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
temp_h2 = figure(fig_num2);
assert(~isempty(get(temp_h2,'Children')))
%% Test 2 no plotting

clf
% load the data
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
fig_num = [];
[VehiclePose, ~] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num));

% get the scanlineRange
station_1 = 1400; 
station_2 = 1450;
range_of_LiDAR = 10;
[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);  
scanLineRange = [station1_minus_range_index station2_plus_range_index]; 

% number of scan lines
Nscans = length(VehiclePose(:,1));

% scalar for shifting the boundar points
shift = 5;

% Find the driven path
fig_num = 2;
fig_num2 = 22;
boundary_points_driven_path = fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose, scanLineRange, Nscans, shift, [], []);

assert(isequal(length(boundary_points_driven_path(1,:)),3))

% Check that the figure plotted
temp_h = figure(fig_num);
temp_h2 = figure(fig_num2);
assert(isempty(get(temp_h,'Children')))
assert(isempty(get(temp_h2,'Children')))
close(fig_num)
close(fig_num2)


%% Test 3 different scanline range

% load the data
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
fig_num = [];
[VehiclePose, ~] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num));

% get the scanlineRange
station_1 = 1400; 
station_2 = 1470;
range_of_LiDAR = 10;
[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);  
scanLineRange = [station1_minus_range_index station2_plus_range_index]; 

% number of scan lines
Nscans = length(VehiclePose(:,1));

% scalar for shifting the boundar points
shift = 5;

% Find the driven path
fig_num = 3;
fig_num2 = 33;
boundary_points_driven_path = fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose, scanLineRange, Nscans, shift, fig_num, fig_num2);

assert(isequal(length(boundary_points_driven_path(1,:)),3))
% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
temp_h2 = figure(fig_num2);
assert(~isempty(get(temp_h2,'Children')))

%% Test 4, different shift scalar  

% load the data
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
fig_num = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num));

% get the scanlineRange
station_1 = 1400; 
station_2 = 1470;
range_of_LiDAR = 10;
[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);  
scanLineRange = [station1_minus_range_index station2_plus_range_index]; 

% number of scan lines
Nscans = length(VehiclePose(:,1));

% scalar for shifting the boundar points
shift = 100;

% Find the driven path
fig_num = 4;
fig_num2 = 44;
boundary_points_driven_path = fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose, scanLineRange, Nscans, shift, fig_num, fig_num2);

assert(isequal(length(boundary_points_driven_path(1,:)),3))
% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
temp_h2 = figure(fig_num2);
assert(~isempty(get(temp_h2,'Children')))