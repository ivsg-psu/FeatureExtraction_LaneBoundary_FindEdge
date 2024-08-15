%script_test_fcn_findEdge_findDrivableSurface
% Exercises the function: fcn_findEdge_findDrivableSurface
%
% Revision History:
% 2024_08_11 - Jiabao Zhao
% -- Wrote the script
%% Test 1 
fig_num = 1;
fig_num2 = 11;
% Load data - Loads in VehiclePose and LiDAR_Scan_ENU_Entire_Loop
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), [], []);

% Set defaults for which scans to extract
scanLineRange = [1400 1450];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
    LIDAR_ENU, ~, ~] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), [],[]);

[LIDAR_ENU_under_vehicle] = fcn_findEdge_findDrivableSurface (LIDAR_ENU, VehiclePose_ENU, VehiclePose_UnitOrthoVectors,fig_num,fig_num2);

assert(isequal(length(LIDAR_ENU_under_vehicle(1,:)),3));
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
temp_h2 = figure(fig_num2);
assert(~isempty(get(temp_h2,'Children')))

%% Test 2, with different scanLineRange
fig_num = 2;
fig_num2 = 22;
% Load data - Loads in VehiclePose and LiDAR_Scan_ENU_Entire_Loop
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), [], []);

% Set defaults for which scans to extract
scanLineRange = [1400 1480];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
    LIDAR_ENU, ~, ~] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), [],[]);

[LIDAR_ENU_under_vehicle] = fcn_findEdge_findDrivableSurface (LIDAR_ENU, VehiclePose_ENU, VehiclePose_UnitOrthoVectors,fig_num,fig_num2);

assert(isequal(length(LIDAR_ENU_under_vehicle(1,:)),3));
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
temp_h2 = figure(fig_num2);
assert(~isempty(get(temp_h2,'Children')))

%% Test 3 no plotting 
fig_num = 3;
fig_num2 = 33;
% Load data - Loads in VehiclePose and LiDAR_Scan_ENU_Entire_Loop
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), [], []);

% Set defaults for which scans to extract
scanLineRange = [1400 1480];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
    LIDAR_ENU, ~, ~] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), [],[]);

[LIDAR_ENU_under_vehicle] = fcn_findEdge_findDrivableSurface (LIDAR_ENU, VehiclePose_ENU, VehiclePose_UnitOrthoVectors,[],[]);

assert(isequal(length(LIDAR_ENU_under_vehicle(1,:)),3));
temp_h = figure(fig_num);
assert(isempty(get(temp_h,'Children')))
temp_h2 = figure(fig_num2);
assert(isempty(get(temp_h2,'Children')))
close(fig_num)
close(fig_num2)
