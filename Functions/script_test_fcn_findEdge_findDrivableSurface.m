%script_test_fcn_findEdge_findDrivableSurface
% Exercises the function: fcn_findEdge_findDrivableSurface
%
% Revision History:
% 2024_08_11 - Jiabao Zhao
% -- Wrote the script
%% Test 1 

% Load data - Loads in VehiclePose and LiDAR_Scan_ENU_Entire_Loop
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

% Set defaults for which scans to extract
scanLineRange = [1400 1450];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
    LIDAR_ENU, ~, ~] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));

[LIDAR_ENU_under_vehicle] = fcn_findEdge_findDrivableSurface (LIDAR_ENU, VehiclePose_ENU, VehiclePose_UnitOrthoVectors);
assert(isequal(length(LIDAR_ENU_under_vehicle(1,:)),3));


%% Test 2, with different scanLineRange
scaling=[];
color_map=[];
marker_size=[];
reference_LLA=[];


% Load data - Loads in VehiclePose and LiDAR_Scan_ENU_Entire_Loop
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

% Set defaults for which scans to extract
scanLineRange = [1400 1480];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
    LIDAR_ENU,  ~, ~] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));

[LIDAR_ENU_under_vehicle] = fcn_findEdge_findDrivableSurface (LIDAR_ENU, VehiclePose_ENU, VehiclePose_UnitOrthoVectors);
assert(isequal(length(LIDAR_ENU_under_vehicle(1,:)),3));