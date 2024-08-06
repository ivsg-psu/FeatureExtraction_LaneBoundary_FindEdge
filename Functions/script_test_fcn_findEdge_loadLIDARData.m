% script_test_fcn_findEdge_loadLIDARData
% Exercises the function: fcn_findEdge_loadLIDARData

% Revision history:
% 2024_07_25 - Jiabao Zhao 
% -- wrote the code
% 2024_07_31 - S. Brennan
% -- rewrote entire code


%% Test 1 simple load using only defaults, no plotting
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData;


% Check sizes
assert(isequal(length(VehiclePose(:,1)),length(LiDAR_Scan_ENU_Entire_Loop)));
assert(isequal(length(VehiclePose(1,:)),6)); % XYZRPY (roll pitch yaw)
assert(isequal(length(LiDAR_Scan_ENU_Entire_Loop{1}(1,:)),6)); % XYZ intensity scanLine deltaT


%% Test 2: simple load using only defaults via blank entries, no plotting
fig_num = []; 

test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];

[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num));


% Check sizes
assert(isequal(length(VehiclePose(:,1)),length(LiDAR_Scan_ENU_Entire_Loop)));
assert(isequal(length(VehiclePose(1,:)),6)); % XYZRPY (roll pitch yaw)
assert(isequal(length(LiDAR_Scan_ENU_Entire_Loop{1}(1,:)),6)); % XYZ intensity scanLine deltaT

%% Test 3: simple load using only defaults via blank entries, with plotting
fig_num = 1; 

test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];

[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num));


% Check sizes
assert(isequal(length(VehiclePose(:,1)),length(LiDAR_Scan_ENU_Entire_Loop)));
assert(isequal(length(VehiclePose(1,:)),6)); % XYZRPY (roll pitch yaw)
assert(isequal(length(LiDAR_Scan_ENU_Entire_Loop{1}(1,:)),6)); % XYZ intensity scanLine deltaT

%% Test 4: forced load
fig_num = 1; 

test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = 1;

[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num));


% Check sizes
assert(isequal(length(VehiclePose(:,1)),length(LiDAR_Scan_ENU_Entire_Loop)));
assert(isequal(length(VehiclePose(1,:)),6)); % XYZRPY (roll pitch yaw)
assert(isequal(length(LiDAR_Scan_ENU_Entire_Loop{1}(1,:)),6)); % XYZ intensity scanLine deltaT

%% Test 5: specify all inputs
fig_num = 1; 

test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];

[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num));


% Check sizes
assert(isequal(length(VehiclePose(:,1)),length(LiDAR_Scan_ENU_Entire_Loop)));
assert(isequal(length(VehiclePose(1,:)),6)); % XYZRPY (roll pitch yaw)
assert(isequal(length(LiDAR_Scan_ENU_Entire_Loop{1}(1,:)),6)); % XYZ intensity scanLine deltaT

%% Speed test

test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];


% Perform the calculation in slow mode
fig_num = [];
REPS = 5; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num));

    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf;
tic;
for i=1:REPS
    tstart = tic;
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num));

    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fcn_findEdge_loadLIDARData without speed setting (slow) and with speed setting (fast):\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Fail conditions
if 1==0
    % FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§






