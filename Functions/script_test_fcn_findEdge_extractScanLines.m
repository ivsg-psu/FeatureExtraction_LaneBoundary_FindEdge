% script_test_fcn_findEdge_extractScanLines
% Exercises the function: fcn_findEdge_extractScanLines

% Revision history:
% 2024_08_06 - S. Brennan
% -- wrote code


%% Test 1: an example load
fig_num = 1;
figure(fig_num);
clf;

% Load data
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

% Set defaults
scanLineRange = [1400 1450];
ringsRange = [];

% Extract scan lines
[VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (fig_num));

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

% Check that the points are filled to same sizes
Npoints = length(VehiclePose_ENU(:,1));
assert(isequal(size(VehiclePose_ENU),[Npoints 3]));
assert(isequal(size(VehiclePose_UnitOrthoVectors),[Npoints 3]));
assert(isequal(size(LIDAR_ENU),[Npoints 3]));
assert(isequal(size(LIDAR_intensity),[Npoints 1]));
assert(isequal(size(LIDAR_scanLineAndRingID),[Npoints 2]));

% Check key values
assert(min(LIDAR_scanLineAndRingID(:,1))==scanLineRange(1))
assert(max(LIDAR_scanLineAndRingID(:,1))==scanLineRange(2))
assert(min(LIDAR_scanLineAndRingID(:,2))==0)
assert(max(LIDAR_scanLineAndRingID(:,2))==15)

%% Test 2: call function with NO plotting
fig_num = 2;
figure(fig_num);
clf;

% Load data
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

% Set defaults
scanLineRange = [1400 1450];
ringsRange = [];

% Extract scan lines
[VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));

% Check that the figure is NOT plotted
temp_h = figure(fig_num);
assert(isempty(get(temp_h,'Children')))
close(fig_num);

% Check that the points are filled to same sizes
Npoints = length(VehiclePose_ENU(:,1));
assert(isequal(size(VehiclePose_ENU),[Npoints 3]));
assert(isequal(size(VehiclePose_UnitOrthoVectors),[Npoints 3]));
assert(isequal(size(LIDAR_ENU),[Npoints 3]));
assert(isequal(size(LIDAR_intensity),[Npoints 1]));
assert(isequal(size(LIDAR_scanLineAndRingID),[Npoints 2]));

% Check key values
assert(min(LIDAR_scanLineAndRingID(:,1))==scanLineRange(1))
assert(max(LIDAR_scanLineAndRingID(:,1))==scanLineRange(2))
assert(min(LIDAR_scanLineAndRingID(:,2))==0)
assert(max(LIDAR_scanLineAndRingID(:,2))==15)


%% Speed test

% Load data
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

% Set defaults
scanLineRange = [1400 1450];
ringsRange = [];


% Perform the calculation in slow mode
fig_num = [];
REPS = 5; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (fig_num));

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
    [VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fcn_findEdge_plotVehicleXY without speed setting (slow) and with speed setting (fast):\n');
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