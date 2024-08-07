% script_test_fcn_findEdge_plotVehicleLLA
% Exercises the function: fcn_findEdge_plotVehicleLLA

% Revision history:
% 2024_08_07 - S. Brennan
% -- wrote code


%% Test 1 - use empty defaults, except for figure number
fig_num = 1;

% Load data
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, ~] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

reference_LLA = [];
zoom_in_location = [];
zoomLevel = [];
LLA_VehiclePose = fcn_findEdge_plotVehicleLLA(VehiclePose, (reference_LLA), (zoom_in_location), (zoomLevel), (fig_num));


% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

% Check that the data came out the right length
assert(isequal(length(LLA_VehiclePose(:,1)), length(VehiclePose(:,1))));


%% Test 2: call function with NO plotting
fig_num = 2;
figure(fig_num);
clf;


% Load data
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, ~] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

reference_LLA = [];
zoom_in_location = [];
zoomLevel = [];
LLA_VehiclePose = fcn_findEdge_plotVehicleLLA(VehiclePose, (reference_LLA), (zoom_in_location), (zoomLevel), ([]));

% Check that the figure plotted
temp_h = figure(fig_num);
assert(isempty(get(temp_h,'Children')))

% Check that the data came out the right length
assert(isequal(length(LLA_VehiclePose(:,1)), length(VehiclePose(:,1))));

close(fig_num);

%% Test 3: call function with specific reference_LLA
fig_num = 3;
figure(fig_num);
clf;


% Load data
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, ~] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

reference_LLA = [40.862686255557058 -77.834608817369883 0];
zoom_in_location = [];
zoomLevel = [];
LLA_VehiclePose = fcn_findEdge_plotVehicleLLA(VehiclePose, (reference_LLA), (zoom_in_location), (zoomLevel), (fig_num));


% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

% Check that the data came out the right length
assert(isequal(length(LLA_VehiclePose(:,1)), length(VehiclePose(:,1))));


%% Test 4: call function with specific zoom_in_location
fig_num = 4;
figure(fig_num);
clf;


% Load data
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, ~] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

reference_LLA = [];
zoom_in_location = [40.862686255557058 -77.834608817369883];
zoomLevel = [];
LLA_VehiclePose = fcn_findEdge_plotVehicleLLA(VehiclePose, (reference_LLA), (zoom_in_location), (zoomLevel), (fig_num));


% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

% Check that the data came out the right length
assert(isequal(length(LLA_VehiclePose(:,1)), length(VehiclePose(:,1))));

%% Test 5: call function with specific zoomLevel
fig_num = 5;
figure(fig_num);
clf;


% Load data
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, ~] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

reference_LLA = [];
zoom_in_location = [40.862686255557058 -77.834608817369883];
zoomLevel = 20.5;
LLA_VehiclePose = fcn_findEdge_plotVehicleLLA(VehiclePose, (reference_LLA), (zoom_in_location), (zoomLevel), (fig_num));


% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

% Check that the data came out the right length
assert(isequal(length(LLA_VehiclePose(:,1)), length(VehiclePose(:,1))));
%% Speed test

% Load data
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
zoomLevel = [];
[VehiclePose, ~] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string),(flag_load_all_data), (-1));

reference_LLA = [];
zoom_in_location = [];



% Perform the calculation in slow mode
fig_num = [];
REPS = 5; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    LLA_VehiclePose = fcn_findEdge_plotVehicleLLA(VehiclePose, (reference_LLA), (zoom_in_location), (zoomLevel), (fig_num));

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
    LLA_VehiclePose = fcn_findEdge_plotVehicleLLA(VehiclePose, (reference_LLA), (zoom_in_location), (zoomLevel), (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fcn_findEdge_plotVehicleXYZ without speed setting (slow) and with speed setting (fast):\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

% Check that the data came out the right length
assert(isequal(length(LLA_VehiclePose(:,1)), length(VehiclePose(:,1))));


%% Fail conditions
if 1==0
    % FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_findEdge_fitSlopeInterceptNPoints(points,fig_num);
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






