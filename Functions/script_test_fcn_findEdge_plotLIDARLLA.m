%script_test_fcn_findEdge_plotLIDARLLA
% Exercises the function: fcn_findEdge_plotLIDARLLA
%
% Revision History:
% 2024_08_06 - Aleksandr Goncharov
% -- Wrote the script
% 2024_08_07 -Jiabao Zhao
% -- reordered and simplified the inputs, allowing variable input arguments
% -- changed the input of each script and added another test script for LLA

%% Test 1 - No Optional Inputs Besides Fig_Num
%Optional Inputs
fig_num=1;
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
scanLineRange = [1400 1450];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));


%plotLIDARLLA
LIDAR_intensity = [];
format = [];
fcn_findEdge_plotLIDARLLA(LIDAR_ENU,(LIDAR_intensity),(scaling),(color_map),(marker_size),(reference_LLA), (format),(fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% Test 2 - No Optional Inputs - No Plot
%Optional Inputs
fig_num=2;
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
scanLineRange = [1400 1450];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));


%plotLIDARLLA
LIDAR_intensity = [];
format = [];
fcn_findEdge_plotLIDARLLA(LIDAR_ENU,(LIDAR_intensity),(scaling),(color_map),(marker_size),(reference_LLA), (format),(-1))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(isempty(get(temp_h,'Children')))
close(fig_num)

%% Test 3 - add LIDAR_intensity input
%Optional Inputs
fig_num=3;
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
scanLineRange = [1400 1450];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[~, ~, LIDAR_ENU, LIDAR_intensity, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));


%plotLIDARLLA

format = [];
fcn_findEdge_plotLIDARLLA(LIDAR_ENU,(LIDAR_intensity),(scaling),(color_map),(marker_size),(reference_LLA), (format),(fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% Test 4 - Scaling Changed to 5 everything else default
%Optional Inputs
fig_num=4;
scaling=5;
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
scanLineRange = [1400 1450];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));


%plotLIDARLLA
%plotLIDARLLA
LIDAR_intensity = [];
format = [];
fcn_findEdge_plotLIDARLLA(LIDAR_ENU,(LIDAR_intensity),(scaling),(color_map),(marker_size),(reference_LLA), (format),(fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% Test 5- Colormap changed to 'autumn' with LIDAR_intensity, rest is default
%Optional Inputs
fig_num=5;
scaling=[];
color_map='autumn';
marker_size=[];
reference_LLA=[];

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
[~, ~, LIDAR_ENU, LIDAR_intensity, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));


%plotLIDARLLA
format = [];
fcn_findEdge_plotLIDARLLA(LIDAR_ENU,(LIDAR_intensity),(scaling),(color_map),(marker_size),(reference_LLA), (format),(fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% Test 6 - Marker Size changed to 10, rest is default
%Optional Inputs
fig_num=6;
scaling=[];
color_map=[];
marker_size=10;
reference_LLA=[];

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
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));


%plotLIDARLLA
format = [];
LIDAR_intensity = [];
fcn_findEdge_plotLIDARLLA(LIDAR_ENU,(LIDAR_intensity),(scaling),(color_map),(marker_size),(reference_LLA), (format),(fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% Test 7 add reference_LLA, rest is default

fig_num=7;
scaling=[];
color_map=[];
marker_size=[];
reference_LLA = [40.862686255557058 -77.834608817369883 0];

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
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));


%plotLIDARLLA
format = [];
LIDAR_intensity = [];
fcn_findEdge_plotLIDARLLA(LIDAR_ENU,(LIDAR_intensity),(scaling),(color_map),(marker_size), (reference_LLA), (format), (fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% Test 8 add format string

fig_num=7;
scaling=[];
color_map=[];
marker_size=[];
reference_LLA = [];

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
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));


%plotLIDARLLA
format = sprintf(' ''.'',''Color'',[0 0 1],''MarkerSize'', 20');
LIDAR_intensity = [];
fcn_findEdge_plotLIDARLLA(LIDAR_ENU,(LIDAR_intensity),(scaling),(color_map),(marker_size), (reference_LLA), (format), (fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
%% Test 8 - Speed Test Default Inputs
% NOTE: Assertion sometimes fails due to average time of fast mode being
% sometimes slower than average time of regular runtime due to variance
%
%Optional Inputs
scaling=[];
color_map=[];
marker_size=[];

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
[~, ~, LIDAR_ENU, LIDAR_intensity, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));

% Load data - GPS Object test track reference

reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;

gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class

%
%
% Speed Test Calculation
%
%
fig_num=[];
REPS=5; minTimeSlow=Inf;
tic;

%slow mode calculation - code copied from plotVehicleXYZ

for i=1:REPS
    tstart=tic;
    
    fcn_findEdge_plotLIDARLLA(LIDAR_ENU,(LIDAR_intensity),(scaling),(color_map),(marker_size),(reference_LLA),(format),(fig_num))
    
    telapsed=toc(tstart);
    minTimeSlow=min(telapsed,minTimeSlow);
end
averageTimeSlow=toc/REPS;

%slow mode END

%Fast Mode Calculation
fig_num = -1;
minTimeFast = Inf;
tic;
for i=1:REPS
    tstart = tic;
    
    fcn_findEdge_plotLIDARLLA(LIDAR_ENU,(LIDAR_intensity),(scaling),(color_map),(marker_size),(reference_LLA),(format),(fig_num))
    
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

%Display Console Comparison
if 1==1
fprintf(1,'\n\nComparison of fcn_findEdge_plotLIDARXYZ without speed setting (slow) and with speed setting (fast):\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);
end

%Assertion on averageTime NOTE: Due to the variance, there is a chance that
%the assertion will fail.
assert(averageTimeFast<averageTimeSlow);

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