% script_test_fcn_findEdge_plotLIDARXY
% Exercises the function: fcn_findEdge_plotLIDARXY

% Revision History:
% 2024_08_06 - Aleksandr Goncharov
% -- Wrote the script
% 2024_08_07 - S. Brennan
% -- Cleaned up comments. Script looks good.


%% TEST 1 - No Optional Inputs Besides Figure Number
% INPUTS
fig_num=1;
color_triplet=[];
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
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (-1));

% Plot LIDARXY

fcn_findEdge_plotLIDARXY(LIDAR_ENU, (color_triplet), (marker_size), (fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% TEST 2 - No Optional Inputs - No Figure
% INPUTS
fig_num=2;
color_triplet=[];
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
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (-1));

% Plot LIDARXY

fcn_findEdge_plotLIDARXY(LIDAR_ENU, (color_triplet), (marker_size), [])

% Check that the figure plotted
temp_h = figure(fig_num);
assert(isempty(get(temp_h,'Children')))
close(fig_num);

%% TEST 3 - Color Triplet Changed to Red
% INPUTS
fig_num=3;
color_triplet=[1 0 0];
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
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (-1));

% Plot LIDARXY

fcn_findEdge_plotLIDARXY(LIDAR_ENU, (color_triplet), (marker_size), (fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% TEST 4 - Marker Size changed to 20, everything else default
% INPUTS
fig_num=4;
color_triplet=[];
marker_size=20;


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
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (-1));

% Plot LIDARXY

fcn_findEdge_plotLIDARXY(LIDAR_ENU, (color_triplet), (marker_size), (fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% TEST 5 - Both Size and Color Changed, Green markers size 10.
% INPUTS
fig_num=5;
color_triplet=[0 1 0];
marker_size=10;


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
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (-1));

% Plot LIDARXY

fcn_findEdge_plotLIDARXY(LIDAR_ENU, (color_triplet), (marker_size), (fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% TEST 6 - SPEED TEST

%Optional Inputs
color_triplet=[];
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
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (-1));

% Speed Test Calculation
fig_num=[];
REPS=5; minTimeSlow=Inf;
tic;

%slow mode calculation - code copied from plotVehicleXYZ

for i=1:REPS
    tstart=tic;
    
    fcn_findEdge_plotLIDARXY(LIDAR_ENU, (color_triplet), (marker_size), (fig_num))
    
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

 fcn_findEdge_plotLIDARXY(LIDAR_ENU, (color_triplet), (marker_size), (fig_num))
   
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
% assert(averageTimeFast<averageTimeSlow*1.2);
