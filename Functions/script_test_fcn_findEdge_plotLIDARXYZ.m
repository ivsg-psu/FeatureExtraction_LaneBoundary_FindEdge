% script_test_fcn_findEdge_plotLIDARXYZ
% Exercises the function: fcn_findEdge_plotLIDARXYZ

% Revision History:
% 2024_08_06 - Aleksandr Goncharov
% -- Wrote the script
% 2024_08_07 - Sean Brennan
% -- changed ordering of variable inputs from most important to least
% -- cleaned up comments and details throughout
% 2024_08_11 - Jiabao Zhao
% -- added format string as a optional input

%% Test 1 - No Optional Inputs Besides Fig_Num
%Optional Inputs
fig_num = 1; 

LIDAR_intensity = [];
scaling=[];
color_map=[];



% Load data - Loads in VehiclePose and LiDAR_Scan_ENU_Entire_Loop
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

% Set defaults for which scans to extract
scanLineRange = [1400 1420];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), -1);

% Call the plotting function
format = [];
fcn_findEdge_plotLIDARXYZ(LIDAR_ENU, (LIDAR_intensity), (scaling), (color_map),(format),(fig_num));

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))


%% Test 2 - No Optional Inputs - No Plots

%Optional Inputs
fig_num = 2; 

LIDAR_intensity = [];
scaling=[];
color_map=[];

% Load data - Loads in VehiclePose and LiDAR_Scan_ENU_Entire_Loop
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

% Set defaults for which scans to extract
scanLineRange = [1400 1420];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), -1);

% Call the plotting function
fcn_findEdge_plotLIDARXYZ(LIDAR_ENU, (LIDAR_intensity), (scaling), (color_map), (format), (-1))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(isempty(get(temp_h,'Children')))
close(fig_num);

%% Test 3 - add intensity input

%Optional Inputs
fig_num = 3; 

scaling=[];
color_map=[];


% Load data - Loads in VehiclePose and LiDAR_Scan_ENU_Entire_Loop
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

% Set defaults for which scans to extract
scanLineRange = [1400 1420];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[~, ~, LIDAR_ENU, LIDAR_intensity, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), -1);

% Call the plotting function
fcn_findEdge_plotLIDARXYZ(LIDAR_ENU, (LIDAR_intensity), (scaling), (color_map), (format), (fig_num));

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% Test 4 - Add color_map as autumn

%Optional Inputs
fig_num = 4; 

scaling=[];
color_map='autumn';

% Load data - Loads in VehiclePose and LiDAR_Scan_ENU_Entire_Loop
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

% Set defaults for which scans to extract
scanLineRange = [1400 1420];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[~, ~, LIDAR_ENU, LIDAR_intensity, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), -1);

% Call the plotting function
format = [];
fcn_findEdge_plotLIDARXYZ(LIDAR_ENU, (LIDAR_intensity), (scaling), (color_map), (format), (fig_num));

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% Test 5 - Add color_map as autumn, change scaling to 10

%Optional Inputs
fig_num = 5; 

scaling=10;
color_map='autumn';

% Load data - Loads in VehiclePose and LiDAR_Scan_ENU_Entire_Loop
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

% Set defaults for which scans to extract
scanLineRange = [1400 1420];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[~, ~, LIDAR_ENU, LIDAR_intensity, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), -1);

% Call the plotting function
format = [];
fcn_findEdge_plotLIDARXYZ(LIDAR_ENU, (LIDAR_intensity), (scaling), (color_map), (format) ,(fig_num));

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% Test 6 - Add format string

%Optional Inputs
fig_num = 6; 

scaling=[];
color_map=[];

% Load data - Loads in VehiclePose and LiDAR_Scan_ENU_Entire_Loop
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

% Set defaults for which scans to extract
scanLineRange = [1400 1420];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), -1);

% Call the plotting function
LIDAR_intensity = [];
format = sprintf(' ''.'',''Color'',[0 0 1],''MarkerSize'', 20');
fcn_findEdge_plotLIDARXYZ(LIDAR_ENU, (LIDAR_intensity), (scaling), (color_map), (format) ,(fig_num));

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% Test 7 - Speed Mode

%Optional Inputs

scaling=[];
color_map=[];



% Load data - Loads in VehiclePose and LiDAR_Scan_ENU_Entire_Loop
test_date_string = [];
vehicle_pose_string = [];
LIDAR_file_string = [];
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (-1));

% Set defaults for which scans to extract
scanLineRange = [1400 1420];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines - Creates a Variable LIDAR_ENU and LIDAR_intensity
[~, ~, LIDAR_ENU, ~, ~] = ...
fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), -1);

% Speed Test Calculation
fig_num=[];
REPS=5; minTimeSlow=Inf;
tic;

%slow mode calculation - code copied from plotVehicleXYZ

for i=1:REPS
    tstart=tic;
    fcn_findEdge_plotLIDARXYZ(LIDAR_ENU, (LIDAR_intensity), (scaling), (color_map), (format), (fig_num));
    
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
    fcn_findEdge_plotLIDARXYZ(LIDAR_ENU, (LIDAR_intensity), (scaling), (color_map), (format), (fig_num));
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

