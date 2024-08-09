%% Introduction to and Purpose of the Find Edge library
% This is a demonstration script to show the primary functionality of the
% FindEdge library. This library finds the road edge from LIDAR data and
% XYZ trajectory data of a mapping vehicle. The road edge is the location
% where the pavement stops being a drivable surface, usually the edge of
% the pavement and the vegetation next to the road.
%
% This script follows the updated surface analysis schema step by step
% following the codes originally developed in:
% script_test_geometry_updatedSurfaceAnalysis
%
% This is the explanation of the code that can be found by running
%       script_demo_FindEdge.m
% 
% This code repo is typically located at:
%
%   https://github.com/ivsg-psu/FeatureExtraction_LaneBoundary_FindEdge
%
% If you have questions or comments, please contact Sean Brennan at
% sbrennan@psu.edu 
% 
% Additional contributers include:
% 2024 - Aneesh Batchu, Jiabao Zhao, and Aleksandr Goncharov


%% Revision History:
% 2024_08_05 - sbrennan@psu.edu
% -- pulled core codes out of Geometry class library into this repo.
% -- Functionalize core data loading steps
% -- pass on to Jiabao
% 2024_08_06 -  jpz5469@psu.edu
% -- document the output in fcn_findEdge_extractScanLines
% -- add a LargeData subdirectory in the main directory. Be sure to
% download large data files from the IVSG OneDrive folder - these files are
% too large to host on GitHub.
% -- Started copying core functions from the Geometry class to findEdge
% -- removed copied functions from the Geometry class.  
% 2024_08_07 -  S. Brennan, sbrennan@psu.edu
% -- cleaned up comments in fcn_findEdge_extractScanLines
% -- added script_test_fcn_findEdge_plotVehicleLLA
% 2024_08_07 -  Jiabao Zhao, jpz5469@psu.edu
% --add test section for plotVehicleLLA to test zoomLevel. Add code in
% --function to implement this as a variable input argument
% --cleaned up more functions and add replace the script with two main
% functions
% 2024_08_08 -  Jiabao Zhao, jpz5469@psu.edu
% -- added format string as a optional input for fcn_findEdge_plotVehicleXY
% and fcn_findEdge_plotLIDARLLA since these two are the "core" function for
% later use. 


%% To-do items
% 2024_08_05 - created by S. Brennan
% -- Aleks and Jiabao - finish funtionalizing starting at step 1 below and
% document each change in the revision history above when done. Update the
% flow chart with functions
% -- Add a LargeData file in the main directory? EVERYONE?\
% -- Should leave this function in the geometry libriary: fcn_geometry_fitPlaneLinearRegression
%
% FOR JIABAO ONLY
% ADD SIMPLE TEST SCRIPT FOR THE FOLLOWING FUNCTIONS: 
% fcn_findEdge_findGridsWithPoints
% fcn_findEdge_plotVehicleLLA
% fcn_findEdge_GridsIntoMappedUnmapped
% fcn_findEdge_findGridsWithPoints\
% fcn_findEdge_plotVehicleLLA (This one needs a new script)
% fcn_geometry_classifyGridsAsDrivable
%
% Delete as many plot commands in the script as we can by just calling the
% "core" plotting commands we already developed.
% add colormap command for fcn_findEdge_plotVehicleXY.

%% Prep the workspace
close all
clc

%% Dependencies and Setup of the Code
% The code requires several other libraries to work, namely the following
%
% * DebugTools - the repo can be found at: https://github.com/ivsg-psu/Errata_Tutorials_DebugTools


% List what libraries we need, and where to find the codes for each
clear library_name library_folders library_url

ith_library = 1;
library_name{ith_library}    = 'DebugTools_v2023_04_22';
library_folders{ith_library} = {'Functions','Data'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/archive/refs/tags/DebugTools_v2023_04_22.zip';

ith_library = ith_library+1;
library_name{ith_library}    = 'PathClass_v2024_03_14';
library_folders{ith_library} = {'Functions'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary/archive/refs/tags/PathClass_v2024_03_14.zip';

ith_library = ith_library+1;
library_name{ith_library}    = 'GPSClass_v2023_06_29';
library_folders{ith_library} = {'Functions'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/FieldDataCollection_GPSRelatedCodes_GPSClass/archive/refs/tags/GPSClass_v2023_06_29.zip';

% ith_library = ith_library+1;
% library_name{ith_library}    = 'LineFitting_v2023_07_24';
% library_folders{ith_library} = {'Functions'};
% library_url{ith_library}     = 'https://github.com/ivsg-psu/FeatureExtraction_Association_LineFitting/archive/refs/tags/LineFitting_v2023_07_24.zip';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'FindCircleRadius_v2023_08_02';
% library_folders{ith_library} = {'Functions'};                                
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_GeomTools_FindCircleRadius/archive/refs/tags/FindCircleRadius_v2023_08_02.zip';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'BreakDataIntoLaps_v2023_08_25';
% library_folders{ith_library} = {'Functions'};                                
% library_url{ith_library}     = 'https://github.com/ivsg-psu/FeatureExtraction_DataClean_BreakDataIntoLaps/archive/refs/tags/BreakDataIntoLaps_v2023_08_25.zip';
% 
% ith_library = ith_library+1;
% library_name{ith_library}    = 'ParseXODR_v2023_10_23';
% library_folders{ith_library} = {'Functions'};                                
% library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_MapTools_ParseXODR/archive/refs/tags/ParseXODR_v2023_10_23.zip';


%% Clear paths and folders, if needed
if 1==0
    clear flag_findEdge_Folders_Initialized;
    fcn_INTERNAL_clearUtilitiesFromPathAndFolders;
end

%% Do we need to set up the work space?
if ~exist('flag_findEdge_Folders_Initialized','var')
    this_project_folders = {'Functions','Data','LargeData'};
    fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders);
    flag_findEdge_Folders_Initialized = 1;
end

%% Set environment flags for input checking
% These are values to set if we want to check inputs or do debugging
% setenv('MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS','1');
% setenv('MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG','1');
setenv('MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS','1');
setenv('MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG','0');

%% Basic Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Functions
% ____            _        ______                _   _                 
% |  _ \          (_)      |  ____|              | | (_)                
% | |_) | __ _ ___ _  ___  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___ 
% |  _ < / _` / __| |/ __| |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
% | |_) | (_| \__ \ | (__  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
% |____/ \__,_|___/_|\___| |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% STEP 1: Load and study the data
% Set the "inputs" to the file loading process - need the date and names
% and date of file creation for the Vehicle Pose data file
fig_num = 1; 
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];

[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num));

% Plot the vehicle in 2D ENU
ENU_XY_fig_num = 2;
figure(ENU_XY_fig_num);
format = sprintf(' ''-'', ''Color'', [0 0 0], ''MarkerSize'', 10, ''LineWidth'', 3');
clf;
fcn_findEdge_plotVehicleXY(VehiclePose, format, fig_num); % -- add string as optional input (format)

Nscans = length(VehiclePose(:,1));
% Set defaults for which scans to extract
scanLineRange = [1400 1450];
ringsRange = []; % If leave empty, it loads all rings

% Extract scan lines
[VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), ([]));

% Plot the LIDAR in 2D ENU
figure(ENU_XY_fig_num);
format = sprintf(' ''-'', ''Color'', [0 0 1], ''MarkerSize'', 1');
fcn_findEdge_plotVehicleXY(LIDAR_ENU(:,1:2),format,ENU_XY_fig_num);


% Plot the vehicle in 3D ENU
ENU_XYZ_fig_num = 3;
figure(ENU_XYZ_fig_num);
clf;
scanLineRange = [1400 1450];
fcn_findEdge_plotVehicleXYZ(VehiclePose,(scanLineRange), (ENU_XYZ_fig_num)) 

% Plot the LIDAR in 3D ENU
scaling = [];
color_map = [];
fcn_findEdge_plotLIDARXYZ(LIDAR_ENU, (LIDAR_intensity), (scaling), (color_map), (ENU_XYZ_fig_num)); 

% Plot the LIDAR in XY, XZ, and YZ
LIDAR_XY_XZ_YZ_fig_num = 4;
color_triplet = [];
marker_size = [];
fcn_findEdge_plotLIDARXY(LIDAR_ENU, (color_triplet), (marker_size), (LIDAR_XY_XZ_YZ_fig_num))

% Plot the vehicle trajectory in LLA
LLA_fig_num = 5; % figure number
reference_LLA = [];
zoom_in_location = [];
zoomLevel = [];
LLA_VehiclePose = fcn_findEdge_plotVehicleLLA(VehiclePose, (reference_LLA), (zoom_in_location), (zoomLevel), (LLA_fig_num));

% plot the LIDAR in LLA 
fig_num = 6;
scaling = [];
color_map = 'autumn';
marker_size = [];
reference_LLA = [];
fcn_findEdge_plotLIDARLLA(LIDAR_ENU,(LIDAR_intensity),(scaling),(color_map),(marker_size),(reference_LLA),(fig_num))
%% CODE ABOVE THIS LINE WORKS, CODE BELOW THIS LINE NEEDS WORK 
%%%
% Jiabao and Alek - start here

%% Find the drivable surface --Jiabao
% This must be done in ENU coordinates because LLA is not an orthogonal
% coordinate system. To find the surface, we find the distance of the lidar
% points in the orthogonal direction by taking the dot product of the LIDAR
% points, relative to the vehicle center with the unit projection vector
% pointed to the left of the vehicle.
% Jiabao
gps_object = GPS();
% Calculate the vectors
vector_from_vehicle_pose_to_LIDAR_points = LIDAR_ENU - VehiclePose_ENU;

% Calculate the transverse distance
transverse_only_LIDAR_points = sum(vector_from_vehicle_pose_to_LIDAR_points.*VehiclePose_UnitOrthoVectors,2);

% Define lane width limits. Take 40 percent of the lane width. Numbers
% obtained by converting 12 ft (standard lane) to meters
lane_half_width = (3.6576/2) * 0.40;  

indicies_under_vehicle = find(abs(transverse_only_LIDAR_points)<lane_half_width);
LIDAR_ENU_under_vehicle = LIDAR_ENU(indicies_under_vehicle,:);
concatenate_LiDAR_LLA_points_under_vehicle = gps_object.ENU2WGSLLA(LIDAR_ENU_under_vehicle(:,1:3));


% Plot the LIDAR data underneath the vehicle in XYZ
ENU_3D_fig_num = 1;
figure(ENU_3D_fig_num);
plot3(LIDAR_ENU_under_vehicle(:,1),LIDAR_ENU_under_vehicle(:,2),LIDAR_ENU_under_vehicle(:,3), '.','Color',[0 1 0],'MarkerSize',1);
%%%%%%%%%%%% The plot above could be done by using fcn_findEdge_plotLIDARXYZ
% LIDAR_intensity = [];
% scaling = [];
% color_map = [];
% ENU_XYZ_fig_num = 2;
% fcn_findEdge_plotLIDARXYZ(LIDAR_ENU_under_vehicle, (LIDAR_intensity), (scaling), (color_map), (ENU_XYZ_fig_num));

% Plot the LIDAR data underneath the vehicle in LLA
figure(LLA_fig_num);
geoplot(concatenate_LiDAR_LLA_points_under_vehicle(:,1),concatenate_LiDAR_LLA_points_under_vehicle(:,2),'g.','MarkerSize',2);
geobasemap satellite
geotickformat -dd  % Sets the tick marks to decimal format, not degrees/minutes/seconds which is default
%%%%%%%%%% The plot above could be done by using fcn_findEdge_plotVehicleLLA
% LLA_fig_num = 5; % figure number
% reference_LLA = [];
% zoom_in_location = [];
% zoomLevel = [];
% LLA_VehiclePose = fcn_findEdge_plotVehicleLLA(LIDAR_ENU_under_vehicle, (reference_LLA), (zoom_in_location), (zoomLevel), (LLA_fig_num));

%% STEP 2: Find the driven path (left and right side points) (Yet to be functionalized)
% Alek

% This is done in "script_test_geometry_boundaryPointsDrivenPath" - Steven 

% This script can be found in "Functions" directory of "Geom Class" repo
% That script needs to be deleted once this library is done

%%% FORMERLY script_test_geometry_boundaryPointsDrivenPath
% This script is written to find the boundary points of driven path
%
% 2024_07_18 - Aneesh Batchu
% -- wrote the code originally

% To-DO: 
% Need to use "shift" (a varible name) to shift the boundary points

%%%% Calculate the vehicle orientation

vehicle_change_in_pose_XY = diff(VehiclePose(:,1:2));

% Repeat the last value again, since diff removes one row. We want the same
% number of vectors as the number of points, and diff removed one point.
vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];

% Convert these to unit vectors
unit_vehicle_change_in_pose_XY = fcn_INTERNAL_calcUnitVector(vehicle_change_in_pose_XY);

% Find orthogonal vetors by rotating by 90 degrees in the CCW direction
unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY*[0 1; -1 0];

% Define lane width limits. Take 40 percent of the lane width. Numbers
% obtained by converting 12 ft (standard lane) to meters
lane_half_width = (3.6576/2) * 0.40;  

% Transverse distance of the left and right boundary points of the driven
% path from vehicle center 
transverse_distance_of_boundary_points = [lane_half_width*unit_ortho_vehicle_vectors_XY, zeros(length(unit_ortho_vehicle_vectors_XY),1)];

shift = 5; 
% Shift
shift_distance = [unit_vehicle_change_in_pose_XY*shift, zeros(length(unit_vehicle_change_in_pose_XY),1)]; 

% Left boundary points of the driven path
left_boundary_points = VehiclePose(:,1:3) + transverse_distance_of_boundary_points - shift_distance; 

% right boundary points of the driven path
right_boundary_points = VehiclePose(:,1:3) - transverse_distance_of_boundary_points - shift_distance; 

% % Left boundary points of the driven path
% left_boundary_points = VehiclePose(:,1:3) + transverse_distance_of_boundary_points; 
% 
% % Left boundary points of the driven path
% right_boundary_points = VehiclePose(:,1:3) - transverse_distance_of_boundary_points; 


scanLineNumber_start = scanLineRange(1);
scanLineNumber_end   = scanLineRange(2);

boundaryLineNumber_start = max(scanLineNumber_start-1,1); 
boundaryLineNumber_end   = min(scanLineNumber_end+1, Nscans);

% lengths_boundary_points = [lane_half_width*unit_ortho_vehicle_vectors_XY(scanLineNumber_start:scanLineNumber_end,:), zeros((scanLineNumber_end-scanLineNumber_start+1),1)];
% 
% left_boundary_points = VehiclePose(scanLineNumber_start:scanLineNumber_end,1:3) + lengths_boundary_points;
% right_boundary_points = VehiclePose(scanLineNumber_start:scanLineNumber_end,1:3) - lengths_boundary_points;

fig_num = 1098;
figure(fig_num);clf;

hold on;
grid on;
axis equal

plot3(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1),VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,2),VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 0 0],'MarkerSize',30,'LineWidth',3);
plot3(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1),VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,2),VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[1 1 0],'MarkerSize',10,'LineWidth',3);

% Show the orthogonal arrows showing vehicle motion directions. Green
% is forward, bLue is Left
quiver3(...
    VehiclePose(boundaryLineNumber_start,1),VehiclePose(boundaryLineNumber_start,2),VehiclePose(boundaryLineNumber_start,3), ...
    unit_vehicle_change_in_pose_XY(boundaryLineNumber_start,1),unit_vehicle_change_in_pose_XY(boundaryLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 0 1]);
quiver3(...
    VehiclePose(boundaryLineNumber_start,1),VehiclePose(boundaryLineNumber_start,2),VehiclePose(boundaryLineNumber_start,3), ...
    unit_ortho_vehicle_vectors_XY(boundaryLineNumber_start,1),unit_ortho_vehicle_vectors_XY(boundaryLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 0 1]);

plot3(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1),left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,2),left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 1 0],'MarkerSize',30,'LineWidth',3);
plot3(right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1),right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,2),right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 0 1],'MarkerSize',30,'LineWidth',3);


xlabel('East position [m]');
ylabel('North position [m]');
zlabel('Up position [m]');
view(3)
hold off

%% -- Alek

ENU_3D_fig_num = 3;
figure(ENU_3D_fig_num);

plot3(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1),left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,2),left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 1 0],'MarkerSize',30,'LineWidth',3);
plot3(right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1),right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,2),right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 0 1],'MarkerSize',30,'LineWidth',3);

boundary_points_driven_path = [right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3);
    flipud(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3));
    right_boundary_points(boundaryLineNumber_start,1:3)];

boundary_points_driven_path_LLA = gps_object.ENU2WGSLLA(boundary_points_driven_path);

LLA_fig_num = 2;
figure(LLA_fig_num);
% Plot the LIDAR data underneath the vehicle in LLA
figure(LLA_fig_num);
hold off
geoplot(boundary_points_driven_path_LLA(:,1),boundary_points_driven_path_LLA(:,2),'b.','MarkerSize',30);


% %%  INPOLYGON
% boundary_points = [right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3); flipud(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3)); right_boundary_points(boundaryLineNumber_start,1:3)];
% 
% points = VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1:3);
% 
% % xv = [0.5;0.2;1.0;0;0.8;0.5];
% % yv = [1.0;0.1;0.7;0.7;0.1;1];
% % 
% % 
% % xq = [0.1;0.5;0.9;0.2;0.4;0.5;0.5;0.9;0.6;0.8;0.7;0.2];
% % yq = [0.4;0.6;0.9;0.7;0.3;0.8;0.2;0.4;0.4;0.6;0.2;0.6];
% 
% 
% [in,on] = inpolygon(points(:,1),points(:,2),boundary_points(:,1),boundary_points(:,2));
% 
% figure(122)
% 
% plot(boundary_points(:,1),boundary_points(:,2)) % polygon
% 
% hold on
% plot(points(in,1),points(in,2),'r+') % points strictly inside
% % plot(xq(on),yq(on),'k*') % points on edge
% plot(points(~in,1),points(~in,2),'bo') % points outside
% hold off

%% STEP 3 & STEP 4: Seperate the data into grids, and classify the grids as the grids with zero points and grids with more than zero points
% Jiabao - convert this out of geometry library

% These are concatenated LiDAR points of chosen scans and cells in the
% first step. 
% LiDAR_allPoints = concatenate_LiDAR_XYZ_points(:,1:3);

LiDAR_allPoints = [LIDAR_ENU, LIDAR_scanLineAndRingID];

% remove NANs
LiDAR_allPoints = LiDAR_allPoints(~isnan(LiDAR_allPoints(:,1)),:);

% scanLine_rings without NaNs
LIDAR_scanLines = LiDAR_allPoints(:,4:5); 


% Input points for seperating the data into grids. The points are in 2D as
% the analysis is carried out in 2D
input_points = LiDAR_allPoints(:,1:2); 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1; %0.8;%1;%1.25; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_findEdge_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries in #D
% grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

% The 2D grid boundaries required for the current analysis
grid_boundaries = [Min_x Max_x Min_y Max_y]; 

% plot all LiDAR points                     ----------- To-Do: ALEKS
%
% Change the function name into fcn_geometry_plotPointsinLLA
% - Add an option for chnaging the marker such as '.' (currently plots),
% '+', '*', 'o' etc
% BUG: This function does not work if Legend options are given as zero.
% Write some test cases in the script
% BUG: This function should work even if the legend_option is empty
% BUG: This function should work even if the Legend_name is empty
% Check this function carefully. Write all the test script clearly. It
% should cover all the cases. 
% 
% marker_size = 10;
% RGB_triplet = [0.8, 0.8, 0.8]; 
% legend_option = 1;
% legend_name = 'LiDAR Points';
% legend_position = [];
% plot_LiDAR_allPoints = LiDAR_allPoints; 
% [~] = fcn_geometry_plotGridCenters(LiDAR_allPoints,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);

% Define GPS object - This should be the input for the above function 
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;
% Define GPS object
gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class

% Use the class to convert ENU to LLA
LIDAR_allPoints_LLA = gps_object.ENU2WGSLLA(LiDAR_allPoints(:,1:3));

% Currently plotting here without uisng any function 
fig_num_LLA = 30;
% Plot the ENU results
figure(fig_num_LLA);clf;

geoplot(LIDAR_allPoints_LLA(:,1),LIDAR_allPoints_LLA(:,2),'mo','MarkerSize',10);
hold on
geoplot(LIDAR_allPoints_LLA(:,1),LIDAR_allPoints_LLA(:,2),'k.','MarkerSize',10);
geoplot(boundary_points_driven_path_LLA(:,1),boundary_points_driven_path_LLA(:,2),'g.','MarkerSize',10);
title('LLA Trace geometry')

geobasemap satellite
geotickformat -dd  % Sets the tick marks to decimal format, not degrees/minutes/seconds which is default
%%%%%%%%%%--the plot above could be replace with this fcn_findEdge_plotVehicleXY
% fig_num = 1514874;
% format = sprintf('''mo'',''MarkerSize'',10');
% fcn_findEdge_plotVehicleXY(LiDAR_allPoints(:,1:2),format,fig_num);
% hold on 
% format = sprintf('''k.'',''MarkerSize'',10');
% fcn_findEdge_plotVehicleXY(LiDAR_allPoints(:,1:2),format,fig_num);
% format = sprintf('''g.'',''MarkerSize'',10');
% fcn_findEdge_plotVehicleXY(boundary_points_driven_path(:,1:2),format,fig_num);
% hold off

% however, this will not do a good job sicne we need to edit the color and
% shape of the points. In addition, you can not overlap this. 

fig_num_first_classification = 40; 
figure(fig_num_first_classification);clf

[gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, grids_with_zero_points,...
    grids_greater_than_zero_points, gridCenters_zero_point_density,...
    gridCenters_greater_than_zero_point_density, gridIndices, grid_AABBs] = fcn_findEdge_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,fig_num_first_classification);

%% STEP 5: Find the driven path grids within the grids more than zero points
% Jiabao
% NEED TO BE FUNCTION


% -----------------------------NOTE------------------------------
% After finding the grids without anypoints, the grids are completely
% removed from the analysis. Only, grids with greater than zero points were
% analyzed from here. 
% -----------------------------NOTE------------------------------

% Plot all the grids greater than zero point density

% Figure number
fig_num_gridLines_greater_than_zero_point_density = 50;
figure(fig_num_gridLines_greater_than_zero_point_density); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Gridlines and points of the grids greater than zero point density')

% Pre-allocation: To find the total number of points in each grid
total_points_in_grids_greater_than_zero = zeros(length(grids_greater_than_zero_points),1);

% allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied 
gridlines_grids_greater_than_zero = zeros(11*length(grids_greater_than_zero_points),2); % length(gridlines) = 11

concatenate_gridPoints_scanLines_rings = [];

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(grids_greater_than_zero_points)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.2 0.2 0.2];

    % Plot current AABB
    current_AABB = grid_AABBs(grids_greater_than_zero_points(ith_grid),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,1) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,2) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,3); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,4); ...
        current_AABB(1,2) current_AABB(1,4); ...
        ];

    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==grids_greater_than_zero_points(ith_grid);
    
    % XY coordinates of the input points that are in ith_grid
    points_in_domain = input_points(rows_in_domain,:);
 
    % Scan lines and rings
    scanLines_and_rings = LIDAR_scanLines(rows_in_domain,:);

    % Find number of LiDAR scan lines in each grid
    scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));
    % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP
    % total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid]; % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP

    % Combine points_in_domain and ScanLines and rings
    gridPoints_scanLines_rings_to_add = [ith_grid*(ones(length(points_in_domain(:,1)),1)),points_in_domain,scanLines_and_rings];
    % gridPoints_scanLines_rings_to_add = [gridPoints_scanLines_rings_to_add; naa nan nan nan]; %#ok<AGROW>

    % Concatenate the points in domain, scan lines and rings
    concatenate_gridPoints_scanLines_rings = [concatenate_gridPoints_scanLines_rings; gridPoints_scanLines_rings_to_add]; %#ok<AGROW>


    % Find the total number of points in each grid
    total_points_in_grids_greater_than_zero(ith_grid) = length(points_in_domain);
    
    % Save the grid lines of all the grids greater than zero density in a
    % matrix
    length_gridlines = length(gridlines);
    gridlines_grids_greater_than_zero(1+(ith_grid-1)*length_gridlines:ith_grid*length_gridlines,:) = gridlines;
   
    % Plot the result
    % plot the grid lines of the grids greater than zero
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    % plot the points in the grid
    plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)

end

% Write the grid number at the grid center for reference. 

for ith_text = 1:length(grids_greater_than_zero_points(:,1))
    current_text = sprintf('%.0d',ith_text);

    % Place the text on the grid center
    text(gridCenters_greater_than_zero_point_density(ith_text,1), gridCenters_greater_than_zero_point_density(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end
%%%%%%%%%%--the plot above could be replace with this fcn_findEdge_plotVehicleXY
% fig_num = 1514874;
% fcn_findEdge_plotVehicleXY(LiDAR_allPoints(:,1:2),fig_num);
% however, this will not do a good job sicne we need to edit the color and
% shape of the points. In addition, you can not overlap this. 

%% Testing: Delete this afterwards --Alek


% Figure number
fig_num_gridLines_greater_than_zero_point_density = 545;
figure(fig_num_gridLines_greater_than_zero_point_density); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Gridlines and points of the grids greater than zero point density')


% allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied 
gridlines_grids_greater_than_zero = zeros(11*length(grids_greater_than_zero_points),2); % length(gridlines) = 11

concatenate_gridPoints_scanLines_rings = [];
orthogonal_dist_each_grid = []; 

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(grids_greater_than_zero_points)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.2 0.2 0.2];

    % Plot current AABB
    current_AABB = grid_AABBs(grids_greater_than_zero_points(ith_grid),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,1) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,2) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,3); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,4); ...
        current_AABB(1,2) current_AABB(1,4); ...
        ];

    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==grids_greater_than_zero_points(ith_grid);
    
    % XY coordinates of the input points that are in ith_grid
    points_in_domain = input_points(rows_in_domain,:);
 
    % Scan lines and rings
    scanLines_and_rings = LIDAR_scanLines(rows_in_domain,:);

    % Find number of LiDAR scan lines in each grid
    scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));
    % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP
    % total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid]; % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP

    % Combine points_in_domain and ScanLines and rings
    gridPoints_scanLines_rings_to_add = [ith_grid*(ones(length(points_in_domain(:,1)),1)),points_in_domain,scanLines_and_rings];

    % Sort gridPoints_scanLines_rings_to_add based on scan line
    [~,sorted_gridPoints_scanLines_rings] = sort(gridPoints_scanLines_rings_to_add(:,4));

    % Sorted gridPoints_scanLines_rings_to_add matrix
    gridPoints_scanLines_rings_to_add = gridPoints_scanLines_rings_to_add(sorted_gridPoints_scanLines_rings,:);

    % Indices first scan line of the matrix as a seperate matrix
    indices_gridPoints_scanLines_first_scan = find(gridPoints_scanLines_rings_to_add(:,4) == gridPoints_scanLines_rings_to_add(1,4)); 

    % Seperate the first scan line of the matrix as a seperate matrix
    gridPoints_scanLines_first_scan = gridPoints_scanLines_rings_to_add(indices_gridPoints_scanLines_first_scan,:);

    %
    if length(gridPoints_scanLines_first_scan(:,1)) > 1 && (gridPoints_scanLines_first_scan(1,5) == gridPoints_scanLines_first_scan(2,5))
        change_in_vector = gridPoints_scanLines_first_scan(2,2:3) - gridPoints_scanLines_first_scan(1,2:3);
        unit_change_in_vector = fcn_INTERNAL_calcUnitVector(change_in_vector);
        orth_unit_change_in_vector = unit_change_in_vector*[0 1; -1 0];
        
        % The remaining number of grids
        remaining_grids = length(indices_gridPoints_scanLines_first_scan)+1:length(gridPoints_scanLines_rings_to_add(:,1)); 
        
        % 
        vector_from_base_point_first_scan_to_points_in_otherScans_rings = gridPoints_scanLines_rings_to_add(remaining_grids,2:3) - ...
            gridPoints_scanLines_first_scan(1,2:3).*(ones(length(remaining_grids),2));

        % Unit orthogonal vector
        repeated_orth_unit_change_in_vector = orth_unit_change_in_vector.*(ones(length(remaining_grids),2));

        % Calculate the transverse distance
        transverse_dist_grid_points_other_scanLines = sum(vector_from_base_point_first_scan_to_points_in_otherScans_rings.*repeated_orth_unit_change_in_vector,2);

        % Mean of absolute values of transverse distances
        mean_dist = mean(abs(transverse_dist_grid_points_other_scanLines));

        orthogonal_dist_each_grid = [orthogonal_dist_each_grid; mean_dist]; %#ok<AGROW>
    else 

        orthogonal_dist_each_grid = [orthogonal_dist_each_grid; nan]; %#ok<AGROW>

    end

   

    % Concatenate the points in domain, scan lines and rings
    concatenate_gridPoints_scanLines_rings = [concatenate_gridPoints_scanLines_rings; gridPoints_scanLines_rings_to_add]; %#ok<AGROW>

    % Save the grid lines of all the grids greater than zero density in a
    % matrix
    length_gridlines = length(gridlines);
    gridlines_grids_greater_than_zero(1+(ith_grid-1)*length_gridlines:ith_grid*length_gridlines,:) = gridlines;
   
    % Plot the result
    % plot the grid lines of the grids greater than zero
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    % plot the points in the grid
    plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)

end

% Write the grid number at the grid center for reference. 

for ith_text = 1:length(grids_greater_than_zero_points(:,1))
    % current_text = sprintf('%.0d',current_mapped_grids(ith_text));
    current_text = sprintf('%.0d',(ith_text));

    % Place the text on the grid center
    text(gridCenters_greater_than_zero_point_density(ith_text,1), gridCenters_greater_than_zero_point_density(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

fig_num = 5837; 
figure(fig_num); clf;
hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(grids_greater_than_zero_points)

    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.2 0.2 0.2];

    % Plot current AABB
    current_AABB = grid_AABBs(grids_greater_than_zero_points(ith_grid),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,1) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,2) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,3); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,4); ...
        current_AABB(1,2) current_AABB(1,4); ...
        ];

    % Plot the result
    % plot the grid lines of the grids greater than zero
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);

end

% Write the grid number at the grid center for reference. 

for ith_text = 1:length(grids_greater_than_zero_points(:,1))
    current_text = sprintf('%.3f',orthogonal_dist_each_grid(ith_text));
    % current_text = sprintf('%.3f',orthogonal_dist_each_grid);

    % Place the text on the grid center
    text(gridCenters_greater_than_zero_point_density(ith_text,1), gridCenters_greater_than_zero_point_density(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end


%% STEP 5: Plot grid centers and boundary points in a same figure in ENU --Alek

% "inpolygon" is used to find the grids within the boundary points 
[in,on] = inpolygon(gridCenters_greater_than_zero_point_density(:,1),gridCenters_greater_than_zero_point_density(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));

% Original grid numbers of driven path
original_grid_numbers_of_driven_path = grids_greater_than_zero_points(in); 

% Current grid numbers in driven path 
current_grid_numbers_of_driven_path = find(in); 

% Total points in each grid in the driven path
total_points_in_each_grid_in_the_driven_path = total_N_points_in_each_grid(original_grid_numbers_of_driven_path); 

% Total points in each grid with points greater than zero
total_points_in_each_grid_with_points_greater_than_zero = total_N_points_in_each_grid(grids_greater_than_zero_points); 

% Grid centers of the driven path
gridCenters_driven_path = [gridCenters_greater_than_zero_point_density(in,1),gridCenters_greater_than_zero_point_density(in,2)];

grid_corners = grid_AABBs(grids_greater_than_zero_points(in),1:4); 

% Extract the coordinates
x_low = grid_corners(:,1); 
x_high = grid_corners(:,2);
y_low = grid_corners(:,3);
y_high =  grid_corners(:,4); 

% Rearrange grid corners
rearranged_grid_corners = [x_low y_low; x_high y_low; x_low y_high; x_high y_high]; 


% "inpolygon" is used to find the grid corners within the boundary points 
[in_bb,on_bb] = inpolygon(rearranged_grid_corners(:,1),rearranged_grid_corners(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));


inliers_corners = rearranged_grid_corners(in_bb,:); 


fig_num_gridCenters_and_boundary_points_greater_than_zero = 51;
figure(fig_num_gridCenters_and_boundary_points_greater_than_zero); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers and boundary points')

plot(gridCenters_greater_than_zero_point_density(:,1), gridCenters_greater_than_zero_point_density(:,2), '.','MarkerSize',30,'Color',[0.8 0.8 0.8]);
plot(boundary_points_driven_path(:,1), boundary_points_driven_path(:,2), '.', 'MarkerSize',30, 'Color',[0 1 0]); 
plot(rearranged_grid_corners(in_bb,1), rearranged_grid_corners(in_bb,2), '.', 'MarkerSize',10, 'Color',[0 1 0]); 

% plot the grids in the driven path
plot(gridCenters_greater_than_zero_point_density(in,1),gridCenters_greater_than_zero_point_density(in,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',1.5) % points strictly inside
aaa = grid_AABBs(grids_greater_than_zero_points(in),1:4);

for ith_text = 1:length(grids_greater_than_zero_points(:,1))
    current_text = sprintf('%.0d',ith_text);
    % Place the text on the grid center
    text(gridCenters_greater_than_zero_point_density(ith_text,1), gridCenters_greater_than_zero_point_density(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 6);
end


%% STEP 6: Statistic 1

% Figure number of histogram
fig_histogram = 52; 
figure(fig_histogram); clf; 
% Add labels and title 
hold on
grid on
xlabel('Points per grid');
ylabel('Frequency');
title('Statistic 1: Determining suitable point density');


% edges = (floor(min(total_points_in_each_grid_with_points_greater_than_zero)/10)*10):10:(ceil(max(total_points_in_each_grid_with_points_greater_than_zero)/10)*10); % Define the bin edges

% Create the histogram
% actual_driving_surface_grids_hist = histogram(total_points_in_each_grid_in_the_driven_path,edges,'Visible','on'); 

% total_grids_hist = histogram(total_points_in_each_grid_with_points_greater_than_zero,edges,'Visible','on'); 
total_grids_greater_than_zero_hist = histogram(total_points_in_each_grid_with_points_greater_than_zero,20,'Visible','on'); 
actual_driven_path_grids_hist = histogram(total_points_in_each_grid_in_the_driven_path,10,'Visible','on'); 

% Extract the counts for both histograms
counts1 = total_grids_greater_than_zero_hist.Values;
[~,index_max_counts2] = max(counts1);
counts2 = actual_driven_path_grids_hist.Values;

binEdges = total_grids_greater_than_zero_hist.BinEdges;

% Calculate the overlapping counts
% % overlapCounts = min(counts2, counts1);

% Find a ratio
 % point_density = sum(binEdges(index_max_counts2:(index_max_counts2+1)))/2; 

% point_density = floor(mean(total_points_in_each_grid_in_the_driven_path) - 7.5*(std(total_points_in_each_grid_in_the_driven_path)));
point_density = floor(mean(total_points_in_each_grid_in_the_driven_path) - 1.5*(std(total_points_in_each_grid_in_the_driven_path)));

% point_density = floor(mean(total_points_in_each_grid_in_the_driven_path));


% Minimum number of points required 
point_density = floor(20*((grid_size^2)/(0.3^2))); 
disp('Chosen point density')
disp(point_density)
% mean(total_points_in_each_grid_in_the_driven_path)/mean(total_points_in_each_grid_with_points_greater_than_zero);

plot(point_density,0, 'k.', 'MarkerSize',20)
current_text = sprintf('point density = %.2d',point_density);
text(650, 200,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');


%% STEP 6: Statistic 2 - Determine number of LiDAR scans in each grid

% This was also done in STEP 5, however, it was done using a for loop. Need to
% do it without using a for loop.

total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(grids_greater_than_zero_points)
     %Get all points in this domain and plot them
    rows_in_domain = gridIndices==grids_greater_than_zero_points(ith_grid);

     %Find number of LiDAR scan lines in each grid
    scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));
  
    % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORK     S WITHOUT A FOR LOOP
    total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid]; % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP
end


% total_scan_lines_in_each_grid
%wrote by Jiabao Zhao

%total_scan_lines_in_each_grid_with_more_than_zero_points = [];
%rows_in_domain = ismember(gridIndices, grids_greater_than_zero_points);
%scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));
%total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid];


%% STEP 7: Classify the grids with more than zero points into mapped and unmapped grids
% Jiabao
% plotting
fig_num = 72; 
figure(fig_num);clf

% Classify the grids with more than zero points into mapped and unmapped grids
[original_grids_with_low_point_density, original_grids_with_required_point_density, original_grids_with_more_than_one_scan_line, original_grids_with_one_scan_line, ...
    original_mapped_grids, original_unmapped_grids, gridCenters_low_point_density, gridCenters_required_point_density, gridCenters_with_more_than_one_scan_line, ...
gridCenters_with_one_scan_line, gridCenters_mapped_grids, gridCenters_unmapped_grids, current_grids_with_low_point_density, current_grids_with_required_point_density, ...
current_grids_with_more_than_one_scan_line, current_grids_with_one_scan_line, current_mapped_grids, current_unmapped_grids]...
= fcn_findEdge_GridsIntoMappedUnmapped(point_density, total_N_points_in_each_grid, total_scan_lines_in_each_grid_with_more_than_zero_points, ...
grids_greater_than_zero_points, gridCenters, fig_num); 

% Plot the grids with low point density and required density 
fig_num = 70; 
figure(fig_num);clf

marker_size = 30;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Grids with low point density';
legend_position = [];
marker_type = [];
plot_gridCenters_low_point_density = [gridCenters_low_point_density, zeros(length(gridCenters_low_point_density(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_low_point_density,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Grids with required density';
legend_position = [];
marker_type = [];
plot_gridCenters_required_point_density = [gridCenters_required_point_density, zeros(length(gridCenters_required_point_density(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_required_point_density,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% Plot grids with one scan line and more than one scan line
fig_num = 71; 
figure(fig_num);clf

marker_size = 30;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Grids with one scan line';
legend_position = [];
marker_type = [];
plot_gridCenters_with_one_scan_line = [gridCenters_with_one_scan_line, zeros(length(gridCenters_with_one_scan_line(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_with_one_scan_line,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot grid centers
marker_size = 30;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Grids with more than one scan line';
legend_position = [];
marker_type = [];
plot_gridCenters_with_more_than_one_scan_line = [gridCenters_with_more_than_one_scan_line, zeros(length(gridCenters_with_more_than_one_scan_line(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_with_more_than_one_scan_line,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% ENU - plotting mapped and unmapped
fig_num_ENU = 73; 
figure(fig_num_ENU)
clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers of mapped and unmapped grids and boundary points')

plot(gridCenters_unmapped_grids(:,1), gridCenters_unmapped_grids(:,2), '.','MarkerSize',40,'Color',[0.8 0.8 0.8]);
plot(gridCenters_mapped_grids(:,1), gridCenters_mapped_grids(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);

for ith_text = 1:length(current_unmapped_grids(:,1))
    current_text = sprintf('%.0d',current_unmapped_grids(ith_text));
    % Place the text on the grid center
    text(gridCenters_unmapped_grids(ith_text,1), gridCenters_unmapped_grids(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 8);
end

for ith_text = 1:length(current_mapped_grids(:,1))
    current_text = sprintf('%.0d',current_mapped_grids(ith_text));
     % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[1 1 0],'HorizontalAlignment','center','FontSize', 8);
end


%% STEP 7.25: Find the orthogonal distances (NEW STATISTIC)

% Figure number
fig_num_gridLines_greater_than_zero_point_density = 50;
figure(fig_num_gridLines_greater_than_zero_point_density); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Gridlines and points of the grids greater than zero point density')


% allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied 
gridlines_grids_greater_than_zero = zeros(11*length(original_mapped_grids),2); % length(gridlines) = 11

concatenate_gridPoints_scanLines_rings = [];
orthogonal_dist_each_grid = []; 
transverse_span_each_grid = [];

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(original_mapped_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.2 0.2 0.2];

    % Plot current AABB
    current_AABB = grid_AABBs(original_mapped_grids(ith_grid),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,1) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,2) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,3); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,4); ...
        current_AABB(1,2) current_AABB(1,4); ...
        ];

    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==original_mapped_grids(ith_grid);
    
    % XY coordinates of the input points that are in ith_grid
    points_in_domain = input_points(rows_in_domain,:);
 
    % Scan lines and rings
    scanLines_and_rings = LIDAR_scanLines(rows_in_domain,:);

    % Find number of LiDAR scan lines in each grid
    scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));
    % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP
    % total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid]; % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP

    % Combine points_in_domain and ScanLines and rings
    gridPoints_scanLines_rings_to_add = [ith_grid*(ones(length(points_in_domain(:,1)),1)),points_in_domain,scanLines_and_rings];

    % Sort gridPoints_scanLines_rings_to_add based on scan line
    [~,sorted_gridPoints_scanLines_rings] = sort(gridPoints_scanLines_rings_to_add(:,4));

    % Sorted gridPoints_scanLines_rings_to_add matrix
    gridPoints_scanLines_rings_to_add = gridPoints_scanLines_rings_to_add(sorted_gridPoints_scanLines_rings,:);

    % Count occurrences of each unique number in scan lines
    [uniqueNumbers, ~, uniqueNum_idx] = unique(gridPoints_scanLines_rings_to_add(:,4));
    counts = histc(uniqueNum_idx, 1:numel(uniqueNumbers));

    % Create the new array with the counts
    count_array = counts;

    % Index of the scan line with more than one occurence
    index_of_scanLines = find(count_array>1, 1);

    % Indices first scan line of the matrix as a seperate matrix
    indices_gridPoints_scanLines_first_scan = find(gridPoints_scanLines_rings_to_add(:,4) == gridPoints_scanLines_rings_to_add(index_of_scanLines,4));

    % Seperate the scan line with more than one occurence of the matrix as a seperate matrix
    gridPoints_scanLines_first_scan = gridPoints_scanLines_rings_to_add(indices_gridPoints_scanLines_first_scan,:);

    % Count occurrences of each unique number in rings
    [uniqueNumbers, ~, uniqueNum_idx] = unique(gridPoints_scanLines_first_scan(:,5));
    counts = histc(uniqueNum_idx, 1:numel(uniqueNumbers));

    % Create the new array with the counts
    count_array = counts;

    % Index of the scan line with more than one occurence
    index_of_rings = find(count_array>1, 1);

    % if length(gridPoints_scanLines_first_scan(:,1)) == 1
    % 
    %     indices_gridPoints_scanLines_first_scan

    if length(gridPoints_scanLines_first_scan(:,1)) > 1 & (gridPoints_scanLines_first_scan(index_of_rings,5) == gridPoints_scanLines_first_scan(index_of_rings+1,5))
        change_in_vector = gridPoints_scanLines_first_scan(2,2:3) - gridPoints_scanLines_first_scan(1,2:3);
        unit_change_in_vector = fcn_INTERNAL_calcUnitVector(change_in_vector);
        orth_unit_change_in_vector = unit_change_in_vector*[0 1; -1 0];
        
        % The remaining number of grids
        remaining_grids = length(indices_gridPoints_scanLines_first_scan)+1:length(gridPoints_scanLines_rings_to_add(:,1)); 
        
        % 
        vector_from_base_point_first_scan_to_points_in_otherScans_rings = gridPoints_scanLines_rings_to_add(remaining_grids,2:3) - ...
            gridPoints_scanLines_first_scan(1,2:3).*(ones(length(remaining_grids),2));

        % Unit orthogonal vector
        repeated_orth_unit_change_in_vector = orth_unit_change_in_vector.*(ones(length(remaining_grids),2));

        % Calculate the transverse distance
        transverse_dist_grid_points_other_scanLines = sum(vector_from_base_point_first_scan_to_points_in_otherScans_rings.*repeated_orth_unit_change_in_vector,2);

        % Positive transverse distances
        positive_transverse_dist_grid_points_other_scanLines = transverse_dist_grid_points_other_scanLines(transverse_dist_grid_points_other_scanLines>=0);

        % Negative transverse distances
        negative_transverse_dist_grid_points_other_scanLines = transverse_dist_grid_points_other_scanLines(transverse_dist_grid_points_other_scanLines<0);

        % maximum span distance
        if ~isempty(positive_transverse_dist_grid_points_other_scanLines) && ~isempty(negative_transverse_dist_grid_points_other_scanLines)

            maximum_span_distance = max(positive_transverse_dist_grid_points_other_scanLines) + max(abs(negative_transverse_dist_grid_points_other_scanLines));

        elseif isempty(negative_transverse_dist_grid_points_other_scanLines)

            maximum_span_distance = max(positive_transverse_dist_grid_points_other_scanLines);

        elseif isempty(positive_transverse_dist_grid_points_other_scanLines)

            maximum_span_distance = max(abs(negative_transverse_dist_grid_points_other_scanLines));

        end
        % Conatenate maximum transverse span
        transverse_span_each_grid = [transverse_span_each_grid; maximum_span_distance]; %#ok<AGROW>

        % Mean of absolute values of transverse distances
        mean_dist = mean(abs(transverse_dist_grid_points_other_scanLines));

        % Concatenate the orthogonal distances
        orthogonal_dist_each_grid = [orthogonal_dist_each_grid; mean_dist]; %#ok<AGROW>
    else

        % Conacatenate orthogonal distance
        orthogonal_dist_each_grid = [orthogonal_dist_each_grid; nan]; %#ok<AGROW>

        % Conatenate maximum transverse span
        transverse_span_each_grid = [transverse_span_each_grid; nan]; %#ok<AGROW>

    end
    % Concatenate the points in domain, scan lines and rings
    concatenate_gridPoints_scanLines_rings = [concatenate_gridPoints_scanLines_rings; gridPoints_scanLines_rings_to_add]; %#ok<AGROW>

    % Save the grid lines of all the grids greater than zero density in a
    % matrix
    length_gridlines = length(gridlines);
    gridlines_grids_greater_than_zero(1+(ith_grid-1)*length_gridlines:ith_grid*length_gridlines,:) = gridlines;
   
    % Plot the result
    % plot the grid lines of the grids greater than zero
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    % plot the points in the grid
    plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)

end

% Write the grid number at the grid center for reference. 

for ith_text = 1:length(original_mapped_grids(:,1))
    % current_text = sprintf('%.0d',current_mapped_grids(ith_text));
    current_text = sprintf('%.0d',(ith_text));

    % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

fig_num = 5837; 
figure(fig_num); clf;
hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(original_mapped_grids)

    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.2 0.2 0.2];

    % Plot current AABB
    current_AABB = grid_AABBs(original_mapped_grids(ith_grid),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,1) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,2) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,3); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,4); ...
        current_AABB(1,2) current_AABB(1,4); ...
        ];

    % Plot the result
    % plot the grid lines of the grids greater than zero
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);

end

% Write the grid number at the grid center for reference. 

for ith_text = 1:length(current_mapped_grids(:,1))
    current_text = sprintf('%.3f',orthogonal_dist_each_grid(ith_text));
    % current_text = sprintf('%.3f',orthogonal_dist_each_grid);

    % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

fig_num = 58309; 
figure(fig_num); clf;
hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Transverse span distances')

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(original_mapped_grids)

    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.2 0.2 0.2];

    % Plot current AABB
    current_AABB = grid_AABBs(original_mapped_grids(ith_grid),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,1) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,2) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,3); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,4); ...
        current_AABB(1,2) current_AABB(1,4); ...
        ];

    % Plot the result
    % plot the grid lines of the grids greater than zero
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);

end

% Write the grid number at the grid center for reference. 

for ith_text = 1:length(current_mapped_grids(:,1))
    current_text = sprintf('%.3f',transverse_span_each_grid(ith_text));
    % current_text = sprintf('%.3f',orthogonal_dist_each_grid);

    % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end


%% Updated mapped grids
% Jiabao
% Threshold of transverse span
threshold_transverse_dist = 0.15;

% Updated original mapped grids
updated_original_mapped_grids = original_mapped_grids(orthogonal_dist_each_grid>threshold_transverse_dist); 

% Updated original unmapped grids
updated_original_unmapped_grids = sort([original_unmapped_grids; original_mapped_grids(orthogonal_dist_each_grid<=threshold_transverse_dist)]); 

% Updated current mapped grids 
updated_current_mapped_grids = current_mapped_grids(orthogonal_dist_each_grid>threshold_transverse_dist);

% Updated current unmapped grids 
updated_current_unmapped_grids = sort([current_unmapped_grids; current_mapped_grids(orthogonal_dist_each_grid<=threshold_transverse_dist)]);

% Updated original mapped grid centers
gridCenters_updated_original_mapped_grids = gridCenters(updated_original_mapped_grids,1:2); 

% Updated original unmapped grid centers
gridCenters_updated_original_unmapped_grids = gridCenters(updated_original_unmapped_grids,1:2); 

% Plot the grids with low point density and required density 
fig_num = 8765; 
figure(fig_num);clf

marker_size = 30;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Updated mapped grids';
legend_position = [];
marker_type = [];
plot_gridCenters_updated_original_mapped_grids = [gridCenters_updated_original_mapped_grids, zeros(length(gridCenters_updated_original_mapped_grids(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_updated_original_mapped_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot grid centers
marker_size = 30;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Updated unmapped grids';
legend_position = [];
marker_type = [];
plot_gridCenters_updated_original_unmapped_grids = [gridCenters_updated_original_unmapped_grids, zeros(length(gridCenters_updated_original_unmapped_grids(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_updated_original_unmapped_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


%% Find mapped grids current number

% ENU - plotting mapped and unmapped
fig_num_ENU = 74; 
figure(fig_num_ENU)
clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers of mapped grids')

plot(gridCenters_mapped_grids(:,1), gridCenters_mapped_grids(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',2) % points strictly inside



for ith_text = 1:length(current_mapped_grids(:,1))
    current_text = sprintf('%.0d',(ith_text));
     % Place the text on the grid center
    text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[1 1 1],'HorizontalAlignment','center','FontSize', 6, 'FontWeight','bold');
end

%% STEP 7.5: Find driven path of the vehicle in mapped grids

% "inpolygon" is used to find the grids within the boundary points 
[in,on] = inpolygon(gridCenters_updated_original_mapped_grids(:,1),gridCenters_updated_original_mapped_grids(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));

% Original grid numbers of driven path
original_grid_numbers_of_driven_path = updated_original_mapped_grids(in); 

% Current grid numbers in driven path 
current_grid_numbers_of_driven_path = updated_current_mapped_grids(in); %find(in); 

% Total points in each grid in the driven path
total_points_in_each_grid_in_the_driven_path = total_N_points_in_each_grid(original_grid_numbers_of_driven_path); 

% Total points in each grid with points greater than zero
total_points_in_each_grid_with_points_greater_than_zero = total_N_points_in_each_grid(updated_current_mapped_grids); 

% Grid centers of the driven path
gridCenters_driven_path = [gridCenters_updated_original_mapped_grids(in,1),gridCenters_updated_original_mapped_grids(in,2)];


fig_num_gridCenters_and_boundary_points_greater_than_zero = 51;
figure(fig_num_gridCenters_and_boundary_points_greater_than_zero); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers and boundary points')

plot(gridCenters_updated_original_mapped_grids(:,1), gridCenters_updated_original_mapped_grids(:,2), '.','MarkerSize',30,'Color',[0.8 0.8 0.8]);
plot(boundary_points_driven_path(:,1), boundary_points_driven_path(:,2), '.', 'MarkerSize',30, 'Color',[0 1 0]); 

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',1.5) % points strictly inside

for ith_text = 1:length(gridCenters_updated_original_mapped_grids(:,1))
    current_text = sprintf('%.0d',ith_text);
    % Place the text on the grid center
    text(gridCenters_updated_original_mapped_grids(ith_text,1), gridCenters_updated_original_mapped_grids(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 6);
end

%% STEP 8: Statistic 3 - Standard deviation in Z

input_points =  LiDAR_allPoints(:,1:3);

% The indices of the mapped grids are extracted and concatenated 
original_mapped_gridIndices_cell = gridIndices_cell_array(updated_original_mapped_grids); 

% Total number of mapped grids
total_mapped_grids = length(original_mapped_gridIndices_cell); 

% Standard deviations in orthogonal distances of points in the grid to
% plane
% standard_deviation_in_plane_orthogonals = zeros(total_mapped_grids,1); 
standard_deviation_in_z = zeros(total_mapped_grids,1); 

for ith_mapped_grid = 1:total_mapped_grids
    % points = input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:);
    % points = points(~isnan(points(:,1)),:);
    [~, standard_deviation_in_z(ith_mapped_grid,:), ~, ~, ~, ~] =...
    fcn_findEdge_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),-1);
    % mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
    % z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
end

fig_num_ENU_statistic_three = 80; 
figure(fig_num_ENU_statistic_three);clf;

hold on
axis on
grid on 
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers mapped grids in ENU')

% % plot(gridCenters_required_point_density(:,1), gridCenters_required_point_density(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);
% plot(gridCenters_drivable_grids(:,1), gridCenters_drivable_grids(:,2), '.','MarkerSize',50,'Color',[0 1 0]);
% plot(gridCenters_non_drivable_grids(:,1), gridCenters_non_drivable_grids(:,2), '.','MarkerSize',50,'Color',[1 0 0]);

plot(gridCenters_updated_original_mapped_grids(:,1), gridCenters_updated_original_mapped_grids(:,2), '.','MarkerSize',45,'Color',[0.2 0.2 0.2]);

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',2) % points strictly inside


for ith_text = 1:length(updated_current_mapped_grids(:,1))
    current_text = sprintf('%.0d',updated_current_mapped_grids(ith_text));
     % Place the text on the grid center
    text(gridCenters_updated_original_mapped_grids(ith_text,1), gridCenters_updated_original_mapped_grids(ith_text,2),current_text,'Color',[1 1 1],'HorizontalAlignment','center','FontSize', 6, 'FontWeight','bold');
end



% Plot grid lines and standard deviation
fig_num = 81;
figure(fig_num)

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Standard deviation of Z of mapped grids')

% Pre-allocation: To find the total number of points in each grid
total_points_in_mapped_grids = zeros(length(updated_original_mapped_grids),1);

% Pre-allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied
gridlines_mapped_grids = zeros(11*length(updated_original_mapped_grids),2); % length(gridlines) = 11

for ith_domain = 1:length(updated_original_mapped_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

    current_color = [0.5 0.5 0.5];
    % Plot current AABB
    current_AABB = grid_AABBs(updated_original_mapped_grids(ith_domain),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,1) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,2) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,3); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,4); ...
        current_AABB(1,2) current_AABB(1,4); ...
        ];

    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==original_mapped_grids(ith_domain);
    points_in_domain = input_points(rows_in_domain,:);

    total_points_in_mapped_grids(ith_domain) = length(points_in_domain);

    length_gridlines = length(gridlines);
    % % Plot the mapped points green
    % p1 = plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);

    % Plot the result
    % p1 = plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    % hold on
    % plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)
    gridlines_mapped_grids(1+(ith_domain-1)*length_gridlines:ith_domain*length_gridlines,:) = gridlines;

end

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',30,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

% [0.4660 0.6740 0.1880] - Green (Not so bright)
standard_deviation_in_z_round = round(standard_deviation_in_z,3);
% mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1));
for ith_text = 1:length(gridCenters_updated_original_mapped_grids(:,1))
    current_text = sprintf('%.3f',standard_deviation_in_z_round(ith_text));
    % Place the text on the grid center
    text(gridCenters_updated_original_mapped_grids(ith_text,1), gridCenters_updated_original_mapped_grids(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 10, 'FontWeight','bold','FontSmoothing','on');
end


% Find driven path indices in current mapped grids
driven_path_grid_indices_in_current_mapped_grids = ismember(updated_current_mapped_grids,current_grid_numbers_of_driven_path);

% Standard deviation in Z of driven path grids
std_in_z_driven_path = standard_deviation_in_z(driven_path_grid_indices_in_current_mapped_grids);

% Standard deviation in Z of other mapped grids
std_in_z_other_mapped_grids = standard_deviation_in_z(~driven_path_grid_indices_in_current_mapped_grids); 


% Plot grid lines and standard deviation
fig_num = 82;
figure(fig_num);clf;

hold on
grid on
xlabel('Mapped grid centers')
ylabel('Standard deviation in Z')
title('Mapped grid centers vs standard deviation in Z ')

plot(updated_current_mapped_grids, standard_deviation_in_z,'.','MarkerSize',30,'Color',[0.2 0.2 0.2])
plot(current_grid_numbers_of_driven_path, std_in_z_driven_path,'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',1.5)
plot(updated_current_mapped_grids(~driven_path_grid_indices_in_current_mapped_grids), std_in_z_other_mapped_grids,'.','MarkerSize',10,'Color',[1 0 0])


% Find mean std in z of driven path
mean_std_in_z_driven_path = mean(std_in_z_driven_path); 

% Find mean std in z of not driven path
mean_std_in_z_not_driven_path = mean(std_in_z_other_mapped_grids(~isnan(std_in_z_other_mapped_grids))); 

% Find max std in z of not driven path
max_std_in_z_not_driven_path = max(std_in_z_other_mapped_grids); 

% Std Threshold
% Instead of choosing 6, try to find a ratio
% ratio: mean_std_in_z_driven_path/mean_std_of_all_grids
% std_threshold = mean_std_in_z_driven_path*6;
% 
% disp(std_threshold)

%% Histogram of standard deviation

figure(123);clf;
hold on
grid on
xlabel('Standard deviation in z error after plane fit')
ylabel('Frequency')
title('Statistic 3: Determining suitable standard deviation in z')
% histogram(standard_deviation_in_z)
histogram(standard_deviation_in_z,100)
histogram(std_in_z_driven_path,5)

% std_threshold = mean_std_in_z_driven_path + 6*std(std_in_z_driven_path); 
% std_threshold = mean_std_in_z_driven_path + 3*std(std_in_z_driven_path); 
std_threshold = 0.1;%0.08; 
% plot(std_threshold,max(std_in_z_driven_path), 'k.', 'MarkerSize',20)
disp('mean of std_threshold of driven path')
disp(mean_std_in_z_driven_path)
disp('chosen std_threshold')
disp(std_threshold)
% std_threshold = 0.05; 
plot(std_threshold,0, 'k.', 'MarkerSize',18)
current_text = sprintf('std threshold = %.4f',std_threshold);
text(0.4, 80,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');

%% STEP 8: Statistic 4 - angle deviation

input_points = LiDAR_allPoints(:,1:3); 

% The indices of the mapped grids are extracted and concatenated 
original_mapped_gridIndices_cell = gridIndices_cell_array(updated_original_mapped_grids); 

% Total number of mapped grids
total_mapped_grids = length(original_mapped_gridIndices_cell); 

% Unit normal vectors of the plane fits of each mapped grid
unit_normal_vectors = zeros(total_mapped_grids,3); 


for ith_mapped_grid = 1:total_mapped_grids
    [~, ~, ~, unit_normal_vectors(ith_mapped_grid,:), ~, ~] =...
    fcn_findEdge_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),-1);
    % mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
    % z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
end

% STEP 2
% Comparing normal vector with vertical direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_normal_vectors.*unit_vector_vertical_direction,2);

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_normals_and_vertical = acos(dot_product);

fig_num_ENU_statistic_four = 83; 
figure(fig_num_ENU_statistic_four);clf;

hold on
axis on
grid on 
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers mapped grids in ENU')

% % plot(gridCenters_required_point_density(:,1), gridCenters_required_point_density(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);
% plot(gridCenters_drivable_grids(:,1), gridCenters_drivable_grids(:,2), '.','MarkerSize',50,'Color',[0 1 0]);
% plot(gridCenters_non_drivable_grids(:,1), gridCenters_non_drivable_grids(:,2), '.','MarkerSize',50,'Color',[1 0 0]);

plot(gridCenters_updated_original_mapped_grids(:,1), gridCenters_updated_original_mapped_grids(:,2), '.','MarkerSize',45,'Color',[0.2 0.2 0.2]);

for ith_text = 1:length(updated_current_mapped_grids(:,1))
    current_text = sprintf('%.0d',updated_current_mapped_grids(ith_text));
     % Place the text on the grid center
    text(gridCenters_updated_original_mapped_grids(ith_text,1), gridCenters_updated_original_mapped_grids(ith_text,2),current_text,'Color',[1 1 0],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',20,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

% Plot grid lines and standard deviation
fig_num = 84;
figure(fig_num)

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Angle deviation of mapped grids')

% Pre-allocation: To find the total number of points in each grid
total_points_in_mapped_grids = zeros(length(updated_original_mapped_grids),1);

% Pre-allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied
gridlines_mapped_grids = zeros(11*length(updated_original_mapped_grids),2); % length(gridlines) = 11

for ith_domain = 1:length(updated_original_mapped_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

    current_color = [0.5 0.5 0.5];
    % Plot current AABB
    current_AABB = grid_AABBs(updated_original_mapped_grids(ith_domain),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,1) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,2) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,3); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,4); ...
        current_AABB(1,2) current_AABB(1,4); ...
        ];

    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==updated_original_mapped_grids(ith_domain);
    points_in_domain = input_points(rows_in_domain,:);

    total_points_in_mapped_grids(ith_domain) = length(points_in_domain);

    length_gridlines = length(gridlines);
    % % Plot the mapped points green
    % p1 = plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);

    % Plot the result
    % p1 = plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    % hold on
    % plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)
    gridlines_mapped_grids(1+(ith_domain-1)*length_gridlines:ith_domain*length_gridlines,:) = gridlines;

end

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',30,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

% [0.4660 0.6740 0.1880] - Green (Not so bright)
angle_btw_unit_normals_and_vertical_round = round((angle_btw_unit_normals_and_vertical*180/pi),3);
% mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1));
for ith_text = 1:length(gridCenters_updated_original_mapped_grids(:,1))
    current_text = sprintf('%.3f',angle_btw_unit_normals_and_vertical_round(ith_text));
    % Place the text on the grid center
    text(gridCenters_updated_original_mapped_grids(ith_text,1), gridCenters_updated_original_mapped_grids(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 10, 'FontWeight','bold','FontSmoothing','on');
end

% Find driven path indices in current mapped grids
driven_path_grid_indices_in_current_mapped_grids = ismember(updated_current_mapped_grids,current_grid_numbers_of_driven_path);

% Standard deviation in Z of driven path grids
angle_btw_unit_normals_and_vertical_driven_path = angle_btw_unit_normals_and_vertical(driven_path_grid_indices_in_current_mapped_grids);

% Standard deviation in Z of other mapped grids
angle_btw_unit_normals_and_vertical_other_mapped_grids = angle_btw_unit_normals_and_vertical(~driven_path_grid_indices_in_current_mapped_grids); 


% Plot grid lines and standard deviation
fig_num = 85;
figure(fig_num)

hold on
grid on
xlabel('Mapped grid centers')
ylabel('Standard deviation in Z')
title('Mapped grid centers vs standard deviation in Z ')

plot(updated_current_mapped_grids, angle_btw_unit_normals_and_vertical*180/pi,'.','MarkerSize',30,'Color',[0.2 0.2 0.2])
plot(current_grid_numbers_of_driven_path, angle_btw_unit_normals_and_vertical_driven_path*180/pi,'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',1.5)
plot(updated_current_mapped_grids(~driven_path_grid_indices_in_current_mapped_grids), angle_btw_unit_normals_and_vertical_other_mapped_grids*180/pi,'.','MarkerSize',10,'Color',[1 0 0])


% Find mean std in z of driven path
mean_angle_btw_unit_normals_and_vertical_driven_path = mean(angle_btw_unit_normals_and_vertical_driven_path); 

% Find mean std in z of not driven path
max_angle_btw_unit_normals_and_vertical_not_driven_path = max(angle_btw_unit_normals_and_vertical_other_mapped_grids); 

% Theta threshold
% RATIO: Find a ratio
% theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 0.04;
% 
% disp(theta_threshold*(180/pi))

%% Histogram of angle deviation

figure(1223);clf;
hold on
grid on
xlabel('Angle deviation in z error after plane fit')
ylabel('Frequency')
title('Statistic 3: Determining suitable angle deviation')
% histogram(standard_deviation_in_z)
histogram(angle_btw_unit_normals_and_vertical,100)
histogram(angle_btw_unit_normals_and_vertical_driven_path,5)

% theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 3.7*std(angle_btw_unit_normals_and_vertical_driven_path); 
% theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 3.3*std(angle_btw_unit_normals_and_vertical_driven_path); 
% theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 3*std(angle_btw_unit_normals_and_vertical_driven_path);

% plot(theta_threshold,max(angle_btw_unit_normals_and_vertical), 'k.', 'MarkerSize',20)
theta_threshold = 0.1745; 
% theta_threshold = 0.1504; 
% theta_threshold = 0.3; 

disp('mean of angle_deviation_driven_path')
disp(mean_angle_btw_unit_normals_and_vertical_driven_path*(180/pi))
disp('theta threshold')
disp(theta_threshold*(180/pi))

plot(theta_threshold,0, 'k.', 'MarkerSize',18)
current_text = sprintf('theta threshold = %.1f',theta_threshold*(180/pi));
text(0.5, 30,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');

%% STEP 9: 3rd Classification
%Jiabao 
input_points = LiDAR_allPoints(:,1:3); 
% std_threshold = 0.05; 
% theta_threshold = 7*pi/180;
% theta_threshold = 30*pi/180;
% gridCenters

fig_num = 6000; 
figure(fig_num);clf

% theta_threshold = [];
% std_threshold = [];

% % Classify mapped grids into drivable and drivable
% [standard_deviation_in_z, angle_btw_unit_normals_and_vertical, ...
%     original_drivable_grids, original_non_drivable_grids, current_drivable_grid_numbers_in_mapped_grids, current_non_drivable_grid_numbers_in_mapped_grids, ...
%     gridCenters_drivable_grids,gridCenters_non_drivable_grids, concatenate_gridCenters_drivable_non_drivable_grids] = ...
%     fcn_geometry_classifyGridsAsDrivable(gridIndices_cell_array, original_mapped_grids, input_points, std_threshold, theta_threshold, gridCenters, fig_num);


% OLD: fixed by SNB on 2024_08_05 to correctly fill right variables
% % Classify mapped grids into drivable and drivable
% [standard_deviation_in_z, angle_btw_unit_normals_and_vertical, ...
%     original_drivable_grids, original_non_drivable_grids, current_drivable_grid_numbers_in_mapped_grids, current_non_drivable_grid_numbers_in_mapped_grids, ...
%     gridCenters_drivable_grids,gridCenters_non_drivable_grids, concatenate_gridCenters_drivable_non_drivable_grids] = ...
%     fcn_geometry_classifyGridsAsDrivable(gridIndices_cell_array, updated_original_mapped_grids, input_points, std_threshold, theta_threshold, gridCenters, fig_num);


[standard_deviation_in_z, ...
    angle_btw_unit_normals_and_vertical, ...
    original_drivable_grids, ...
    original_non_drivable_grids, ...
    current_drivable_grid_numbers_in_mapped_grids, ...
    current_non_drivable_grid_numbers_in_mapped_grids, ...
    ~, .... % current_failed_grid_numbers_in_mapped_grids, ...
    ~, ...  % current_uncertain_grid_numbers_in_mapped_grids, ...
    ~, ...  % gridCenters_failed_grids, ...
    ~, ...  % gridCenters_uncertain_grids, ...
    gridCenters_drivable_grids, ...
    gridCenters_non_drivable_grids, ...
    concatenate_gridCenters_drivable_non_drivable_grids] = ...
    fcn_findEdge_classifyGridsAsDrivable(gridIndices_cell_array, updated_original_mapped_grids, input_points, std_threshold, theta_threshold, gridCenters, fig_num);

std_threshold_failed_gridCenters = gridCenters_updated_original_mapped_grids(standard_deviation_in_z>std_threshold,:); 

theta_threshold_failed_gridCenters = gridCenters_updated_original_mapped_grids(angle_btw_unit_normals_and_vertical>theta_threshold,:); 

% plot(std_threshold_failed_gridCenters(:,1), std_threshold_failed_gridCenters(:,2), 'o','MarkerSize',20,'Color',[0 1 1], 'LineWidth',2) 

% fig_num = 98989; 
% figure(fig_num); clf; 

marker_size = 10;
RGB_triplet = [1, 1, 0]; 
legend_option = 1;
legend_name = 'Std threshold failed grids';
marker_type = 'o';
legend_position = [];
plot_std_threshold_failed_gridCenters = [std_threshold_failed_gridCenters, zeros(length(std_threshold_failed_gridCenters(:,1)),1)]; 
[LLA_data_std_threshold_failed_grids] = fcn_findEdge_plotPointsinLLA(plot_std_threshold_failed_gridCenters,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
% fcn_geometry_plotPointsinLLA(plot_gridCenters_updated_original_unmapped_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% fcn_geometry_plotPointsinLLA(plot_gridCenters_non_drivable_grids,marker_size,RGB_triplet,[],legend_option,legend_name,[],[],[],[],fig_num);

% geoplot(LLA_data_std_threshold_failed_grids(:,1), LLA_data_std_threshold_failed_grids(:,2), 'o','MarkerSize',10,'Color',[1 1 0], 'LineWidth',2) 

% plot grid centers
marker_size = 15;
RGB_triplet = [1, 0, 1]; 
legend_option = 1;
legend_name = 'Theta threshold failed grids';
marker_type = 'o';
legend_position = [];
plot_theta_threshold_failed_gridCenters = [theta_threshold_failed_gridCenters, zeros(length(theta_threshold_failed_gridCenters(:,1)),1)]; 
[LLA_data_theta_threshold_failed_grids] = fcn_findEdge_plotPointsinLLA(plot_theta_threshold_failed_gridCenters,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% geoplot(LLA_data_theta_threshold_failed_grids(:,1), LLA_data_theta_threshold_failed_grids(:,2), 'o','MarkerSize',15,'Color',[1 0 1], 'LineWidth',2) 

%% STEP 10: Find the boundary points of drivable and non-drivable grids
% Dr B
% Part1 - Find the boundary points of mapped and unmapped grids

% Revision History
% Funtionalized this code
% Added plotting options

% INPUTS - gridCenters_low_point_density,
% gridCenters_required_point_density, figure num
% OUTPUTS - X, Y, Z 

fig_num_mapped_unmapped = 767787; 

XYZ_matrix_mapped_grids = [gridCenters_updated_original_mapped_grids(:,1:2) ones(length(gridCenters_updated_original_mapped_grids(:,1)),1)]; 

XYZ_matrix_unmapped_grids = [gridCenters_updated_original_unmapped_grids(:,1:2) zeros(length(gridCenters_updated_original_unmapped_grids(:,1)),1)]; 

XYZ_matrix_mapped_unmapped_gridcenters = [XYZ_matrix_mapped_grids; XYZ_matrix_unmapped_grids]; 

% Find the unique elements
XYZ_matrix = unique(XYZ_matrix_mapped_unmapped_gridcenters,'rows'); 

% Reverse the matrix
XYZ_matrix_reversed = flipud(XYZ_matrix);

% Find unique rows based on the first two columns (X and Y)
[~,XYZ_matrix_indices] = unique(XYZ_matrix_reversed(:,1:2),'rows'); 

% Create the new matrix with only unique rows, keeping the last occurrence of each duplicate
XYZ_matrix_reversed = XYZ_matrix_reversed(XYZ_matrix_indices,:);  

% Reverse the matrix back to the original order
XYZ_matrix = flipud(XYZ_matrix_reversed);

XYZ_matrix = round(XYZ_matrix,4);

x_range = unique(XYZ_matrix(:,1))'; 
y_range = unique(XYZ_matrix(:,2))';


% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = NaN(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

XY_pairs = round(XY_pairs,4); 

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% C = setdiff(XY_pairs, XY_pairs(idx==0,:), 'rows');

% Fill Z matrix with corresponding Z values
Z(idx~=0) = XYZ_matrix(idx(idx~=0), 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

%%%%%%%%%%%%%%---------------------------------------------------------------------------
% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
figure(fig_num_mapped_unmapped);
clf;
boundary_points_mapped_unmapped = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num_mapped_unmapped);

% assert(length(drivable_grids)>=1)
% assert(length(non_drivable_grids)>=1)
% % assert(length(unmapped_grids)==zeros(0,1))
% assert(length(gridCenters_mapped_grids(:,1))>=1)
% assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(boundary_points_mapped_unmapped)>=1)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundary points of drivable and non-drivable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 2 - Find boundary points of drivable and non-drivable grids

fig_num_drivable_non_drivable = 98898;

% XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)];
XYZ_matrix = concatenate_gridCenters_drivable_non_drivable_grids;

XYZ_matrix = unique(XYZ_matrix,'rows'); 

% Reverse the matrix
XYZ_matrix_reversed = flipud(XYZ_matrix);

% Find unique rows based on the first two columns (X and Y)
[~,XYZ_matrix_indices] = unique(XYZ_matrix_reversed(:,1:2),'rows'); 

% Create the new matrix with only unique rows, keeping the last occurrence of each duplicate
XYZ_matrix_reversed = XYZ_matrix_reversed(XYZ_matrix_indices,:);  

% Reverse the matrix back to the original order
XYZ_matrix = flipud(XYZ_matrix_reversed);

XYZ_matrix = round(XYZ_matrix,4); 
% Given x_range and y_range
% x_range = min(XYZ_matrix(:,1)):grid_size:max(XYZ_matrix(:,1)); % min_x:gridSize:max_x 
% y_range = min(XYZ_matrix(:,2)):grid_size:max(XYZ_matrix(:,2)); % min_y:gridSize:max_y 

x_range = unique(XYZ_matrix(:,1))'; 
y_range = unique(XYZ_matrix(:,2))';
% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = NaN(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

XY_pairs = round(XY_pairs,4); 

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% C = setdiff(XY_pairs, XY_pairs(idx==0,:), 'rows');

% Fill Z matrix with corresponding Z values
Z(idx~=0) = XYZ_matrix(idx(idx~=0), 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points

figure(fig_num_drivable_non_drivable)
clf;
boundary_points = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num_drivable_non_drivable);

% assert(length(drivable_grids)>=1)
% assert(length(non_drivable_grids)>=1)
% % assert(length(unmapped_grids)==zeros(0,1))
% assert(length(gridCenters_mapped_grids(:,1))>=1)
% assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(boundary_points)>=1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finding true boundary points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 3 - Find true boundary points by eliminating the boundary points of
% mapped and unmapped from drivable and non-drivable boubndary points. 
% Jiabao
fig_num_bd_pts_ENU = 1000; 

[members, id_x] = ismember(boundary_points,boundary_points_mapped_unmapped,'rows'); 

not_boundary_points = boundary_points(members,:);

true_boundary_points = boundary_points(members==0,:);

figure(fig_num_bd_pts_ENU)
clf;
% figure(10001)
hold on
plot(true_boundary_points(:,1), true_boundary_points(:,2), 'c.', 'MarkerSize',40)
plot(true_boundary_points(:,1), true_boundary_points(:,2), 'b.', 'MarkerSize',30)

fig_num = 6000;
% plot computed boundary points
marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 1;
legend_name = 'Computed boundary points';
legend_position = [];
marker_type = [];
plot_true_boundary_points = [true_boundary_points, zeros(length(true_boundary_points),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot computed boundary points
marker_size = 10;
RGB_triplet = [0 0 1]; 
legend_option = 0;
legend_name = 'Computed boundary points';
legend_position = [];
marker_type = [];
[~] = fcn_findEdge_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot driven path
marker_size = 30;
RGB_triplet = [0 0 1]; 
legend_option = 0;
legend_name = 'Driven path grids';
legend_position = [];
marker_type = [];
plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% Find the nearest boundary points in ENU - updated

% grid_width = ceil(max(true_boundary_points(:,1)) - min(true_boundary_points(:,1))); 
% grid_height = ceil(max(true_boundary_points(:,2)) - min(true_boundary_points(:,2))); 
% plot computed boundary points

fig_num = 6020;
% plot computed boundary points
marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 1;
legend_name = 'Computed boundary points';
legend_position = [];
plot_true_boundary_points = [true_boundary_points, zeros(length(true_boundary_points),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
marker_size = 10;
RGB_triplet = [0 0 1]; 
legend_option = 0;
legend_name = 'Computed boundary points';
legend_position = [];

[~] = fcn_findEdge_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


%% plot the true borders
% Dr B
fig_num = 1224; 
figure(fig_num); clf;

% after running script_test_geometry_updatedSurfaceAnalysis
% [true_borders,true_borders_x,true_borders_y] = fcn_geometry_findNearestBoundaryPoints(true_boundary_points,...
%      gridCenters_driven_path, fig_num);

[~, ~, true_borders]  = fcn_findEdge_findNearestBoundaryPoints(true_boundary_points, ...
     gridCenters_driven_path, grid_size, grid_boundaries, fig_num);


%%


fig_num = 87909; 
figure(fig_num); clf; 

marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 0;
legend_name = 'Nearest Boundary Points';
legend_position = [];
marker_type = []; 
nearest_boundary_points = [ true_borders(:,1), true_borders(:,2), zeros(length(true_borders),1)]; 

[~] = fcn_findEdge_plotPointsinLLA(nearest_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

marker_size = 10;
RGB_triplet = [0 0 1]; 
legend_option = 0;
legend_name = 'Nearest Boundary Points';
legend_position = [];
marker_type = []; 
nearest_boundary_points = [true_borders(:,1), true_borders(:,2),  zeros(length(true_borders),1)]; 


[~] = fcn_findEdge_plotPointsinLLA(nearest_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


%% Find the nearest boundary points in ENU 

% Instructions
% Run this script, after running the scripts in STEP 1 and STEP 2. 


% Write the grid number at the grid center for reference. 

fig_num = 3332; 
figure(fig_num);clf;

hold on
axis on
grid on 
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers of drivable and non-drivable grids in ENU')

% gridCenters_drivable_grids
% 
% gridCenters_non_drivable_grids

% plot(gridCenters_mapped_grids(:,1), )

% plot(gridCenters_mapped_grids(:,1), gridCenters_mapped_grids(:,2), '.','MarkerSize',45,'Color',[0.2 0.2 0.2]);
p1 = plot(gridCenters_drivable_grids(:,1), gridCenters_drivable_grids(:,2), '.','MarkerSize',45,'Color',[0.4660 0.6740 0.1880],'DisplayName','Drivable grids');
p2 = plot(gridCenters_non_drivable_grids(:,1), gridCenters_non_drivable_grids(:,2), '.','MarkerSize',45,'Color',[0.6350 0.0780 0.1840], 'DisplayName','Non-drivable grids');

% for ith_text = 1:length(original_mapped_grids(:,1))
%     current_text = sprintf('%.0d',original_mapped_grids(ith_text));
% 
%     % Place the text on the grid center
%     text(gridCenters_mapped_grids(ith_text,1), gridCenters_mapped_grids(ith_text,2),current_text,'Color',[1, 1, 0],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
% end

% plot true boundary points
plot(nearest_boundary_points(:,1), nearest_boundary_points(:,2), 'c.', 'MarkerSize',40, 'DisplayName','Boundary points');
p3 = plot(nearest_boundary_points(:,1), nearest_boundary_points(:,2), 'b.', 'MarkerSize',30, 'DisplayName','Boundary points');

% % plot the grids in the driven path
p4 = plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',20,'Color',[0 1 0], 'LineWidth',2, 'DisplayName','Driven path grids'); % points strictly inside
% % plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',20,'Color',[0 0 0], 'LineWidth',0.5) % points strictly inside

%% Right boundary points

VehiclePose_current = VehiclePose(scanLineNumber_start:scanLineNumber_end,1:2); 
% Calculate the vehicle orientation
vehicle_change_in_pose_XY = diff(VehiclePose_current(:,1:2));

% Repeat the last value again, since diff removes one row. We want the same
% number of vectors as the number of points, and diff removed one point.
vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];

% Convert these to unit vectors
unit_vehicle_change_in_pose_XY = fcn_INTERNAL_calcUnitVector(vehicle_change_in_pose_XY);

shift = 5; 
% Shift
shift_distance = unit_vehicle_change_in_pose_XY*shift; 

unit_vehicle_change_in_pose_XY_shifted = unit_vehicle_change_in_pose_XY - shift_distance;


% Shift the vehicle pose
VehiclePose_current_shifted = VehiclePose_current - shift_distance; 

% Get the number of rows 
rows_nearest_boundary_points = size(nearest_boundary_points, 1);
rows_VehiclePose_current_shifted = size(VehiclePose_current_shifted, 1);

% Calculate how many rows need to be added to VehiclePose_current_shifted
rows_to_add = rows_nearest_boundary_points - rows_VehiclePose_current_shifted;


if 0 <= rows_to_add

    % Create the additional rows by repeating the last row of VehiclePose_current_shifted
    additional_rows = repmat(VehiclePose_current_shifted(rows_VehiclePose_current_shifted, :), rows_to_add, 1);

    % Concatenate the original VehiclePose_current_shifted and the additional rows
    updated_VehiclePose_current_shifted = [VehiclePose_current_shifted; additional_rows];
    
else
    % Find the total number of rows required 
    current_rows = length(VehiclePose_current_shifted(:,1)) + rows_to_add; 
    
    % Remove some rows from the original VehiclePose_current_shifted 
    updated_VehiclePose_current_shifted = VehiclePose_current_shifted(1:current_rows,:); 

end

% Calculate the vectors
vector_from_vehicle_pose_to_boundary_points = nearest_boundary_points(:,1:2) - updated_VehiclePose_current_shifted;

% Find orthogonal vetors by rotating by 90 degrees in the CCW direction
unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY_shifted*[0 1; -1 0];


if 0 <= rows_to_add

    % Create the additional rows by repeating the last row of unit_ortho_vehicle_vectors_XY
    additional_rows_unit_ortho_vehicle_vectors_XY = repmat(unit_ortho_vehicle_vectors_XY(rows_VehiclePose_current_shifted, :), rows_to_add, 1);

    % Concatenate the original unit_ortho_vehicle_vectors_XY and the additional_rows_unit_ortho_vehicle_vectors_XY
    updated_unit_ortho_vehicle_vectors_XY = [unit_ortho_vehicle_vectors_XY; additional_rows_unit_ortho_vehicle_vectors_XY];

else

    % Find the total number of rows required 
    current_rows = length(VehiclePose_current_shifted(:,1)) + rows_to_add; 
    
    % Remove some rows from the original VehiclePose_current_shifted 
    updated_unit_ortho_vehicle_vectors_XY = unit_ortho_vehicle_vectors_XY(1:current_rows,:); 

end

% Calculate the transverse distance
transverse_dist_boundary_points = sum(vector_from_vehicle_pose_to_boundary_points.*updated_unit_ortho_vehicle_vectors_XY,2);

% Transverse distances of the right boundaries
transverse_dist_right_boundary_points = transverse_dist_boundary_points(transverse_dist_boundary_points>0,:);

% Boundary points on the right
boundary_points_right = nearest_boundary_points(transverse_dist_boundary_points>0,:);

figure(31543);clf; 
hold on
axis on
grid on 
xlabel('X[m]')
ylabel('Y[m]')
title('Right Points')

plot(VehiclePose_current_shifted(:,1), VehiclePose_current_shifted(:,2),'.','Color',[0 0 0],'MarkerSize',30) 

plot(nearest_boundary_points(:,1), nearest_boundary_points(:,2), 'c.', 'MarkerSize',40, 'DisplayName','Boundary points');
plot(nearest_boundary_points(:,1), nearest_boundary_points(:,2), 'b.', 'MarkerSize',30, 'DisplayName','Boundary points');

plot(boundary_points_right(:,1), boundary_points_right(:,2), 'ro', 'MarkerSize',30, 'DisplayName','Boundary points');

quiver(...
    VehiclePose_current(:,1),VehiclePose_current(:,2),...
    unit_vehicle_change_in_pose_XY(:,1),unit_vehicle_change_in_pose_XY(:,2),'-','LineWidth',3,'Color',[0 1 0]);

quiver(...
    VehiclePose_current(:,1),VehiclePose_current(:,2),...
    unit_ortho_vehicle_vectors_XY(:,1),unit_ortho_vehicle_vectors_XY(:,2),'-','LineWidth',3,'Color',[0 0 1]);

% % Unit orthogonal vector
% repeated_unit_ortho_vehicle_vectors_XYZ = unit_ortho_vehicle_vectors_XY(28,:).*(ones(length(nearest_boundary_points(:,1)),2)); 

%% Plot right boundary points

fig_num = 71959; 
figure(fig_num); clf; 

marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 0;
legend_name = 'Right Boundary Points';
legend_position = [];
marker_type = []; 


[~] = fcn_findEdge_plotPointsinLLA(boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

marker_size = 12;
RGB_triplet = [1 0 0]; 
legend_option = 0;
legend_name = 'Right Boundary Points';
legend_position = [];
marker_type = []; 


[~] = fcn_findEdge_plotPointsinLLA(boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%%
% %% Remove boundary points on the left
% 
% % scanLineNumber_start = 1400; 
% % scanLineNumber_end = 1410; 
% 
% boundaryLineNumber_start = scanLineNumber_start - 8;
% boundaryLineNumber_end = scanLineNumber_end - 6; 
% 
% VehiclePose_current = VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1:2); 
% % Calculate the vehicle orientation
% vehicle_change_in_pose_XY = diff(VehiclePose_current(:,1:2));
% 
% % Repeat the last value again, since diff removes one row. We want the same
% % number of vectors as the number of points, and diff removed one point.
% vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];
% 
% % Convert these to unit vectors
% unit_vehicle_change_in_pose_XY = fcn_INTERNAL_calcUnitVector(vehicle_change_in_pose_XY);
% 
% % Find orthogonal vetors by rotating by 90 degrees in the CCW direction
% unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY*[0 1; -1 0];
% 
% % Calculate the vectors
% vector_from_vehicle_pose_to_boundary_points = true_boundary_points - VehiclePose_current(28,:).*(ones(length(true_boundary_points(:,1)),2)); 
% 
% % Unit orthogonal vector
% repeated_unit_ortho_vehicle_vectors_XYZ = unit_ortho_vehicle_vectors_XY(28,:).*(ones(length(true_boundary_points(:,1)),2)); 
% 
% 
% % Calculate the transverse distance
% transverse_dist_boundary_points = sum(vector_from_vehicle_pose_to_boundary_points.*repeated_unit_ortho_vehicle_vectors_XYZ,2);
% 
% % Transverse distances of the right boundaries
% transverse_dist_right_boundary_points = transverse_dist_boundary_points(transverse_dist_boundary_points<0,:);
% 
% % Boundary points on the right
% boundary_points_right = true_boundary_points(transverse_dist_boundary_points<0,:);
% 
% % nearest boundary points to the right
% nearest_boundary_points_right = boundary_points_right(abs(transverse_dist_right_boundary_points)<2.8,:);
% 
% figure(3543);clf; 
% hold on
% axis on
% grid on 
% xlabel('X[m]')
% ylabel('Y[m]')
% title('Right Points')
% 
% plot(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1), VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,2),'.','Color',[0 0 0],'MarkerSize',30) 
% 
% plot(true_boundary_points(:,1), true_boundary_points(:,2), 'c.', 'MarkerSize',40, 'DisplayName','Boundary points');
% plot(true_boundary_points(:,1), true_boundary_points(:,2), 'b.', 'MarkerSize',30, 'DisplayName','Boundary points');
% 
% plot(nearest_boundary_points_right(:,1), nearest_boundary_points_right(:,2), 'ro', 'MarkerSize',30, 'DisplayName','Boundary points');
% 
% 
% quiver(...
%     VehiclePose_current(:,1),VehiclePose_current(:,2),...
%     unit_vehicle_change_in_pose_XY(:,1),unit_vehicle_change_in_pose_XY(:,2),'-','LineWidth',3,'Color',[0 1 0]);
% 
% quiver(...
%     VehiclePose_current(:,1),VehiclePose_current(:,2),...
%     unit_ortho_vehicle_vectors_XY(:,1),unit_ortho_vehicle_vectors_XY(:,2),'-','LineWidth',3,'Color',[0 0 1]);
% 
% % quiver3(...
% %     VehiclePose(scanLineNumber_start,1),VehiclePose(scanLineNumber_start,2),VehiclePose(scanLineNumber_start,3), ...
% %     unit_ortho_vehicle_vectors_XY(scanLineNumber_start,1),unit_ortho_vehicle_vectors_XY(scanLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 0 1]);
% 
% 
% % 
% %% plot boundary points to the right
% 
% 
% % grid_width = ceil(max(true_boundary_points(:,1)) - min(true_boundary_points(:,1))); 
% % grid_height = ceil(max(true_boundary_points(:,2)) - min(true_boundary_points(:,2))); 
% % plot computed boundary points
% 
% fig_num = 70067;clf;
% % plot computed boundary points
% marker_size = 30;
% RGB_triplet = [0 1 1]; 
% legend_option = 1;
% legend_name = 'Computed boundary points';
% legend_position = [];
% plot_nearest_boundary_points_right = [nearest_boundary_points_right, zeros(length(nearest_boundary_points_right),1)];
% [~] = fcn_geometry_plotGridCenters(plot_nearest_boundary_points_right,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);
% 
% marker_size = 10;
% RGB_triplet = [0 0 1]; 
% legend_option = 0;
% legend_name = 'Computed boundary points';
% legend_position = [];
% 
% [~] = fcn_geometry_plotGridCenters(plot_nearest_boundary_points_right,marker_size,RGB_triplet,legend_option,legend_name,legend_position,fig_num);


%%

% oriGIN = [0 0 0]; 
% unit_vertical = [0 0 1]; 
% input_points = LiDAR_allPoints;
% total_non_drivable_grids = length(original_non_drivable_grids);
% for ith_mapped_grid = 7
%     figure(ith_mapped_grid*10000); clf;
% 
%     [~, standard_deviation_in_z(ith_mapped_grid,:), ~, unit_normal_vectors, base_point, ~] =...
%         fcn_geometry_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),ith_mapped_grid*10000);
%     xlabel('X[m]')
%     ylabel('Y[m]')
%     zlabel('Z[m]')
%      title('Plane fit of a mapped grid')
%     % title(sprintf('Plane fit of  %dth mapped grid',ith_mapped_grid))
%     figure(ith_mapped_grid*11100)
%     hold on;
%     grid on;
%     xlabel('X[m]')
%     ylabel('Y[m]')
%     zlabel('Z[m]')
%     title(sprintf('angle btw %dth mapped grid (Non-drivable)',ith_mapped_grid))
%     % Plot the base point
%     view(3)
%     % plot3(base_point(1,1),base_point(1,2),base_point(1,3),'r.','MarkerSize',50);
%     quiver3(oriGIN(1,1),oriGIN(1,2),oriGIN(1,3), unit_normal_vectors(1,1),unit_normal_vectors(1,2),unit_normal_vectors(1,3),0,'g','Linewidth',3);
%     quiver3(oriGIN(1,1),oriGIN(1,2),oriGIN(1,3), unit_vertical(1,1),unit_vertical(1,2),unit_vertical(1,3),0,'r','Linewidth',3);
%     % mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
% %     % z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
% end



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

%% function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
% Clear out the variables
clear global flag* FLAG*
clear flag*
clear path

% Clear out any path directories under Utilities
path_dirs = regexp(path,'[;]','split');
utilities_dir = fullfile(pwd,filesep,'Utilities');
for ith_dir = 1:length(path_dirs)
    utility_flag = strfind(path_dirs{ith_dir},utilities_dir);
    if ~isempty(utility_flag)
        rmpath(path_dirs{ith_dir});
    end
end

% Delete the Utilities folder, to be extra clean!
if  exist(utilities_dir,'dir')
    [status,message,message_ID] = rmdir(utilities_dir,'s');
    if 0==status
        error('Unable remove directory: %s \nReason message: %s \nand message_ID: %s\n',utilities_dir, message,message_ID);
    end
end

end % Ends fcn_INTERNAL_clearUtilitiesFromPathAndFolders

%% fcn_INTERNAL_initializeUtilities
function  fcn_INTERNAL_initializeUtilities(library_name,library_folders,library_url,this_project_folders)
% Reset all flags for installs to empty
clear global FLAG*

fprintf(1,'Installing utilities necessary for code ...\n');

% Dependencies and Setup of the Code
% This code depends on several other libraries of codes that contain
% commonly used functions. We check to see if these libraries are installed
% into our "Utilities" folder, and if not, we install them and then set a
% flag to not install them again.

% Set up libraries
for ith_library = 1:length(library_name)
    dependency_name = library_name{ith_library};
    dependency_subfolders = library_folders{ith_library};
    dependency_url = library_url{ith_library};
    
    fprintf(1,'\tAdding library: %s ...',dependency_name);
    fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url);
    clear dependency_name dependency_subfolders dependency_url
    fprintf(1,'Done.\n');
end

% Set dependencies for this project specifically
fcn_DebugTools_addSubdirectoriesToPath(pwd,this_project_folders);

disp('Done setting up libraries, adding each to MATLAB path, and adding current repo folders to path.');
end % Ends fcn_INTERNAL_initializeUtilities


function fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url, varargin)
%% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES - MATLAB package installer from URL
%
% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES installs code packages that are
% specified by a URL pointing to a zip file into a default local subfolder,
% "Utilities", under the root folder. It also adds either the package
% subfoder or any specified sub-subfolders to the MATLAB path.
%
% If the Utilities folder does not exist, it is created.
%
% If the specified code package folder and all subfolders already exist,
% the package is not installed. Otherwise, the folders are created as
% needed, and the package is installed.
%
% If one does not wish to put these codes in different directories, the
% function can be easily modified with strings specifying the
% desired install location.
%
% For path creation, if the "DebugTools" package is being installed, the
% code installs the package, then shifts temporarily into the package to
% complete the path definitions for MATLAB. If the DebugTools is not
% already installed, an error is thrown as these tools are needed for the
% path creation.
%
% Finally, the code sets a global flag to indicate that the folders are
% initialized so that, in this session, if the code is called again the
% folders will not be installed. This global flag can be overwritten by an
% optional flag input.
%
% FORMAT:
%
%      fcn_DebugTools_installDependencies(...
%           dependency_name, ...
%           dependency_subfolders, ...
%           dependency_url)
%
% INPUTS:
%
%      dependency_name: the name given to the subfolder in the Utilities
%      directory for the package install
%
%      dependency_subfolders: in addition to the package subfoder, a list
%      of any specified sub-subfolders to the MATLAB path. Leave blank to
%      add only the package subfolder to the path. See the example below.
%
%      dependency_url: the URL pointing to the code package.
%
%      (OPTIONAL INPUTS)
%      flag_force_creation: if any value other than zero, forces the
%      install to occur even if the global flag is set.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      This code will automatically get dependent files from the internet,
%      but of course this requires an internet connection. If the
%      DebugTools are being installed, it does not require any other
%      functions. But for other packages, it uses the following from the
%      DebugTools library: fcn_DebugTools_addSubdirectoriesToPath
%
% EXAMPLES:
%
% % Define the name of subfolder to be created in "Utilities" subfolder
% dependency_name = 'DebugTools_v2023_01_18';
%
% % Define sub-subfolders that are in the code package that also need to be
% % added to the MATLAB path after install; the package install subfolder
% % is NOT added to path. OR: Leave empty ({}) to only add
% % the subfolder path without any sub-subfolder path additions.
% dependency_subfolders = {'Functions','Data'};
%
% % Define a universal resource locator (URL) pointing to the zip file to
% % install. For example, here is the zip file location to the Debugtools
% % package on GitHub:
% dependency_url = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_18.zip?raw=true';
%
% % Call the function to do the install
% fcn_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url)
%
% This function was written on 2023_01_23 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_01_23:
% -- wrote the code originally
% 2023_04_20:
% -- improved error handling
% -- fixes nested installs automatically

% TO DO
% -- Add input argument checking

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_check_inputs
    % Are there the right number of inputs?
    narginchk(3,4);
end

%% Set the global variable - need this for input checking
% Create a variable name for our flag. Stylistically, global variables are
% usually all caps.
flag_varname = upper(cat(2,'flag_',dependency_name,'_Folders_Initialized'));

% Make the variable global
eval(sprintf('global %s',flag_varname));

if nargin==4
    if varargin{1}
        eval(sprintf('clear global %s',flag_varname));
    end
end

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~exist(flag_varname,'var') || isempty(eval(flag_varname))
    % Save the root directory, so we can get back to it after some of the
    % operations below. We use the Print Working Directory command (pwd) to
    % do this. Note: this command is from Unix/Linux world, but is so
    % useful that MATLAB made their own!
    root_directory_name = pwd;
    
    % Does the directory "Utilities" exist?
    utilities_folder_name = fullfile(root_directory_name,'Utilities');
    if ~exist(utilities_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(root_directory_name,'Utilities');
        
        % Did it work?
        if ~success_flag
            error('Unable to make the Utilities directory. Reason: %s with message ID: %s\n',error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The Utilities directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',error_message, message_ID);
        end
        
    end
    
    % Does the directory for the dependency folder exist?
    dependency_folder_name = fullfile(root_directory_name,'Utilities',dependency_name);
    if ~exist(dependency_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(utilities_folder_name,dependency_name);
        
        % Did it work?
        if ~success_flag
            error('Unable to make the dependency directory: %s. Reason: %s with message ID: %s\n',dependency_name, error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The %s directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',dependency_name, error_message, message_ID);
        end
        
    end
    
    % Do the subfolders exist?
    flag_allFoldersThere = 1;
    if isempty(dependency_subfolders{1})
        flag_allFoldersThere = 0;
    else
        for ith_folder = 1:length(dependency_subfolders)
            subfolder_name = dependency_subfolders{ith_folder};
            
            % Create the entire path
            subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
            
            % Check if the folder and file exists that is typically created when
            % unzipping.
            if ~exist(subfunction_folder,'dir')
                flag_allFoldersThere = 0;
            end
        end
    end
    
    % Do we need to unzip the files?
    if flag_allFoldersThere==0
        % Files do not exist yet - try unzipping them.
        save_file_name = tempname(root_directory_name);
        zip_file_name = websave(save_file_name,dependency_url);
        % CANT GET THIS TO WORK --> unzip(zip_file_url, debugTools_folder_name);
        
        % Is the file there?
        if ~exist(zip_file_name,'file')
            error(['The zip file: %s for dependency: %s did not download correctly.\n' ...
                'This is usually because permissions are restricted on ' ...
                'the current directory. Check the code install ' ...
                '(see README.md) and try again.\n'],zip_file_name, dependency_name);
        end
        
        % Try unzipping
        unzip(zip_file_name, dependency_folder_name);
        
        % Did this work? If so, directory should not be empty
        directory_contents = dir(dependency_folder_name);
        if isempty(directory_contents)
            error(['The necessary dependency: %s has an error in install ' ...
                'where the zip file downloaded correctly, ' ...
                'but the unzip operation did not put any content ' ...
                'into the correct folder. ' ...
                'This suggests a bad zip file or permissions error ' ...
                'on the local computer.\n'],dependency_name);
        end
        
        % Check if is a nested install (for example, installing a folder
        % "Toolsets" under a folder called "Toolsets"). This can be found
        % if there's a folder whose name contains the dependency_name
        flag_is_nested_install = 0;
        for ith_entry = 1:length(directory_contents)
            if contains(directory_contents(ith_entry).name,dependency_name)
                if directory_contents(ith_entry).isdir
                    flag_is_nested_install = 1;
                    install_directory_from = fullfile(directory_contents(ith_entry).folder,directory_contents(ith_entry).name);
                    install_files_from = fullfile(directory_contents(ith_entry).folder,directory_contents(ith_entry).name,'*'); % BUG FIX - For Macs, must be *, not *.*
                    install_location_to = fullfile(directory_contents(ith_entry).folder);
                end
            end
        end
        
        if flag_is_nested_install
            [status,message,message_ID] = movefile(install_files_from,install_location_to);
            if 0==status
                error(['Unable to move files from directory: %s\n ' ...
                    'To: %s \n' ...
                    'Reason message: %s\n' ...
                    'And message_ID: %s\n'],install_files_from,install_location_to, message,message_ID);
            end
            [status,message,message_ID] = rmdir(install_directory_from);
            if 0==status
                error(['Unable remove directory: %s \n' ...
                    'Reason message: %s \n' ...
                    'And message_ID: %s\n'],install_directory_from,message,message_ID);
            end
        end
        
        % Make sure the subfolders were created
        flag_allFoldersThere = 1;
        if ~isempty(dependency_subfolders{1})
            for ith_folder = 1:length(dependency_subfolders)
                subfolder_name = dependency_subfolders{ith_folder};
                
                % Create the entire path
                subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
                
                % Check if the folder and file exists that is typically created when
                % unzipping.
                if ~exist(subfunction_folder,'dir')
                    flag_allFoldersThere = 0;
                end
            end
        end
        % If any are not there, then throw an error
        if flag_allFoldersThere==0
            error(['The necessary dependency: %s has an error in install, ' ...
                'or error performing an unzip operation. The subfolders ' ...
                'requested by the code were not found after the unzip ' ...
                'operation. This suggests a bad zip file, or a permissions ' ...
                'error on the local computer, or that folders are ' ...
                'specified that are not present on the remote code ' ...
                'repository.\n'],dependency_name);
        else
            % Clean up the zip file
            delete(zip_file_name);
        end
        
    end
    
    
    % For path creation, if the "DebugTools" package is being installed, the
    % code installs the package, then shifts temporarily into the package to
    % complete the path definitions for MATLAB. If the DebugTools is not
    % already installed, an error is thrown as these tools are needed for the
    % path creation.
    %
    % In other words: DebugTools is a special case because folders not
    % added yet, and we use DebugTools for adding the other directories
    if strcmp(dependency_name(1:10),'DebugTools')
        debugTools_function_folder = fullfile(root_directory_name, 'Utilities', dependency_name,'Functions');
        
        % Move into the folder, run the function, and move back
        cd(debugTools_function_folder);
        fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        cd(root_directory_name);
    else
        try
            fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        catch
            error(['Package installer requires DebugTools package to be ' ...
                'installed first. Please install that before ' ...
                'installing this package']);
        end
    end
    
    
    % Finally, the code sets a global flag to indicate that the folders are
    % initialized.  Check this using a command "exist", which takes a
    % character string (the name inside the '' marks, and a type string -
    % in this case 'var') and checks if a variable ('var') exists in matlab
    % that has the same name as the string. The ~ in front of exist says to
    % do the opposite. So the following command basically means: if the
    % variable named 'flag_CodeX_Folders_Initialized' does NOT exist in the
    % workspace, run the code in the if statement. If we look at the bottom
    % of the if statement, we fill in that variable. That way, the next
    % time the code is run - assuming the if statement ran to the end -
    % this section of code will NOT be run twice.
    
    eval(sprintf('%s = 1;',flag_varname));
end


%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plots
    
    % Nothing to do!
    
    
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function fcn_DebugTools_installDependencies


%% fcn_INTERNAL_calcUnitVector
function unit_vectors = fcn_INTERNAL_calcUnitVector(input_vectors)
vector_length = sum(input_vectors.^2,2).^0.5;
unit_vectors = input_vectors./vector_length;
end