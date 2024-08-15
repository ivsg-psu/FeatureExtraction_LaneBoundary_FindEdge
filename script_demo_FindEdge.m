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
% 2024_08_11 -  Jiabao Zhao, jpz5469@psu.edu
% -- added format string as a optional input for fcn_findEdge_plotLIDARLLA 
% -- functionalize drivable surface 
% 2024_08_12 - Aneesh Batchu
% -- Added more stpes to check if enough data is captured.
% 2024_08_12 -  Jiabao Zhao, jpz5469@psu.edu
% -- added a fcn_findEdge_plotLIDARLLA_Aneesh that is similar to Dr B code 
% but the plot is different, which one we use?
% --functionalize STEP 2.5: Find the LIDAR_ENU and LIDAR_scanLineAndRingID 
% in domain
% --functionalize STEP 2: Find the driven path (left and right side points)
% 2024_08_13 - Jiabao Zhao, jpz5469@psu.edu
% --functionalize %% STEP 4: Find the driven path grids within the grids
% more than zero points.
% --fixed some issues with fig
% 2024_08_13 - Aneesh Batchu
% -- Calculated number of scan lines in each grid without a for loop using
% accumarray
% 2024_08_13 - Aleksandr Goncharov, opg5041@psu.edu
% -- Added Step 6 functions into the main demo.
% 2024_08_15 - Aneesh Batchu
% -- Recalculate driven grid path grids is functionalized



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
%
% script_test for fcn_findEdge_plotLIDARLLA_Aneesh and fcn_findEdge_findPointsInDomain
% script_test for fcn_findEdgefindGridsWithPoints
%
% Write Assertions for all the sections/steps

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

ith_library = ith_library+1;
library_name{ith_library}    = 'PlotRoad_v2024_08_14';
library_folders{ith_library} = {'Functions', 'Data'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/FieldDataCollection_VisualizingFieldData_PlotRoad/archive/refs/tags/PlotRoad_v2024_08_14.zip';


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


%% STEP 1: Load the data
% fcn_findEdge_loadLIDARData:
% [VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num), (fig_num2));
%
% Set the "inputs" to the file loading process - need the date and names
% and date of file creation for the Vehicle Pose data file

fig_num = 1; 
fig_num2 = 2; 
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num),(fig_num2));

%% STEP 2 (part 1): Find the scan lines that are "range of LiDAR" meters away from station 1 and station 2 and find LiDAR ENU and LIDAR_scanLineAndRingID
% fcn_findEdge_pointsAtRangeOfLiDARFromStation:
% [station1_minus_range_index, station2_plus_range_index]= fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,starting_index,ending_index,(range))

station_1 = 1400; 
station_2 = 1450;

range_of_LiDAR = 10;

[station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,station_1,station_2,range_of_LiDAR);  

%% STEP 2 (part 2): find LiDAR ENU and LIDAR_scanLineAndRingID in domain
% fcn_findEdge_pointsAtRangeOfLiDARFromStation:
% [station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,starting_index,ending_index,(range))

scanLineRange = [station1_minus_range_index station2_plus_range_index]; 

Nscans = length(VehiclePose(:,1));
% Set defaults for which scans to extract
% scanLineRange = [1400 1450];
ringsRange = []; % If leave empty, it loads all rings

ENU_XYZ_fig_num = 3;
fig_num = 4;
% Extract scan lines
[VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (ENU_XYZ_fig_num),(fig_num));

%% STEP 2 (part 3): Find the LIDAR_ENU and LIDAR_scanLineAndRingID in domain
%fcn_findEdge_findPointsInDomain:
% [concatenate_LiDAR_XYZ_points_new, boundary_points_of_domain, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2,(fig_num))


fig_num = 7;
[concatenate_LiDAR_XYZ_points_new, boundary_points_of_domain, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2, LIDAR_intensity,(fig_num));

%% Vehicle pose - STEP 1: Load the vehicle pose and the find the driven path.
% fcn_findEdge_findDrivableSurface:
% [LIDAR_ENU_under_vehicle] = fcn_findEdge_findDrivableSurface (LIDAR_ENU, VehiclePose_ENU, VehiclePose_UnitOrthoVectors,(fig_num),(fig_num2))  


ENU_3D_fig_num = 8;
fig_num = 9; % figure number
[LIDAR_ENU_under_vehicle] = fcn_findEdge_findDrivableSurface (LIDAR_ENU, VehiclePose_ENU, VehiclePose_UnitOrthoVectors,(ENU_3D_fig_num),(fig_num));

%% Vehicle pose - STEP 2: Find the driven path (left and right side points)
% fcn_findEdge_findDrivenPathBoundaryPoints:
% fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose, scanLineRange, Nscans, shift, (fig_num))

fig_num = 10;
figure(fig_num);
shift = 5;
fig_num2 = 11;
boundary_points_driven_path = fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose, scanLineRange, Nscans, shift, (fig_num), (fig_num2));

%% STEP 3: Seperate the data into grids
% fcn_findEdge_findMaxMinOfXYZ:
% [Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_findEdge_findMaxMinOfXYZ(N_points,(fig_num))
%
% fcn_findEdge_findGridsWithPoints:
% [gridIndices_cell_array, total_N_points_in_each_grid, gridCenters,
%  grids_with_zero_point, grids_greater_than_zero_points,
%  gridCenters_zero_point_density,gridCenters_greater_than_zero_point_density] = fcn_findEdge_findGridsWithPoints(input_points,grid_size,grid_boundaries,(fig_num))
%
fig_num = 40; 
figure(fig_num);clf

% These are concatenated LiDAR points of chosen scans and cells in the
% first step. 
LiDAR_allPoints = [LIDAR_ENU(in_domain,:), LIDAR_scanLineAndRingID(in_domain,:)];

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

% Seperate the data into grids
[gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, grids_with_zero_points,...
    grids_greater_than_zero_points, gridCenters_zero_point_density,...
    gridCenters_greater_than_zero_point_density, gridIndices, grid_AABBs] = fcn_findEdge_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,fig_num);

%% STEP 4: Find the driven path grids within the grids more than zero points
% fcn_findEdge_findDrivenPathGrids:
%
% [total_points_in_each_grid_in_the_driven_path, total_points_in_each_grid_with_points_greater_than_zero]...
%           = fcn_findEdge_findDrivenPathGrids(gridCenters_greater_than_zero_point_density, boundary_points_driven_path,grids_greater_than_zero_points, (fig_num))

fig_num = 13;
ENU_3D_fig_num = 14;

format = sprintf('''.'',''MarkerSize'',30,''Color'',[0.8 0.8 0.8]');
format1 = sprintf('''.'',''MarkerSize'',30,''Color'',[0 1 0]');


current_grids_greater_than_zero_points = 1:length(grids_greater_than_zero_points); 

[~, ~,total_points_in_each_grid_in_the_driven_path, total_points_in_each_grid_with_points_greater_than_zero]...
    = fcn_findEdge_findDrivenPathGrids(gridCenters_greater_than_zero_point_density, boundary_points_driven_path,...
    grids_greater_than_zero_points, current_grids_greater_than_zero_points, total_N_points_in_each_grid, (format), (format1),[],[], (fig_num), (ENU_3D_fig_num));


%% STEP 5: Grid conditions - Point density (Do not need to run after determining point density)
% fcn_findEdge_determineGridPointDensity:
%
% [point_density] = 
%       fcn_findEdge_determineGridPointDensity(total_points_in_each_grid_with_points_greater_than_zero ...
%       ,total_points_in_each_grid_in_the_driven_path,grid_size,(N_bins_grid_with_points_greater_than_zero)...
%       (N_bins_grid_in_the_driven_path),(fig_num))
%
%fcn_findEdge_classifyGridsBasedOnDensity:
%
% [grid_indices_with_required_point_density, gridCenters_low_point_density] = 
%       fcn_findEdge_determineGridPointDensity(total_points_in_each_grid_with_points_greater_than_zero ...
%       ,total_points_in_each_grid_in_the_driven_path,grid_size,(N_bins_grid_with_points_greater_than_zero)...
%       (N_bins_grid_in_the_driven_path),(fig_num))

fig_num = 52; 
figure(fig_num); clf; 

[point_density] = fcn_findEdge_determineGridPointDensity(total_points_in_each_grid_with_points_greater_than_zero,total_points_in_each_grid_in_the_driven_path,grid_size,[],[],fig_num);

% Minimum number of points required 
point_density = floor(20*((grid_size^2)/(0.3^2)));

% Find the grids that contain enough point density
fig_num = 70; 
figure(fig_num); clf; 

[grid_indices_with_required_point_density, gridCenters_low_point_density] = fcn_findEdge_classifyGridsBasedOnDensity(grids_greater_than_zero_points,total_N_points_in_each_grid,point_density,gridCenters,[],[],fig_num);

%% STEP 5: Grid conditions - Determining number of scan lines in each grid greater than zero points Point density 
% fcn_findEdge_calcNumberOfGridScanLines:
%   [total_scan_lines_in_each_grid_with_more_than_zero_points] = fcn_findEdge_calcNumberOfGridScanLines(gridIndices,LIDAR_scanLines,grids_greater_than_zero_points)
%
% fcn_findEdge_classifyGridsBasedOnScanLines:
%   [grid_indices_with_more_than_one_scan_line, gridCenters_with_one_scan_line] = fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points...
%                                                        ,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters, (format_1),(format_2), (fig_num))
    
[total_scan_lines_in_each_grid_with_more_than_zero_points] = fcn_findEdge_calcNumberOfGridScanLines(gridIndices,LIDAR_scanLines,grids_greater_than_zero_points);


fig_num = 71; 
figure(fig_num); clf;

% Find the grids that contain more than one scan line
[grid_indices_with_more_than_one_scan_line, gridCenters_with_one_scan_line] = fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters,[],[],fig_num);

%% STEP 5: Grid conditions - Finding the orthogonal distances of points in the remaining rings by projecting orthogonally from one ring
%   fcn_findEdge_determineTransverseSpanThreshold:
% [transverse_span_threshold,transverse_span_each_grid]= 
%       fcn_findEdge_determineTransverseSpanThreshold(grids_greater_than_zero_points, grid_AABBs, grid_size, gridIndices, input_points, LIDAR_scanLines, fig_num, fig_num2, fig_num3)
%
%   fcn_findEdge_classifyGridsBasedOnTransverseSpan:
% [grid_indices_with_more_than_transverse_span_threshold, gridCenters_with_less_than_transverse_span_threshold] =
%       fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid,transverse_span_threshold,grids_greater_than_zero_points, gridCenters, (format_1),(format_2),(fig_num))

fig_num_1 = 50; 
figure(fig_num_1); clf; 

fig_num_2 = 5837; 
figure(fig_num_2); clf; 

fig_num_3 = 58309;
figure(fig_num_3); clf; 

[transverse_span_threshold,transverse_span_each_grid] = fcn_findEdge_determineTransverseSpanThreshold...
    (grids_greater_than_zero_points, grid_AABBs, grid_size, gridIndices, input_points, LIDAR_scanLines,...
    gridCenters_greater_than_zero_point_density,fig_num_1,fig_num_2,fig_num_3);


% Threshold of transverse span
transverse_span_threshold = 0.15; 

fig_num = 72; 
figure(fig_num); clf;

% Find the grids
[grid_indices_with_more_than_transverse_span_threshold, gridCenters_with_less_than_transverse_span_threshold] = fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid,transverse_span_threshold,grids_greater_than_zero_points,gridCenters,[],[],fig_num);

%% STEP 6: Qualified and unqualified grids: grids that pass all three conditions above are qualified


fig_num=27;
figure(fig_num); clf;

[current_qualified_grids,current_unqualified_grids,original_qualified_grids,gridCenters_qualified_grids,gridCenters_unqualified_grids] = ...
    fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,grid_indices_with_more_than_transverse_span_threshold,...
    grids_greater_than_zero_points,gridCenters,[],[],fig_num);

%% STEP 6: Plot circles corresponding to each fail condition
% less than required point density - Red (small circle)
% less than one scan line - Green (medium circle)
% less than transverse span threshold - Blue (Large circle)

% Plot qualified and unqualified grids

fig_num = 987; 
figure(fig_num);clf

marker_size = 10;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Unqualified grids';
legend_position = [];
marker_type = [];
plot_gridCenters_unqualified_grids = [gridCenters_unqualified_grids, zeros(length(gridCenters_unqualified_grids(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_unqualified_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot grid centers
marker_size = 10;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Qualified grids';
legend_position = [];
marker_type = [];
plot_gridCenters_qualified_grids = [gridCenters_qualified_grids, zeros(length(gridCenters_qualified_grids(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_qualified_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% Plot the fail conditions
marker_size = 5;
RGB_triplet = [1, 0, 0]; 
legend_option = 1;
legend_name = 'Grids with low point density';
legend_position = [];
marker_type = 'o';
plot_gridCenters_low_point_density = [gridCenters_low_point_density, zeros(length(gridCenters_low_point_density(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_low_point_density,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


marker_size = 7.5;
RGB_triplet = [0, 1, 0]; 
legend_option = 1;
legend_name = 'Grids with one scan line';
legend_position = [];
marker_type = 'o';
plot_gridCenters_with_one_scan_line = [gridCenters_with_one_scan_line, zeros(length(gridCenters_with_one_scan_line(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_with_one_scan_line,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


marker_size = 10;
RGB_triplet = [0, 0, 1]; 
legend_option = 1;
legend_name = 'transverse span threshold fail';
legend_position = [];
marker_type = 'o';
plot_gridCenters_with_less_than_transverse_span_threshold = [gridCenters_with_less_than_transverse_span_threshold, zeros(length(gridCenters_with_less_than_transverse_span_threshold(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_with_less_than_transverse_span_threshold,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 7: Recalculate the driven path grids among qualified grids
% :Find the driven paths in the qualified grids:

% (Function is broken) -- NEED TO FIX IT

% fig_num = 131;
% ENU_3D_fig_num = 114;
% format = sprintf('''.'',''MarkerSize'',40,''Color'',[0.2 0.2 0.2]');
% format1 = sprintf('''o'',''MarkerSize'',10,''Color'',[0 1 0], ''LineWidth'',2');
% legend_name = 'Qualified grids';
% legend_name2 = 'Driven path grids';
% 
% [gridCenters_driven_path, current_grid_numbers_of_driven_path,total_points_in_each_grid_in_the_driven_path, total_points_in_each_grid_with_points_greater_than_zero]...
%     = fcn_findEdge_findDrivenPathGrids(gridCenters_qualified_grids, boundary_points_driven_path,...
%     original_qualified_grids, total_N_points_in_each_grid,(format), (format1),(legend_name),(legend_name2), (fig_num), (ENU_3D_fig_num));
% 
% % [~, ~,total_points_in_each_grid_in_the_driven_path, total_points_in_each_grid_with_points_greater_than_zero]...
% %     = fcn_findEdge_findDrivenPathGrids(gridCenters_greater_than_zero_point_density, boundary_points_driven_path,...
% %     grids_greater_than_zero_points, total_N_points_in_each_grid, (format), (format1),[],[], (fig_num), (ENU_3D_fig_num));
% 

fig_num = 517;
figure(fig_num); clf;

ENU_3D_fig_num = 804; 
figure(ENU_3D_fig_num);clf

[gridCenters_driven_path, current_grid_numbers_of_driven_path,~, ~]...
    = fcn_findEdge_findDrivenPathGrids(gridCenters_qualified_grids, boundary_points_driven_path,...
    original_qualified_grids, current_qualified_grids, total_N_points_in_each_grid, (format), (format1),'Qualified grids',[], (fig_num), (ENU_3D_fig_num));


% % "inpolygon" is used to find the grids within the boundary points 
% [in_qg,on_qg] = inpolygon(gridCenters_qualified_grids(:,1),gridCenters_qualified_grids(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));
% 
% % Original grid numbers of driven path
% original_grid_numbers_of_driven_path = original_qualified_grids(in_qg); 
% 
% % Current grid numbers in driven path 
% current_grid_numbers_of_driven_path = current_qualified_grids(in_qg); %find(in); 
% 
% % Total points in each grid in the driven path
% total_points_in_each_grid_in_the_driven_path = total_N_points_in_each_grid(original_grid_numbers_of_driven_path); 
% 
% % Total points in each grid with points greater than zero
% total_points_in_each_grid_with_points_greater_than_zero = total_N_points_in_each_grid(current_qualified_grids); 
% 
% % Grid centers of the driven path
% gridCenters_driven_path = [gridCenters_qualified_grids(in_qg,1),gridCenters_qualified_grids(in_qg,2)];
% 

% fig_num = 517;
% figure(fig_num); clf;
% 
% hold on
% grid on
% xlabel('X[m]')
% ylabel('Y[m]')
% title('Grid centers and boundary points')
% 
% plot(gridCenters_qualified_grids(:,1), gridCenters_qualified_grids(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);
% % plot(boundary_points_driven_path(:,1), boundary_points_driven_path(:,2), '.', 'MarkerSize',30, 'Color',[0 1 0]); 
% 
% % plot the grids in the driven path
% plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',2) % points strictly inside
% 
% % for ith_text = 1:length(current_qualified_grids(:,1))
% %     current_text = sprintf('%.0d',ith_text);
% %     % Place the text on the grid center
% %     text(gridCenters_qualified_grids(ith_text,1), gridCenters_qualified_grids(ith_text,2),current_text,'Color',[1 1 1],'HorizontalAlignment','center','FontSize', 6, 'FontWeight','bold');
% % end
% 
% 
% % Plot the grids with 
% fig_num = 804; 
% figure(fig_num);clf
% 
% % plot computed boundary points
% marker_size = 10;
% RGB_triplet = [0 0 0]; 
% legend_option = 1;
% legend_name = 'Qualified grids';
% legend_position = [];
% marker_type = [];
% plot_gridCenters_qualified_grids = [gridCenters_qualified_grids(:,1:2), zeros(length(gridCenters_qualified_grids(:,1)),1)];
% [~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_qualified_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
% 
% 
% % plot driven path
% marker_size = 25;
% RGB_triplet = [0 1 0]; 
% legend_option = 1;
% legend_name = 'Driven path grids';
% legend_position = [];
% marker_type = [];
% plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
% [~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
% 

%% STEP 8: Qualified grid conditions - Standard deviation in Z (Do not need to run after determining a std_threshold)
% :Fit a plane to each grid using linear regression:
% (NEED TO FIX THIS) - STD THRESHOLD - OUTPUT

fig_num_1 = 81;
figure(fig_num_1); clf;

fig_num_2 = 82;
figure(fig_num_2); clf;

fig_num_3 = 83;
figure(fig_num_3); clf;

[~,original_mapped_gridIndices_cell,total_mapped_grids,total_points_in_mapped_grids,standard_deviation_in_z,gridlines_mapped_grids,...
    driven_path_grid_indices_in_current_mapped_grids,std_in_z_driven_path,std_in_z_other_mapped_grids,mean_std_in_z_driven_path,mean_std_in_z_not_driven_path,max_std_in_z_not_driven_path] = fcn_findEdge_determineSTDInZError(LiDAR_allPoints,gridIndices_cell_array,original_qualified_grids,gridCenters_qualified_grids,gridCenters_driven_path,...
    current_qualified_grids,grid_AABBs,grid_size,gridIndices,current_grid_numbers_of_driven_path,fig_num_1,fig_num_2,fig_num_3);

% Determined std_threshold
std_threshold = 0.1; 

%% STEP 9: Histogram of standard deviation - (Do not need to run after determining a std_threshold)
% :Find the standard deviations of the z-fit errors for all qualified grids, and compare these standard deviations with those of the grids in the driven path to determine the suitable standard deviation threshold:
 
fig_num = 123;
figure(fig_num); clf;
N_bins_stdz=100;
N_bins_std_drivenpath=5;

fcn_findEdge_histogramSTDinZError(standard_deviation_in_z,N_bins_stdz,std_in_z_driven_path,N_bins_std_drivenpath,mean_std_in_z_driven_path,std_threshold,(fig_num));

%% STEP 8: Qualified grid conditions - angle deviation (Do not need to run after determining a theta_threshold)
% :Fit a plane to each grid using linear regression:

fig_num_1 = 883;
figure(fig_num_1); clf;

fig_num_2 = 84;
figure(fig_num_2); clf;

fig_num_3 = 85;
figure(fig_num_3); clf;

[angle_btw_unit_normals_and_vertical, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path,...
    angle_btw_unit_normals_and_vertical_driven_path] = fcn_findEdge_determineAngleDeviation...
    (LiDAR_allPoints, gridIndices_cell_array, original_qualified_grids,...
    gridCenters_qualified_grids,current_qualified_grids,gridCenters_driven_path, ...
    grid_AABBs, grid_size, gridIndices, current_grid_numbers_of_driven_path, fig_num_1,fig_num_2,fig_num_3);

%% STEP 9: Histogram of angle deviation - (Do not need to run after determining a theta_threshold)
% :Find the angles between the unit normal vectors of the fitted planes and vertical ([0 0 1]) for all qualified grids, and compare these angles with those of the grids in the driven path to determine the suitable theta threshold:

fig_num = 1223;
figure(fig_num); clf;

theta_threshold = fcn_findEdge_histogramAngleDeviation(angle_btw_unit_normals_and_vertical, ...
    angle_btw_unit_normals_and_vertical_driven_path, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path, fig_num);

%% STEP 10: Voting - Drivable, Non-drivable and Uncertain
% :Voting:

% NOTE: The code refers non-drivable grids as failed grids and non-drivable
% as (failed grids and uncertain grids) - NEED to change

input_points = LiDAR_allPoints(:,1:3); 

% Give the thresholds if already found
% std_threshold = 0.05; 
% theta_threshold = 7*pi/180;
% theta_threshold = 30*pi/180;
% gridCenters

fig_num = 6000; 
figure(fig_num);clf

% theta_threshold = [];
% std_threshold = [];


% Classify mapped grids into drivable and drivable
[standard_deviation_in_z, ... 
    angle_btw_unit_normals_and_vertical, ...
    original_drivable_grids, ...
    original_non_drivable_grids, ...
    current_drivable_grid_numbers_in_mapped_grids, ... 
    current_non_drivable_grid_numbers_in_mapped_grids, ...
    current_failed_grid_numbers_in_mapped_grids, ... 
    current_uncertain_grid_numbers_in_mapped_grids, ... 
    gridCenters_failed_grids, ... 
    gridCenters_uncertain_grids, ...
    gridCenters_drivable_grids, ...
    gridCenters_non_drivable_grids, ... 
    concatenate_gridCenters_drivable_non_drivable_grids] = ...
    fcn_findEdge_classifyGridsAsDrivable(gridIndices_cell_array, original_qualified_grids, input_points, std_threshold, theta_threshold, gridCenters, fig_num);

%% STEP 10: Plot the drivable, non-drivable and uncertain grid centers 
% :Voting:

% (The code refers non-drivable as failed) Non-drivable: Everything (uncertain and failed) other than drivable
 
% Plot qualified and unqualified grids

fig_num = 987; 
figure(fig_num);clf

marker_size = 25;
RGB_triplet = [0 1 0]; 
legend_option = 1;
legend_name = 'Drivable grids';
legend_position = [];
marker_type = [];
plot_gridCenters_drivable_grids = [gridCenters_drivable_grids(:,1:2), zeros(length(gridCenters_drivable_grids(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_drivable_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot grid centers
marker_size = 25;
RGB_triplet = [1 0 0]; 
legend_option = 1;
legend_name = 'Non-drivable grids';
legend_position = [];
marker_type = [];
plot_gridCenters_failed_grids = [gridCenters_failed_grids(:,1:2), zeros(length(gridCenters_failed_grids(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_failed_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot grid centers
marker_size = 25;
RGB_triplet = [0 0 1];%[0.9290 0.6940 0.1250]; 
legend_option = 1;
legend_name = 'Uncertain grids';
legend_position = [];
marker_type = [];
plot_gridCenters_uncertain_grids = [gridCenters_uncertain_grids(:,1:2), zeros(length(gridCenters_uncertain_grids(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_uncertain_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

if ~isempty(std_threshold)
    % Find the indices of std threshold failed indices
    current_std_threshold_failed_indices = current_uncertain_grid_numbers_in_mapped_grids(standard_deviation_in_z(current_uncertain_grid_numbers_in_mapped_grids)>std_threshold);

    % Find the grid centers of the std_threshold failed grid centers
    std_threshold_failed_gridCenters = gridCenters_qualified_grids(current_std_threshold_failed_indices,1:2);
end

marker_size = 10;
RGB_triplet = [1 1 0]; 
legend_option = 1;
legend_name = 'Std threshold failed uncertain grids';
marker_type = 'o';
legend_position = [];
plot_std_threshold_failed_gridCenters = [std_threshold_failed_gridCenters, zeros(length(std_threshold_failed_gridCenters(:,1)),1)]; 
[LLA_data_std_threshold_failed_grids] = fcn_findEdge_plotPointsinLLA(plot_std_threshold_failed_gridCenters,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],-1);
% fcn_geometry_plotPointsinLLA(plot_gridCenters_updated_original_unmapped_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% fcn_geometry_plotPointsinLLA(plot_gridCenters_non_drivable_grids,marker_size,RGB_triplet,[],legend_option,legend_name,[],[],[],[],fig_num);

geoplot(LLA_data_std_threshold_failed_grids(:,1), LLA_data_std_threshold_failed_grids(:,2), 'o','MarkerSize',10,'Color',RGB_triplet, 'LineWidth',2, 'DisplayName','Std threshold failed grids') 


if ~isempty(theta_threshold)
    % Find the indices of std threshold failed indices
    current_theta_threshold_failed_indices = current_uncertain_grid_numbers_in_mapped_grids(angle_btw_unit_normals_and_vertical(current_uncertain_grid_numbers_in_mapped_grids)>theta_threshold);

    % Find the grid centers of the theta_threshold failed grid centers
    theta_threshold_failed_gridCenters = gridCenters_qualified_grids(current_theta_threshold_failed_indices,1:2);
end


marker_size = 12.5;
RGB_triplet = [0 1 1]; 
legend_option = 1;
legend_name = 'Theta threshold failed uncertain grids';
marker_type = 'o';
legend_position = [];
plot_theta_threshold_failed_gridCenters = [theta_threshold_failed_gridCenters, zeros(length(theta_threshold_failed_gridCenters(:,1)),1)]; 
[LLA_data_theta_threshold_failed_grids] = fcn_findEdge_plotPointsinLLA(plot_theta_threshold_failed_gridCenters,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],-1);

geoplot(LLA_data_theta_threshold_failed_grids(:,1), LLA_data_theta_threshold_failed_grids(:,2), 'o','MarkerSize',marker_size,'Color',RGB_triplet, 'LineWidth',2, 'DisplayName','Theta threshold failed grids') 



%% STEP 11: Find all boundary points: When uncertain grids are assumed as drivable
% :Boundary points:

% Part1 - Find the boundary points of mapped and unmapped grids

% Revision History
% Funtionalized this code
% Added plotting options

% INPUTS - gridCenters_low_point_density,
% gridCenters_required_point_density, figure num
% OUTPUTS - X, Y, Z 

% prepGridCentersForBoundaryDetection


fig_num_qualified_unqualified = 767787; 

[X, Y, Z] = fcn_findEdge_prepGridCentersForBoundaryDetection...
    (gridCenters_qualified_grids, gridCenters_unqualified_grids);

%%%%%%%%%%%%%%---------------------------------------------------------------------------
% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
figure(fig_num_qualified_unqualified);
clf;
boundary_points_mapped_unmapped = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num_qualified_unqualified);

% assert(length(drivable_grids)>=1)
% assert(length(non_drivable_grids)>=1)
% % assert(length(unmapped_grids)==zeros(0,1))
% assert(length(gridCenters_mapped_grids(:,1))>=1)
% assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(boundary_points_mapped_unmapped)>=1)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundary points of drivable and non-drivable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 2 - Find boundary points of drivable and non-drivable grids

fig_num_drivable_non_drivable = 98898;

[X, Y, Z] = fcn_findEdge_prepGridCentersForBoundaryDetection...
    (gridCenters_drivable_grids, gridCenters_uncertain_grids);

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

% Part 3 - Find true boundary points by eliminating the boundary points of
% mapped and unmapped from drivable and non-drivable boubndary points. 

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

fig_num = 452;
figure(fig_num); clf; 
% plot computed boundary points
marker_size = 25;
RGB_triplet = [0 1 1]; 
legend_option = 1;
legend_name = 'Computed boundary points';
legend_position = [];
marker_type = [];
plot_true_boundary_points = [true_boundary_points, zeros(length(true_boundary_points),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot computed boundary points

marker_size = 25;
RGB_triplet = [0 1 0]; 
legend_option = 1;
legend_name = 'Drivable grids';
legend_position = [];
marker_type = [];
plot_gridCenters_drivable_grids = [gridCenters_drivable_grids(:,1:2), zeros(length(gridCenters_drivable_grids(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_drivable_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot grid centers
marker_size = 25;
RGB_triplet = [1 0 0]; 
legend_option = 1;
legend_name = 'Non-drivable grids';
legend_position = [];
marker_type = [];
plot_gridCenters_failed_grids = [gridCenters_failed_grids(:,1:2), zeros(length(gridCenters_failed_grids(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_failed_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot grid centers
marker_size = 25;
RGB_triplet = [0 0 1];%[0.9290 0.6940 0.1250]; 
legend_option = 1;
legend_name = 'Uncertain grids';
legend_position = [];
marker_type = [];
plot_gridCenters_uncertain_grids = [gridCenters_uncertain_grids(:,1:2), zeros(length(gridCenters_uncertain_grids(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_uncertain_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% marker_size = 10;
% RGB_triplet = [0 0 1]; 
% legend_option = 0;
% legend_name = 'Computed boundary points';
% legend_position = [];
% marker_type = [];
% [~] = fcn_geometry_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot driven path
marker_size = 25;
RGB_triplet = [0 0.5 0]; 
legend_option = 0;
legend_name = 'Driven path grids';
legend_position = [];
marker_type = [];
plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot driven path
marker_size = 10;
RGB_triplet = [0 0.8 0]; 
legend_option = 0;
legend_name = 'Driven path grids';
legend_position = [];
marker_type = [];
plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 11: Find boundary points: When uncertain grids are assumed as non drivable
% :Boundary points:

% Part1 - Find the boundary points of mapped and unmapped grids

% Revision History
% Funtionalized this code
% Added plotting options

% INPUTS - gridCenters_low_point_density,
% gridCenters_required_point_density, figure num
% OUTPUTS - X, Y, Z 

fig_num_qualified_unqualified = 767787; 

[X, Y, Z] = fcn_findEdge_prepGridCentersForBoundaryDetection...
    (gridCenters_qualified_grids, gridCenters_unqualified_grids);

%%%%%%%%%%%%%%---------------------------------------------------------------------------
% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
figure(fig_num_qualified_unqualified);
clf;
boundary_points_mapped_unmapped = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num_qualified_unqualified);

% assert(length(drivable_grids)>=1)
% assert(length(non_drivable_grids)>=1)
% % assert(length(unmapped_grids)==zeros(0,1))
% assert(length(gridCenters_mapped_grids(:,1))>=1)
% assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(boundary_points_mapped_unmapped)>=1)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundary points of drivable and non-drivable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 2 - Find boundary points of drivable and non-drivable grids

fig_num_drivable_non_drivable = 98898;


[X, Y, Z] = fcn_findEdge_prepGridCentersForBoundaryDetection...
    (gridCenters(original_drivable_grids,1:2), gridCenters(original_non_drivable_grids,1:2));

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

% Part 3 - Find true boundary points by eliminating the boundary points of
% mapped and unmapped from drivable and non-drivable boubndary points. 

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

fig_num = 434;
figure(fig_num); clf;

marker_size = 25;
RGB_triplet = [0 1 0]; 
legend_option = 1;
legend_name = 'Drivable grids';
legend_position = [];
marker_type = [];
plot_gridCenters_drivable_grids = [gridCenters_drivable_grids(:,1:2), zeros(length(gridCenters_drivable_grids(:,1)),1)]; 
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_drivable_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot grid centers
marker_size = 25;
RGB_triplet = [1 0 0]; 
legend_option = 1;
legend_name = 'NOT drivable grids';
legend_position = [];
marker_type = [];
plot_gridCenters_non_drivable_grids = [gridCenters_non_drivable_grids(:,1:2), zeros(length(gridCenters_non_drivable_grids(:,1)),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_non_drivable_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot computed boundary points
marker_size = 25;
RGB_triplet = [0 1 1]; 
legend_option = 1;
legend_name = 'Computed boundary points';
legend_position = [];
marker_type = [];
plot_true_boundary_points = [true_boundary_points, zeros(length(true_boundary_points),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% % plot computed boundary points
% marker_size = 10;
% RGB_triplet = [0 0 1]; 
% legend_option = 0;
% legend_name = 'Computed boundary points';
% legend_position = [];
% marker_type = [];
% [~] = fcn_geometry_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
% 

% plot driven path
marker_size = 25;
RGB_triplet = [0 0.5 0]; 
legend_option = 0;
legend_name = 'Driven path grids';
legend_position = [];
marker_type = [];
plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot driven path
marker_size = 10;
RGB_triplet = [0 0.8 0]; 
legend_option = 0;
legend_name = 'Driven path grids';
legend_position = [];
marker_type = [];
plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% plot the boundary points
figure(7639);clf;

% plot computed boundary points
marker_size = 25;
RGB_triplet = [0 1 1]; 
legend_option = 1;
legend_name = 'Computed boundary points';
legend_position = [];
marker_type = [];
plot_true_boundary_points = [true_boundary_points, zeros(length(true_boundary_points),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 12: Find the nearest boundary points
% :The nearest boundary points:

fig_num = 7676;

% Find the nearest boundaries
[~, nearestBorderIndicies, nearestBorderXY] = fcn_findEdge_findNearestBoundaryPoints(true_boundary_points, ...
    gridCenters_driven_path, grid_size, grid_boundaries, fig_num);

% Figure number
fig_num = 1224; 
figure(fig_num); clf;

% plot the nearest boundary points
marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 0;
legend_name = 'Nearest Boundary Points';
legend_position = [];
marker_type = []; 
nearest_boundary_points = [nearestBorderXY(:,1), nearestBorderXY(:,2), zeros(length(nearestBorderXY(:,1)),1)]; 

[~] = fcn_findEdge_plotPointsinLLA(nearest_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

marker_size = 10;
RGB_triplet = [1 1 1]; 
legend_option = 0;
legend_name = 'Nearest Boundary Points';
legend_position = [];
marker_type = []; 
nearest_boundary_points = [nearestBorderXY(:,1), nearestBorderXY(:,2), zeros(length(nearestBorderXY(:,1)),1)]; 

[~] = fcn_findEdge_plotPointsinLLA(nearest_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 13: Seperate the right and left boundaries from the nearest boundaries 
% :Right nearest boundary points: & :Left nearest boundary points:

% Transverse shift 
transverse_shift = 6*3.6576; 
fig_num = 7876; 
[boundary_points_left, boundary_points_right] = fcn_findEdge_seperateLeftRightBoundaries...
    (VehiclePose,station_1,station_2,nearest_boundary_points, grid_size, transverse_shift, fig_num);

%% Plot right boundary points

fig_num = 12824; 
figure(fig_num); 

marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 0;
legend_name = 'Right Boundary Points';
legend_position = [];
marker_type = []; 

plot_boundary_points_right = [boundary_points_right, zeros(length(boundary_points_right),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

marker_size = 12;
RGB_triplet = [1 0 0]; 
legend_option = 0;
legend_name = 'Right Boundary Points';
legend_position = [];
marker_type = []; 

plot_boundary_points_right = [boundary_points_right, zeros(length(boundary_points_right),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% Plot right boundary points

fig_num = 1241; 
figure(fig_num); 

marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 0;
legend_name = 'Left Boundary Points';
legend_position = [];
marker_type = []; 

plot_boundary_points_left = [boundary_points_left, zeros(length(boundary_points_left),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_boundary_points_left,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

marker_size = 12;
RGB_triplet = [0 0 1]; 
legend_option = 0;
legend_name = 'Left Boundary Points';
legend_position = [];
marker_type = []; 

plot_boundary_points_left = [boundary_points_left, zeros(length(boundary_points_left),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_boundary_points_left,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

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