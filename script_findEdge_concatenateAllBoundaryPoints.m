
% script_findEdge_concatenateAllBoundaryPoints
%
% This script finds the boundary points by dividing the scan lines into
% small segments and concatenates all the boundary points at the end

% Revision History
%
% 2024_08_28 - Aneesh Batchu
% -- wrote the code originally
% 2024_10_07 - Aneesh Batchu
% -- Fixed scanLineRange error
%%
clc
clear
close all

%%
savefile = fullfile(pwd, 'Data', 'point6_ComputedBoundaryPoints_Sample1.mat');
load(savefile, 'hand_labeled_path');


%% STEP 1: Load the LiDAR and Vehicle Pose data
fig_num = 101; 
% fig_num2 = 102; 

fig_num2 = -1; 
% test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
% vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
% LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data

test_date_string = '2024_08_05'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
vehicle_pose_string = 'VehiclePose_Seg3-2-1_CCW_Run6.mat'; % The name of the file containing VehiclePose
LIDAR_file_string   = 'VelodyneLiDARScan_In_ENU_Seg3-2-1_CCW_Run6.mat'; % The name of the file containing the LIDAR data


flag_load_all_data = [];
[VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num),(fig_num2));

% % Check sizes
% assert(isequal(length(VehiclePose(:,1)),length(LiDAR_Scan_ENU_Entire_Loop)));
% assert(isequal(length(VehiclePose(1,:)),6)); % XYZRPY (roll pitch yaw)
% assert(isequal(length(LiDAR_Scan_ENU_Entire_Loop{1}(1,:)),6)); % XYZ intensity scanLine deltaT
%%
% VehiclePose = VehicleOutput(1:length(LiDAR_Scan_Transformed_cell),1:6); 
% LiDAR_Scan_ENU_Entire_Loop = LiDAR_Scan_Transformed_cell; 
% 
% % Check for empty rows
% emptyRows = cellfun(@isempty, LiDAR_Scan_ENU_Entire_Loop);
% 
% 
% VehiclePose = VehiclePose(~emptyRows,1:6); 
% LiDAR_Scan_ENU_Entire_Loop = LiDAR_Scan_ENU_Entire_Loop(~emptyRows,:);

%% Find the scan line ranges 

% Intialize: Empty matrix
boundary_points_test_track_right = []; 
boundary_points_test_track_left = []; 

% Choose a grid size
grid_size = 1.25; %0.5; %0.8;%1;%1.25; %1.26
theta_threshold = 0.1745; %0.1745;%0.25; 
% Enter the start scan line
start_scan_line = 1;
% Enter the end scan line
end_scan_line = length(LiDAR_Scan_ENU_Entire_Loop); % Total scan lines
scanLine_segment_length = 50; % Divide the scan lines based on this interval
tic
% Sets range of LiDAR to zero if you are processing the whole test track
% data. 
if end_scan_line == length(LiDAR_Scan_ENU_Entire_Loop) || start_scan_line == 1
    range_of_LiDAR = 0;
else
    range_of_LiDAR = 100;
end

% Create the start points
start_scanLines = start_scan_line:scanLine_segment_length:end_scan_line;

% Create the end points
end_scanLines = start_scanLines + scanLine_segment_length - 1;
end_scanLines(end) = end_scan_line; % Adjust the last end point to match end_scan_line

% Combine the start and end points into a matrix
scanLineRanges = [start_scanLines(:), end_scanLines(:)];

% Delete the range if start and end range are same
scanLineRanges = scanLineRanges(~(scanLineRanges(:,1) == scanLineRanges(:,2)),:);

for scanLineRange_id = 1:length(scanLineRanges)

    % STEP 2: Find the scan lines that are "range of LiDAR" meters away from station 1 and station 2
    
    fig_num = -1;
    scanLineStart = scanLineRanges(scanLineRange_id,1);
    scanLineEnd = scanLineRanges(scanLineRange_id,2);

    % Sets range of LiDAR to zero if you are processing the whole test track
    % data.
    if scanLineEnd == length(LiDAR_Scan_ENU_Entire_Loop) || scanLineStart == 1
        range_of_LiDAR = 0;
    else
        range_of_LiDAR = 100;
    end


    disp(['The scan line range of iteration ', num2str(scanLineRange_id), ' is [', num2str(scanLineStart),' ',num2str(scanLineEnd), '].']);

    h_waitbar = waitbar(0,'Grid preparation ...');
   

    [scanLineStart_minus_range_index, scanLineEnd_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,scanLineStart,scanLineEnd,range_of_LiDAR,fig_num);

    if isempty(scanLineStart_minus_range_index)
        scanLineStart_minus_range_index = scanLineStart;
    elseif isempty(scanLineEnd_plus_range_index)
        scanLineEnd_plus_range_index = scanLineEnd;
    end
    
    %STEP 3: Extracts vehicle pose ENU, vehicle pose unit orthogonal vectors, LiDAR ENU scans, LiDAR Intensity, LiDAR scan line and Ring IDs of the scan line range
    scanLineRange = [scanLineStart_minus_range_index scanLineEnd_plus_range_index];

    ringsRange = []; % If leave empty, it loads all rings

    ENU_XYZ_fig_num = -1;
    fig_num = -1;

    % Extract scan lines
    [VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
        LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
        fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (ENU_XYZ_fig_num),(fig_num));

    % STEP 3.1: Find the boundary points of the driven path to create a bounding box for finding the driven path grids
    fig_num = -1;
    fig_num2 = -1;

    Nscans = length(VehiclePose(:,1));
    shift = 5;

    boundary_points_driven_path = fcn_findEdge_findDrivenPathBoundaryPoints(VehiclePose, scanLineRange, Nscans, shift, (fig_num), (fig_num2));

    % STEP 4: Find the points in the domain from LiDAR ENU scans of the scan line range
    fig_num = -1;
    [concatenate_LiDAR_XYZ_points_new, boundary_points_of_domain, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, scanLineStart, scanLineEnd, LIDAR_intensity,(fig_num));

    % STEP 5: Find the grid boundaries and separate the data into grids to find the empty and non-empty grids
    fig_num = -1;

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
    % grid_size = 1; %0.5; %0.8;%1;%1.25; %1.26

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

    % STEP 6: Find the driven path grids from the non-empty grids
    fig_num = -1;
    ENU_3D_fig_num = -1;

    format = [];
    format1 = [];

    current_grids_greater_than_zero_points = 1:length(grids_greater_than_zero_points);

    [~, ~,total_points_in_each_grid_in_the_driven_path, total_points_in_each_grid_with_points_greater_than_zero]...
        = fcn_findEdge_findDrivenPathGrids(gridCenters_greater_than_zero_point_density, boundary_points_driven_path,...
        grids_greater_than_zero_points, current_grids_greater_than_zero_points, total_N_points_in_each_grid, (format), (format1),[],[], (fig_num), (ENU_3D_fig_num));
    
    waitbar(0.25,h_waitbar,'Drivability of grids ...');

    fig_num = -1;
   
    % STEP 7.1: Grid conditions - Point density
    % Minimum number of points required
    point_density = floor(20*((grid_size^2)/(0.3^2)));

    [grid_indices_with_required_point_density, ...
        current_grids_with_low_point_density, ...
        current_grids_with_required_point_density, ...
        gridCenters_low_point_density, ...
        gridCenters_required_point_density]  = fcn_findEdge_classifyGridsBasedOnDensity(grids_greater_than_zero_points,total_N_points_in_each_grid,point_density,gridCenters,[],[],fig_num);

    % STEP 7.2: Grid conditions - Determining number of scan lines in each grid greater than zero points Point density
    [total_scan_lines_in_each_grid_with_more_than_zero_points] = fcn_findEdge_calcNumberOfGridScanLines(gridIndices,LIDAR_scanLines,grids_greater_than_zero_points);

    fig_num = -1;

    % Find the grids that contain more than one scan line
    [grid_indices_with_more_than_one_scan_line, ...
        current_grids_with_more_than_one_scan_line, ...
        current_grids_with_one_scan_line, ...
        gridCenters_with_more_than_one_scan_line, ...
        gridCenters_with_one_scan_line] = fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters,[],[],fig_num);

    % STEP 7.3: Grid conditions - Finding the orthogonal distances of points in the remaining rings by projecting orthogonally from one ring
    fig_num_1 = -1;

    fig_num_2 = -1;

    fig_num_3 = -1;

    transverse_span_each_grid = fcn_findEdge_determineTransverseSpanThreshold...
        (grids_greater_than_zero_points, grid_AABBs, grid_size, gridIndices, input_points, LIDAR_scanLines,...
        gridCenters_greater_than_zero_point_density,fig_num_1,fig_num_2,fig_num_3);

    % Threshold of transverse span
    transverse_span_threshold = 0.15;

    % Find the grids
    [grid_indices_with_more_than_transverse_span_threshold, gridCenters_with_less_than_transverse_span_threshold] = fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid,transverse_span_threshold,grids_greater_than_zero_points,gridCenters,[],[],fig_num);

    % STEP 8: Qualified and unqualified grids: grids that pass all three conditions above are qualified
    fig_num = -1;

    [current_qualified_grids, ...
        current_unqualified_grids, ...
        original_qualified_grids, ...
        original_unqualified_grids, ...
        gridCenters_qualified_grids, ...
        gridCenters_unqualified_grids] = fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,grid_indices_with_more_than_transverse_span_threshold,...
        grids_greater_than_zero_points,gridCenters,[],[],fig_num);

    % STEP 9: Recalculate the driven path grids among qualified grids
    fig_num = -1;

    ENU_3D_fig_num = -1;

    [gridCenters_driven_path, current_grid_numbers_of_driven_path,~, ~] = fcn_findEdge_findDrivenPathGrids(gridCenters_qualified_grids, boundary_points_driven_path,...
        original_qualified_grids, current_qualified_grids, total_N_points_in_each_grid, (format), (format1),'Qualified grids',[], (fig_num), (ENU_3D_fig_num));
    
    waitbar(0.50,h_waitbar,'Grid voting ...');

    % STEP 11: Voting - Drivable, Non-drivable and Uncertain

    input_points = LiDAR_allPoints(:,1:3);

    fig_num = -1;

    % theta_threshold = 0.1745;%0.25;%0.1745; % (9.98/180)*pi
    std_threshold = 0.1;

    % theta_threshold = 0.15; % (9.98/180)*pi
    % std_threshold = 0.08;


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

waitbar(0.75,h_waitbar,'Post processing ...');

% STEP 12.1: Prepare the grid centers of qualified and unqualified for boundary detection
[Xcoord_gridCenters, Ycoord_gridCenters, Zcoord_gridCenters] = fcn_findEdge_prepGridCentersForBoundaryDetection...
    (gridCenters_qualified_grids, gridCenters_unqualified_grids);

% STEP 13.1: Find the boundary points of qualified and unqualified grids
x_limits = [];  
y_limits = []; 

fig_num_qualified_unqualified = -1; 

boundary_points_qualified_unqualified = fcn_findEdge_findBoundaryPoints(Xcoord_gridCenters,Ycoord_gridCenters,Zcoord_gridCenters,grid_size,x_limits,y_limits,fig_num_qualified_unqualified);

% STEP 12.2: Prepare the grid centers of drivable and non-drivable for boundary detection
[Xcoord_gridCenters, Ycoord_gridCenters, Zcoord_gridCenters] = fcn_findEdge_prepGridCentersForBoundaryDetection...
    (gridCenters_drivable_grids, gridCenters_uncertain_grids);

% STEP 13.2: Find the boundary points of drivable and non-drivable grids
x_limits = [];  
y_limits = []; 

fig_num_drivable_non_drivable = -1; 
boundary_points = fcn_findEdge_findBoundaryPoints(Xcoord_gridCenters,Ycoord_gridCenters,Zcoord_gridCenters,grid_size,x_limits,y_limits,fig_num_drivable_non_drivable);

% STEP 14: Find the true boundary points by removing the boundary points of qualified and unqualified from drivable and non-drivable boundary points
fig_num_bd_pts_ENU = -1; 

[true_boundary_points] = fcn_findEdge_findTrueBoundaryPoints(boundary_points,boundary_points_qualified_unqualified,fig_num_bd_pts_ENU);

% STEP 15: Find the nearest boundary points to the driven path 
fig_num = -1; 

% Find the nearest boundaries
[~, nearestBorderIndicies, nearestBorderXY] = fcn_findEdge_findNearestBoundaryPoints(true_boundary_points, ...
    gridCenters_driven_path, grid_size, grid_boundaries, fig_num);

if ~isempty(nearestBorderXY)
    nearest_boundary_points = [nearestBorderXY(:,1), nearestBorderXY(:,2), zeros(length(nearestBorderXY(:,1)),1)];

    % STEP 16: Seperate the right and left boundaries from the nearest boundaries
    fig_num = -1;

    % Transverse shift
    % transverse_shift = 6*3.6576;

    % [boundary_points_left, boundary_points_right] = fcn_findEdge_seperateLeftRightBoundaries...
    %     (VehiclePose, scanLineStart, scanLineEnd, nearest_boundary_points, grid_size, transverse_shift, fig_num);

    [boundary_points_left, boundary_points_right] = fcn_findEdge_seperateLeftRightBoundaries...
        (VehiclePose, scanLineStart, scanLineEnd, nearest_boundary_points, fig_num);


    boundary_points_test_track_right = [boundary_points_test_track_right; boundary_points_right]; %#ok<AGROW>
    boundary_points_test_track_left = [boundary_points_test_track_left; boundary_points_left]; %#ok<AGROW>
end


plot_boundary_pts_of_segments = 0; 

if(plot_boundary_pts_of_segments == 1)
    fig_num = figure(scanLineRange_id);
    figure(fig_num);

    marker_size = 30;
    RGB_triplet = [0 1 1];
    legend_option = 0;
    legend_name = 'Right Boundary Points';
    legend_position = [];
    marker_type = [];

    plot_boundary_points_right = [boundary_points_right, zeros(length(boundary_points_right),1)];
    [~] = fcn_findEdge_plotPointsinLLA(plot_boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

    if (fig_num>0)
        temp_h = figure(fig_num);
        assert(~isempty(get(temp_h,'Children')))
    end


    marker_size = 12;
    RGB_triplet = [1 0 0];
    legend_option = 0;
    legend_name = 'Right Boundary Points';
    legend_position = [];
    marker_type = [];

    plot_boundary_points_right = [boundary_points_right, zeros(length(boundary_points_right),1)];
    [~] = fcn_findEdge_plotPointsinLLA(plot_boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
end
waitbar(1,h_waitbar,'Finishing the current iteration ...');
close(h_waitbar)

end

% save('boundary_points_test_track.mat', 'boundary_points_test_track'); 
toc
%%

% % Specify the full path to the target location
% filePath = '/Users/aneeshbatchu/Documents/IVSG/FeatureExtraction_LaneBoundary_FindEdge/Data/point5_ComputedBoundaryPoints_Sample1.mat';
% 
% % Save the file
% save(filePath, 'boundary_points_test_track_right');

savefile = fullfile(pwd, 'Data', 'onePoint25_ComputedBoundaryPoints_Sample7.mat');
save(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')

%% Plot the boundary points

fig_num = 408;
figure(fig_num); clf; 

marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 0;
legend_name = 'Right Boundary Points';
legend_position = [];
marker_type = []; 

plot_boundary_points_right = [boundary_points_test_track_right, zeros(length(boundary_points_test_track_right),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


marker_size = 12;
RGB_triplet = [1 0 0]; 
legend_option = 0;
legend_name = 'Right Boundary Points';
legend_position = [];
marker_type = []; 

plot_boundary_points_right = [boundary_points_test_track_right, zeros(length(boundary_points_test_track_right),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

saveas(figure(fig_num), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_boundaryPointPlots/onePoint25_Sample7.jpg')
saveas(figure(fig_num), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_boundaryPointPlots/MATLAB_figs/onePoint25_Sample7.fig')


% %% Plot the boundary points
% if ~isempty(boundary_points_test_track_left)
%     fig_num = 409;
%     figure(fig_num); clf;
% 
%     marker_size = 30;
%     RGB_triplet = [0 1 1];
%     legend_option = 0;
%     legend_name = 'Left Boundary Points';
%     legend_position = [];
%     marker_type = [];
% 
%     plot_boundary_points_left = [boundary_points_test_track_left, zeros(size(boundary_points_test_track_left,1),1)];
%     [~] = fcn_findEdge_plotPointsinLLA(plot_boundary_points_left,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
% 
% 
%     marker_size = 12;
%     RGB_triplet = [0 0 1];
%     legend_option = 0;
%     legend_name = 'Left Boundary Points';
%     legend_position = [];
%     marker_type = [];
% 
%     plot_boundary_points_left = [boundary_points_test_track_left, zeros(size(boundary_points_test_track_left,1),1)];
%     [~] = fcn_findEdge_plotPointsinLLA(plot_boundary_points_left,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
% end