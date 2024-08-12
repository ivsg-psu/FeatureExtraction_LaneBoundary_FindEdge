%% script_test_geometry_findPointsInDomain

% INPUTS
% station_1, station_2, concatenated matrices, vehicle pose
% 
% Run this after running 
% "script_test_plot_sample_data"

% Revision History
% 
% 2024_08_02 - Aneesh Batchu
% -- wrote the code originally

vehicle_change_in_pose_XY = diff(VehiclePose(:,1:2));

% Repeat the last value again, since diff removes one row. We want the same
% number of vectors as the number of points, and diff removed one point.
vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];

% Convert these to unit vectors
unit_vehicle_change_in_pose_XY = fcn_geometry_calcUnitVector(vehicle_change_in_pose_XY);

% Find orthogonal vetors by rotating by 90 degrees in the CCW direction
unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY*[0 1; -1 0];

% Define lane width limits. Take 40 percent of the lane width. Numbers
% obtained by converting 12 ft (standard lane) to meters
% lane_half_width = (3.6576/2) * 0.40; 
% Right transverse shift 
right_transverse_shift = 6*3.6576;  

% Transverse distance of the right boundary points from vehicle center 
right_transverse_distance_of_boundary_points = [right_transverse_shift*unit_ortho_vehicle_vectors_XY, zeros(length(unit_ortho_vehicle_vectors_XY),1)];

% Left transverse shift 
left_transverse_shift = 6*3.6576;  

% Transverse distance of the right boundary points from vehicle center 
left_transverse_distance_of_boundary_points = [left_transverse_shift*unit_ortho_vehicle_vectors_XY, zeros(length(unit_ortho_vehicle_vectors_XY),1)];


longitudinal_shift = 5; 
% Shift
longitudinal_shift_distance = [unit_vehicle_change_in_pose_XY*longitudinal_shift, zeros(length(unit_vehicle_change_in_pose_XY),1)]; 

% Left boundary points of the driven path
left_boundary_points = VehiclePose(:,1:3) + left_transverse_distance_of_boundary_points - longitudinal_shift_distance; 

% right boundary points of the driven path
right_boundary_points = VehiclePose(:,1:3) - right_transverse_distance_of_boundary_points - longitudinal_shift_distance; 

% Find the boundary points
boundary_points_of_domain = [right_boundary_points(station_1:station_2,1:3);
    flipud(left_boundary_points(station_1:station_2,1:3));
    right_boundary_points(station_1,1:3)];

% "inpolygon" is used to find the concatenated points within the boundary
[in_domain,on] = inpolygon(concatenate_LiDAR_XYZ_points(:,1),concatenate_LiDAR_XYZ_points(:,2),boundary_points_of_domain(:,1),boundary_points_of_domain(:,2));

concatenate_LiDAR_XYZ_points_new = concatenate_LiDAR_XYZ_points(in_domain,:);

%% Plotting

LLA_fig_num = 13342;
figure(LLA_fig_num);clf;

LLA_VehiclePose = gps_object.ENU2WGSLLA(VehiclePose(:,1:3));

% Do plot, and set it up so that zoom is correct, tick marks are correct,
% map is centered where we want, etc. To see options, do the geoplot, then
% zoom into the location we want, and then do:
%
% format long % This sets the format so we can see all the digits
% temp = gca  % This Gets Current Axis (GCA) and saves it as temp
%
% Next, click on "show all properties" and copy out the ZoomLevel and
% MapCenter values into the correct locations below. This way, every time
% we run this script, it automatically zooms and centers onto the correct
% location.

h_geoplot = geoplot(LLA_VehiclePose(:,1),LLA_VehiclePose(:,2),'-','Color',[0 0 1],'MarkerSize',10);
hold on;
h_parent =  get(h_geoplot,'Parent');
set(h_parent,'ZoomLevel',20.5,'MapCenter',[40.865718697633348 -77.830965127435817]);
geobasemap satellite
geotickformat -dd  % Sets the tick marks to decimal format, not degrees/minutes/seconds which is default


% NOTE: transition the geoplotting to use the PlotTestTrack library 
% Plot start and end points
geoplot(LLA_VehiclePose(1,1),LLA_VehiclePose(1,2),'.','Color',[0 1 0],'MarkerSize',10);
geoplot(LLA_VehiclePose(end,1),LLA_VehiclePose(end,2),'o','Color',[1 0 0],'MarkerSize',10);

% Plot the LIDAR in LLA
% Use the class to convert LLA to ENU
concatenate_LiDAR_LLA_points = gps_object.ENU2WGSLLA(concatenate_LiDAR_XYZ_points(:,1:3));
concatenate_LiDAR_LLA_points_new = gps_object.ENU2WGSLLA(concatenate_LiDAR_XYZ_points_new(:,1:3));

% Plot the LLA of LIDAR points
figure(LLA_fig_num);

if 1==0
    % Plot the LIDAR data simply as magenta and black points
    geoplot(concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),'mo','MarkerSize',10);
    geoplot(concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),'k.','MarkerSize',10);
else
    scaling = 3;
    intensity_fraction = scaling*concatenate_LIDAR_intensity(in_domain,:)/(intensity_max - intensity_min);
      
    % Use user-defined colormap_string to map intensity to colors. For a
    % full example, see fcn_geometry_fillColorFromNumberOrName
    old_colormap = colormap;
    % color_ordering = colormap('hot');
    color_ordering = flipud(colormap('sky'));
    colormap(old_colormap);
    N_colors = length(color_ordering(:,1));

    % Make sure the plot number is a fraction between 0 and 1
    plot_number = min(max(0,intensity_fraction),1);

    % Convert the plot number to a row
    color_row = floor((N_colors-1)*plot_number) + 1;


    % Plot the LIDAR data with intensity
    for ith_color = min(color_row):max(color_row)
        % Find the color
        color_vector = color_ordering(ith_color,:);

        % Find all the points that are in this color
        index_in_this_color = find(color_row==ith_color);

        % geoplot(concatenate_LiDAR_LLA_points(index_in_this_color,1),concatenate_LiDAR_LLA_points(index_in_this_color,2), '.','Color',color_vector,'MarkerSize',5);

        geoplot(concatenate_LiDAR_LLA_points_new(index_in_this_color,1),concatenate_LiDAR_LLA_points_new(index_in_this_color,2), '.','Color',color_vector,'MarkerSize',5);
    end
end

% Plot the vehicle pose on top of this
geoplot(LLA_VehiclePose(station_1:station_2,1),LLA_VehiclePose(station_1:station_2,2),'.','Color',[1 1 0],'MarkerSize',10);

%%

boundary_points_of_domian_LLA = gps_object.ENU2WGSLLA(boundary_points_of_domain);

geoplot(boundary_points_of_domian_LLA(:,1),boundary_points_of_domian_LLA(:,2),'r.','MarkerSize',30);