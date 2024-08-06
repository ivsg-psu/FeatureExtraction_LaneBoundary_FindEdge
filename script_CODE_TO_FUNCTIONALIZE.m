



%% Plot the LIDAR in ENU
figure(ENU_fig_num);

plot(LIDAR_ENU(:,1),LIDAR_ENU(:,2),'.','Color',[0 0 1],'MarkerSize',1);

% % Plot start and end points
% plot(VehiclePose(1,1),VehiclePose(1,2),'.','Color',[0 1 0],'MarkerSize',10);
% plot(VehiclePose(end,1),VehiclePose(end,2),'o','Color',[1 0 0],'MarkerSize',10);

%% Plot the LIDAR in 3D ENU
ENU_3D_fig_num = 3;
figure(ENU_3D_fig_num);
clf;

% Calculate the intensity ranges
intensity_min = min(LIDAR_intensity);
intensity_max = max(LIDAR_intensity);

hold on;
grid on;
axis equal

if 1==0
    % Plot the vehicle's trajectory
    plot3(VehiclePose(:,1),VehiclePose(:,2),VehiclePose(:,3),'-','Color',[1 0 1],'MarkerSize',10,'LineWidth',3);

    % Plot start and end points of trajectory
    plot3(VehiclePose(1,1),VehiclePose(1,2),VehiclePose(1,3),'.','Color',[1 0 0],'MarkerSize',10,'LineWidth',3);
    plot3(VehiclePose(end,1),VehiclePose(end,2),VehiclePose(end,3),'o','Color',[0 1 0],'MarkerSize',10,'LineWidth',3);

else
    plot3(VehiclePose(scanLineNumber_start:scanLineNumber_end,1),VehiclePose(scanLineNumber_start:scanLineNumber_end,2),VehiclePose(scanLineNumber_start:scanLineNumber_end,3),'.','Color',[0 0 0],'MarkerSize',30,'LineWidth',3);
    plot3(VehiclePose(scanLineNumber_start:scanLineNumber_end,1),VehiclePose(scanLineNumber_start:scanLineNumber_end,2),VehiclePose(scanLineNumber_start:scanLineNumber_end,3),'.','Color',[1 1 0],'MarkerSize',10,'LineWidth',3);

    % Show the orthogonal arrows showing vehicle motion directions. Green
    % is forward, bLue is Left
    quiver3(...
        VehiclePose(scanLineNumber_start,1),VehiclePose(scanLineNumber_start,2),VehiclePose(scanLineNumber_start,3), ...
        unit_vehicle_change_in_pose_XY(scanLineNumber_start,1),unit_vehicle_change_in_pose_XY(scanLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 1 0]);
    quiver3(...
        VehiclePose(scanLineNumber_start,1),VehiclePose(scanLineNumber_start,2),VehiclePose(scanLineNumber_start,3), ...
        unit_ortho_vehicle_vectors_XY(scanLineNumber_start,1),unit_ortho_vehicle_vectors_XY(scanLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 0 1]);
end

xlabel('East position [m]');
ylabel('North position [m]');
zlabel('Up position [m]');
view(3)

if 1==0
    % Plot the LIDAR data simply as blue points
    plot3(LIDAR_ENU(:,1),LIDAR_ENU(:,2),LIDAR_ENU(:,3), '.','Color',[0 0 1],'MarkerSize',1);
else
    scaling = 3;
    intensity_fraction = scaling*LIDAR_intensity/(intensity_max - intensity_min);
    
    % Use user-defined colormap_string to map intensity to colors. For a
    % full example, see fcn_geometry_fillColorFromNumberOrName
    old_colormap = colormap;
    color_ordering = colormap('hot');
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
        plot3(...
            LIDAR_ENU(index_in_this_color,1),...
            LIDAR_ENU(index_in_this_color,2),...
            LIDAR_ENU(index_in_this_color,3), '.','Color',color_vector,'MarkerSize',5);
    end
end

% set(gca,'CameraViewAngle',6)

%% Plot the LIDAR in XZ ENU
ENU_XZ_fig_num = 4;
figure(ENU_XZ_fig_num);
clf;

hold on;
grid on;
axis equal

xlabel('East position [m]');
ylabel('Up position [m]');

% Plot the LIDAR data
plot(LIDAR_ENU(:,1),LIDAR_ENU(:,3), '.','Color',[0 0 1],'MarkerSize',5);

% plot the vehicle
plot(VehiclePose(scanLineNumber_start,1),VehiclePose(scanLineNumber_start,3),'.','Color',[0 0 0],'MarkerSize',30,'LineWidth',3);
plot(VehiclePose(scanLineNumber_start,1),VehiclePose(scanLineNumber_start,3),'.','Color',[1 1 0],'MarkerSize',10,'LineWidth',1);

%% Plot the LIDAR in YZ ENU
ENU_YZ_fig_num = 5;
figure(ENU_YZ_fig_num);
clf;

hold on;
grid on;
axis equal

xlabel('North position [m]');
ylabel('Up position [m]');

% Plot the LIDAR data
plot(LIDAR_ENU(:,2),LIDAR_ENU(:,3), '.','Color',[0 0 1],'MarkerSize',5);

% plot the vehicle
plot(VehiclePose(scanLineNumber_start,2),VehiclePose(scanLineNumber_start,3),'.','Color',[0 0 0],'MarkerSize',30,'LineWidth',3);
plot(VehiclePose(scanLineNumber_start,2),VehiclePose(scanLineNumber_start,3),'.','Color',[1 1 0],'MarkerSize',10,'LineWidth',1);

%% Plot the vehicle pose in LLA

%% Define GPS object 
% Define base coordinates - THIS ONLY WORKS FOR THE TEST TRACK!!! Change for
% other sites as needed.
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;

gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class

LLA_fig_num = 2;
figure(LLA_fig_num);
clf;

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
concatenate_LiDAR_LLA_points = gps_object.ENU2WGSLLA(LIDAR_ENU(:,1:3));


% Plot the LLA of LIDAR points
figure(LLA_fig_num);

if 1==0
    % Plot the LIDAR data simply as magenta and black points
    geoplot(concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),'mo','MarkerSize',10);
    geoplot(concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),'k.','MarkerSize',10);
else
    scaling = 3;
    intensity_fraction = scaling*LIDAR_intensity/(intensity_max - intensity_min);
      
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

        geoplot(concatenate_LiDAR_LLA_points(index_in_this_color,1),concatenate_LiDAR_LLA_points(index_in_this_color,2), '.','Color',color_vector,'MarkerSize',5);
    end
end

% Plot the vehicle pose on top of this
geoplot(LLA_VehiclePose(scanLineNumber_start:scanLineNumber_end,1),LLA_VehiclePose(scanLineNumber_start:scanLineNumber_end,2),'.','Color',[1 1 0],'MarkerSize',10);

%% Find the drivable surface
% This must be done in ENU coordinates because LLA is not an orthogonal
% coordinate system. To find the surface, we find the distance of the lidar
% points in the orthogonal direction by taking the dot product of the LIDAR
% points, relative to the vehicle center with the unit projection vector
% pointed to the left of the vehicle.

% Calculate the vectors
vector_from_vehicle_pose_to_LIDAR_points = LIDAR_ENU - VehiclePose_ENU_points;

% Calculate the transverse distance
transverse_only_LIDAR_points = sum(vector_from_vehicle_pose_to_LIDAR_points.*VehiclePose_UnitOrthoVectors,2);

% Define lane width limits. Take 40 percent of the lane width. Numbers
% obtained by converting 12 ft (standard lane) to meters
lane_half_width = (3.6576/2) * 0.40;  

indicies_under_vehicle = find(abs(transverse_only_LIDAR_points)<lane_half_width);
LIDAR_ENU_under_vehicle = LIDAR_ENU(indicies_under_vehicle,:);
concatenate_LiDAR_LLA_points_under_vehicle = gps_object.ENU2WGSLLA(LIDAR_ENU_under_vehicle(:,1:3));


% Plot the LIDAR data underneath the vehicle in XYZ
figure(ENU_3D_fig_num);
plot3(LIDAR_ENU_under_vehicle(:,1),LIDAR_ENU_under_vehicle(:,2),LIDAR_ENU_under_vehicle(:,3), '.','Color',[0 1 0],'MarkerSize',1);

% Plot the LIDAR data underneath the vehicle in LLA
figure(LLA_fig_num);
geoplot(concatenate_LiDAR_LLA_points_under_vehicle(:,1),concatenate_LiDAR_LLA_points_under_vehicle(:,2),'g.','MarkerSize',2);












