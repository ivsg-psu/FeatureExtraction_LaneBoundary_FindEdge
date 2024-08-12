
%% script_test_geometry_boundaryPointsDrivenPath
% This script is written to find the boundary points of driven path
%
% 2024_07_18 - Aneesh Batchu
% -- wrote the code originally
% 2024_07_20 - Aneesh Batchu
% --  used "shift" (a varible name) to shift the boundary points

%% Calculate the vehicle orientation

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
lane_half_width = (3.6576/2) * 0.40;  
% lane_half_width = 2*3.6576;  

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

scanLineNumber_start = 1400;
scanLineNumber_end = 1450;

boundaryLineNumber_start = scanLineNumber_start-1;  %scanLineNumber_start - 8; 
boundaryLineNumber_end = scanLineNumber_end+1;      %scanLineNumber_end - 6; 

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

%%

ENU_3D_fig_num = 3;
figure(ENU_3D_fig_num);

plot3(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1),left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,2),left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 1 0],'MarkerSize',30,'LineWidth',3);
plot3(right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1),right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,2),right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 0 1],'MarkerSize',30,'LineWidth',3);

boundary_points_driven_path = [right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3);
    flipud(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3));
    right_boundary_points(boundaryLineNumber_start,1:3)];

boundary_points_driven_path_LLA = gps_object.ENU2WGSLLA(boundary_points_driven_path);

LLA_fig_num = 13342;
figure(LLA_fig_num);
% Plot the LIDAR data underneath the vehicle in LLA
figure(LLA_fig_num);
geoplot(boundary_points_driven_path_LLA(:,1),boundary_points_driven_path_LLA(:,2),'g.','MarkerSize',30);


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