function [concatenate_LiDAR_XYZ_points_new, boundary_points_of_domain, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2)
%% fcn_findEdge_findPointsInDomain
% Find the LIDAR_ENU and LIDAR_scanLineAndRingID in domain
% 
% FORMAT:
%
%      [concatenate_LiDAR_XYZ_points_new, boundary_points_of_domain, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2)
%
% INPUTS:     
%       
%      VehiclePose: the position of the mapping vehicle during mapping
%
%      LIDAR_ENU: LiDAR data points corresponding to the intensity
%
%      station_1: the starting scan line of the LIDAR data
%
%      station_2: the ending scan line of the LIDAR data
%      
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      concatenate_LiDAR_XYZ_points_new: LIDAR data that are within the
%      range of the boundary. This is [Nx3] matrix 
%
%      boundary_points_of_domain: boundary points that are within the
%      domain.
%
%      in_domain: For points within the boundary, a value of 1 indicates 
%      that the point is inside the boundary, while a value of 0 indicates 
%      it is outside. 
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_findEdge_findPointsInDomain for a full
%       test suite.
%
% This function was written on 2024_08_12 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu
%
%%% Revision history
%

%% flag_check_inputs

flag_check_inputs=1;
%% check input arguments?
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
    narginchk(3,4); 
end

%% Solve for the Maxs and Mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vehicle_change_in_pose_XY = diff(VehiclePose(:,1:2));

% Repeat the last value again, since diff removes one row. We want the same
% number of vectors as the number of points, and diff removed one point.
vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];

% Convert these to unit vectors
unit_vehicle_change_in_pose_XY = fcn_findEdge_calcUnitVector(vehicle_change_in_pose_XY);

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
[in_domain,~] = inpolygon(LIDAR_ENU(:,1),LIDAR_ENU(:,2),boundary_points_of_domain(:,1),boundary_points_of_domain(:,2));

concatenate_LiDAR_XYZ_points_new = LIDAR_ENU(in_domain,:);


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

% no plot
%%