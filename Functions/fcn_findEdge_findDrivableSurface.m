function [LIDAR_ENU_under_vehicle] = fcn_findEdge_findDrivableSurface (LIDAR_ENU, VehiclePose_ENU, VehiclePose_UnitOrthoVectors)
%% fcn_findEdge_findDrivableSurface
% Find drivable points of mapping vehicle in XYZ coordinates
% 
% FORMAT:
%
%      [LIDAR_ENU_under_vehicle] = fcn_findEdge_findDrivableSurface (LIDAR_ENU, VehiclePose_ENU, VehiclePose_UnitOrthoVectors)
%
% INPUTS:     
%       
%      LIDAR_ENU: Dataset that will be plotted as an [Nx3] vector of XYZ
%      points.
%
%      VehiclePose_ENU: the position of the mapping vehicle during mapping
%      in ENU coordinates as an [Nx3] vector. The first column is
%      EAST: the direction pointing towards the geographic east. The second
%      column is NORTH: the direction pointing towards the geographic
%      north. The third column is UP (U): the direction pointing vertically
%      upwards, perpendicular to the Earth's surface.
%
%      VehiclePose_UnitOrthoVectors: Orthogonal unit vectors of the
%      position of the mapping vehicle during mapping as a [Nx3] vector,
%      with one orthogonal vector per position. The first two columns are
%      the X and Y components of the unit orthogonal vectors, e.g. the
%      direction defined as positive - at that location - in the transverse
%      or "left" direction. The third column contains all zero unit vectors
%      because the elevation of the mapping vehicle is not considered when
%      calculating whether objects are to the right/left of the vehicle.
%      Thus, this output is used typically to determine whether coordinates
%      are to the right or left of the vehicle.
%
% OUTPUTS:
%
%       LIDAR_ENU_under_vehicle: drivable points of mapping vehicle in 
%       XYZ coordinates. Data is an [Nx3] vector. 
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_findEdge_findDrivableSurface for a full
%       test suite.
%
% This function was written on 2024_08_11 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu
%
%%% Revision history

%% Debugging and Input checks
% section is empty since there is no plot


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
% section is empty since there is no plot

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

% Calculate the vectors
vector_from_vehicle_pose_to_LIDAR_points = LIDAR_ENU - VehiclePose_ENU;

% Calculate the transverse distance
transverse_only_LIDAR_points = sum(vector_from_vehicle_pose_to_LIDAR_points.*VehiclePose_UnitOrthoVectors,2);

% Define lane width limits. Take 40 percent of the lane width. Numbers
% obtained by converting 12 ft (standard lane) to meters
lane_half_width = (3.6576/2) * 0.40;  

indicies_under_vehicle = abs(transverse_only_LIDAR_points)<lane_half_width;
LIDAR_ENU_under_vehicle = LIDAR_ENU(indicies_under_vehicle,:);

end

% This must be done in ENU coordinates because LLA is not an orthogonal
% coordinate system. To find the surface, we find the distance of the lidar
% points in the orthogonal direction by taking the dot product of the LIDAR
% points, relative to the vehicle center with the unit projection vector
% pointed to the left of the vehicle.
% Jiabao
% gps_object = GPS();
% % Calculate the vectors
% vector_from_vehicle_pose_to_LIDAR_points = LIDAR_ENU - VehiclePose_ENU;
% 
% % Calculate the transverse distance
% transverse_only_LIDAR_points = sum(vector_from_vehicle_pose_to_LIDAR_points.*VehiclePose_UnitOrthoVectors,2);
% 
% % Define lane width limits. Take 40 percent of the lane width. Numbers
% % obtained by converting 12 ft (standard lane) to meters
% lane_half_width = (3.6576/2) * 0.40;  
% 
% indicies_under_vehicle = find(abs(transverse_only_LIDAR_points)<lane_half_width);
% LIDAR_ENU_under_vehicle = LIDAR_ENU(indicies_under_vehicle,:);