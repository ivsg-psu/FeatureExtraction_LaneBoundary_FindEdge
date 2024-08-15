function [concatenate_LiDAR_XYZ_points_new, boundary_points_of_domain, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2, LIDAR_intensity, varargin)
%% fcn_findEdge_findPointsInDomain
% Find the LIDAR_ENU and LIDAR_scanLineAndRingID in domain
% 
% FORMAT:
%
%      [concatenate_LiDAR_XYZ_points_new, boundary_points_of_domain, in_domain] = fcn_findEdge_findPointsInDomain(VehiclePose, LIDAR_ENU, station_1, station_2,LIDAR_intensity,(fig_num))
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
%      LIDAR_intensity
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

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS");
    MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG = getenv("MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG); 
        flag_check_inputs  = str2double(MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
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

if 0==flag_max_speed
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(5,6);

        % % Check the points input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points, '2column_of_numbers',[2 3]);
        %
        % % Check the transverse_tolerance input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);
        %
        % % Check the station_tolerance input is a positive single number
        % if ~isempty(station_tolerance)
        %     fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);
        % end
    end
end



% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (5<=nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
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
if flag_do_plots
    % plot vehicle trajectory in LLA
    reference_LLA = [];
    zoom_in_location = [40.865718697633348 -77.830965127435817];
    zoomLevel = 20.5;
    [~] = fcn_findEdge_plotVehicleLLA(VehiclePose(:,1:3), (reference_LLA), (zoom_in_location), (zoomLevel), (fig_num));

    % plot LLA of LIDAR points
    marker_size = [];
    scaling = [];
    color_map = 'sky';
    fcn_findEdge_plotLIDARLLA_Aneesh(LIDAR_ENU, LIDAR_intensity, in_domain, concatenate_LiDAR_XYZ_points_new, (scaling),(color_map),(marker_size),(reference_LLA),(fig_num))

    % Plot the vehicle pose
    format = sprintf('''.'',''Color'',[1 1 0],''MarkerSize'',10');
    LIDAR_intensity1 = [];
    fcn_findEdge_plotLIDARLLA(VehiclePose(station_1:station_2,1:3),(LIDAR_intensity1),(scaling),(color_map),(marker_size),(reference_LLA),(format),(fig_num))


    % Plot the boundary points
    format = sprintf('''r.'',''MarkerSize'',30');
    LIDAR_intensity1 = [];
    fcn_findEdge_plotLIDARLLA(boundary_points_of_domain,(LIDAR_intensity1),(scaling),(color_map),(marker_size),(reference_LLA),(format),(fig_num))
end
% no plot
%%