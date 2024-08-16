function [LIDAR_ENU_under_vehicle] = fcn_findEdge_findDrivableSurface (LIDAR_ENU, VehiclePose_ENU, VehiclePose_UnitOrthoVectors,varargin)
%% fcn_findEdge_findDrivableSurface
%  Find the driven path points in LIDAR scans
% 
% FORMAT:
%
%      [LIDAR_ENU_under_vehicle] = fcn_findEdge_findDrivableSurface (LIDAR_ENU, VehiclePose_ENU, VehiclePose_UnitOrthoVectors,(fig_num),(fig_num2))
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

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
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
        narginchk(3,5);

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
if (0==flag_max_speed) &&  (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

% Does user want to specify another fig_num?
if (0==flag_max_speed) &&  (4<=nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num2 = temp;
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

% Calculate the vectors
vector_from_vehicle_pose_to_LIDAR_points = LIDAR_ENU - VehiclePose_ENU;

% Calculate the transverse distance
transverse_only_LIDAR_points = sum(vector_from_vehicle_pose_to_LIDAR_points.*VehiclePose_UnitOrthoVectors,2);

% Define lane width limits. Take 40 percent of the lane width. Numbers
% obtained by converting 12 ft (standard lane) to meters
lane_half_width = (3.6576/2) * 0.40;  

indicies_under_vehicle = abs(transverse_only_LIDAR_points)<lane_half_width;
LIDAR_ENU_under_vehicle = LIDAR_ENU(indicies_under_vehicle,:);

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
    % Plot the LIDAR data underneath the vehicle in XYZ
    LIDAR_intensity = [];
    scaling = [];
    color_map = [];
    format = sprintf(' ''.'',''Color'',[0 1 0],''MarkerSize'', 20');
    fcn_findEdge_plotLIDARXYZ(LIDAR_ENU_under_vehicle, (LIDAR_intensity), (scaling), (color_map), (format), (fig_num));
    daspect([1 1 0.1]); % Adjust aspect ratio

    % Plot the LIDAR data underneath the vehicle in LLA
    reference_LLA = [];
    marker_size = [];
    format = sprintf(' ''.'',''Color'',[0 1 0],''MarkerSize'', 2');
    fcn_findEdge_plotLIDARLLA(LIDAR_ENU_under_vehicle,(LIDAR_intensity),(scaling),(color_map),(marker_size),(reference_LLA),(format),(fig_num2));
end
end
