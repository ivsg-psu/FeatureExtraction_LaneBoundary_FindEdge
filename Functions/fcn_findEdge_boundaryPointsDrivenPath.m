function [left_boundary_points,right_boundary_points,boundaryLineNumber_start,boundaryLineNumber_end] = fcn_findEdge_boundaryPointsDrivenPath(VehiclePose,scanLineRange,shift,varargin)
%% fcn_findEdge_boundaryPointsDrivenPath
% This function will find the driven path (left and right side points).
% Some parts of plotting for this function is done separate from this
% function. This function only plots the quiver functions.
% 
% FORMAT: 
% [left_boundary_points,right_boundary_points,boundaryLineNumber_start,boundaryLineNumber_end]...
% = fcn_findEdge_boundaryPointsDrivenPath(VehiclePose,scanLineRange,shift,(fig_num))
%
% INPUTS:
%       VehiclePose: Vehicle Position Data
%
%       scanLineRange: Range of Lines in format: [lowEnd highEnd]
%
%       shift
%
% (OPTIONAL INPUTS):
%
%       fig_num: number of the figure
%
%
% DEPENDENCIES:
%  (None)
% 
% EXAMPLES:
%   
%   See the script:
%   
%   script_test_fcn_findEdge_boundaryPointsDrivenPath
%
%   for a full test suite
%
% This function was written on 2024_07_18 by Aneesh Batchu
%
% REVISION HISTORY:
% 2024_08_12 -- Aleksandr Goncharov - opg5041@psu.edu - 267-304-8354
% -- Functionalized the code w/Test Script

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
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

% Check correct number of inputs.
if flag_max_speed == 0
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(3,4);
    end
end 

%default shift value if none given
if isempty(shift)
    shift=5;
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (4<=nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

%%
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
unit_vehicle_change_in_pose_XY = fcn_INTERNAL_calcUnitVector(vehicle_change_in_pose_XY);

% Find orthogonal vetors by rotating by 90 degrees in the CCW direction
unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY*[0 1; -1 0];

% Define lane width limits. Take 40 percent of the lane width. Numbers
% obtained by converting 12 ft (standard lane) to meters
lane_half_width = (3.6576/2) * 0.40;  

% Transverse distance of the left and right boundary points of the driven
% path from vehicle center 
transverse_distance_of_boundary_points = [lane_half_width*unit_ortho_vehicle_vectors_XY, zeros(length(unit_ortho_vehicle_vectors_XY),1)];

%Shift 
shift_distance = [unit_vehicle_change_in_pose_XY*shift, zeros(length(unit_vehicle_change_in_pose_XY),1)]; 

% Left boundary points of the driven path
left_boundary_points = VehiclePose(:,1:3) + transverse_distance_of_boundary_points - shift_distance; 

% right boundary points of the driven path
right_boundary_points = VehiclePose(:,1:3) - transverse_distance_of_boundary_points - shift_distance;

scanLineNumber_start = scanLineRange(1);
scanLineNumber_end   = scanLineRange(2);

Nscans = length(VehiclePose(:,1));

boundaryLineNumber_start = max(scanLineNumber_start-1,1); 
boundaryLineNumber_end   = min(scanLineNumber_end+1, Nscans);

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
    figure(fig_num)
   
    %plot the quivers
    % Show the orthogonal arrows showing vehicle motion directions. Green
    % is forward, Blue is Left
    
    quiver3(...
        VehiclePose(boundaryLineNumber_start,1),VehiclePose(boundaryLineNumber_start,2),VehiclePose(boundaryLineNumber_start,3), ...
        unit_vehicle_change_in_pose_XY(boundaryLineNumber_start,1),unit_vehicle_change_in_pose_XY(boundaryLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 0 1]);
    quiver3(...
        VehiclePose(boundaryLineNumber_start,1),VehiclePose(boundaryLineNumber_start,2),VehiclePose(boundaryLineNumber_start,3), ...
        unit_ortho_vehicle_vectors_XY(boundaryLineNumber_start,1),unit_ortho_vehicle_vectors_XY(boundaryLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 0 1]);
end

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

%% fcn_INTERNAL_calcUnitVector
function unit_vectors = fcn_INTERNAL_calcUnitVector(input_vectors)
vector_length = sum(input_vectors.^2,2).^0.5;
unit_vectors = input_vectors./vector_length;
end

end