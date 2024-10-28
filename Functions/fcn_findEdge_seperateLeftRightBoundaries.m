function [boundary_points_left, boundary_points_right] = fcn_findEdge_seperateLeftRightBoundaries...
    (VehiclePose,scanLineNumber_start,scanLineNumber_end,nearest_boundary_points, varargin)
%% fcn_findEdge_seperateLeftRightBoundaries    Seperate the right and left boundaries from the nearest boundaries
% 
% FORMAT:
%
%    [boundary_points_left, boundary_points_right] = fcn_findEdge_seperateLeftRightBoundaries...
%    (VehiclePose,station_1,station_2,nearest_boundary_points, grid_size,
%    transverse_shift, (fig_num)).
%
% INPUTS:
%      VehiclePose: the position of the mapping vehicle during mapping.    
%
%      scanLineNumber_start: the starting scan line of the LIDAR data
%
%      scanLineNumber_end: the end of scan line of the LIDAR data
%
%      nearest_boundary_points: the nearest boundary points close to the
%      driven path of of the mapping vehicle. 
%
%
% (OPTIONAL INPUTS):
%
%
%      (fig_num): figure number, can be set to -1 for fast mode
%
% OUTPUTS: 
%
%      boundary_points_left: the boundary points on the left side of the of
%      the mapping vehcile 
%
%      boundary_points_right: the boundary points on the right side of the of
%      the mapping vehcile 
%
%
% DEPENDENCIES:
%
% EXAMPLES:
%   See the script:
%   
%   script_test_fcn_findEdge_seperateLeftRightBoundaries
%
%   for a full test suite
%
% This function was written by Jiabao Zhao on 2024_08_14
%
% Revision history
% 2024_10_24 - Aneesh Batchu
% -- wrote the code originally 


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

if flag_max_speed == 0
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(4,5);
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

% Path library does not have fast mode options.
if flag_max_speed == 1
    fig_num = [];
end

St_points = fcn_Path_convertXY2St(VehiclePose(scanLineNumber_start:scanLineNumber_end,1:2),nearest_boundary_points(:,1:2),[],fig_num);
transverse_error = real(St_points(:,2));

% Right Boundary Points
boundary_points_right = nearest_boundary_points(transverse_error < 0,1:2);

% Left Boundary Points
boundary_points_left = nearest_boundary_points(transverse_error >= 0,1:2); 


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
% No plots here
end
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
