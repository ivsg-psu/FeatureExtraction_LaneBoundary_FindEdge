function [scanLineStart_minus_range_index, scanLineEnd_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,starting_index,ending_index,varargin)
%% fcn_findEdge_pointsAtRangeOfLiDARFromStation
%
% This code finds the index of a point that is a certain range before the
% start scan line and the index of a point that is a certain range after
% the end scan line.
%
% FORMAT:
% [station1_minus_range_index, station2_plus_range_index]... 
% = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,...
% starting_index,ending_index,(range))
% 
% INPUTS:
%       VehiclePose: Position and oriendtation of the vehicle
%
%       starting_index: Start scan line
%
%       ending_index: End scan line
%
% (OPTIONAL INPUTS):
%
%       range_of_LiDAR: range before/after the indexes in meters that is
%       looked at. Default is 100.
%
% OUTPUTS:
% 
%       station1_minus_range_index: Starting index after range transformation
%
%       station2_plus_range_index: Ending index ater range transformation
%
%
% DEPENDENCIES:
% (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_findEdge_pointsAtRangeOfLiDARFromStation
% for a full test suite.
%
% Revision History:
% 2024_08_02 - Aneesh Batchu
% -- wrote the code originally
% 2024_08_12 - Aleksandr Goncharov
% -- functionalized the code
% 2024_08_20 - Aneesh Batchu
% -- Modified instructions

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
    narginchk(3,5); 
end

% Does user want to specify the range_of_LiDAR?
range_of_LiDAR = 100;
if 4 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        range_of_LiDAR = temp;
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

% Define the vehicle position
vehicle_positionsXY = VehiclePose(:,1:2); 

%define the start/end points
scanLineStart = starting_index; 
scanLineEnd = ending_index;

% Calculate the differences between consecutive points
differences = diff(vehicle_positionsXY);

% Compute the Euclidean distances for each pair of consecutive points
distances = sqrt(sum(differences.^2, 2));

% Compute the cumulative sum of the distances
cumulative_distances = [0; cumsum(distances)];

% Find the cumulative distance of the start scan line and the end scan line
scanLineStart_cumulative_distance = cumulative_distances(scanLineStart);
scanLineEnd_cumulative_distance = cumulative_distances(scanLineEnd);

% Find the index of the point before station 1 whose distance is
% approximately equal to the range of the LiDAR.
scanLineStart_minus_range_index = find(cumulative_distances <= scanLineStart_cumulative_distance - range_of_LiDAR, 1, 'last');
%station1_minus_range_point = vehicle_positionsXY(station1_minus_range_index, :);

% Find the index of the point after station 2 whose distance is
% approximately equal to the range of the LiDAR.
scanLineEnd_plus_range_index = find(cumulative_distances >= scanLineEnd_cumulative_distance + range_of_LiDAR, 1, 'first');
%station2_plus_range_point = vehicle_positionsXY(station2_plus_range_index, :);

%%
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
    % (NOTHING TO PLOT)
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function

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
