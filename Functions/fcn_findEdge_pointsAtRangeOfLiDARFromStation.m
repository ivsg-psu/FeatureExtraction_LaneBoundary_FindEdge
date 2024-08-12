function [station1_minus_range_index, station2_plus_range_index] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,starting_index,ending_index,varargin)
%% fcn_findEdge_pointsAtRangeOfLiDARFromStation
% This code finds the index of a point that is a certain range before
% station 1 and the index of a point that is a certain range after
% station 2. 
%
% FORMAT:
% [station1_minus_range_index, station2_plus_range_index]... 
% = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,...
% starting_index,ending_index,(range))
% 
% INPUTS:
%       Vehicle Pose - Position of the Vehicle
%
%       starting_index - Index of the starting point of station 1
%
%       ending_index - Index of the ending of station 2
%
% (OPTIONAL INPUTS):
%
%       range_of_LiDAR - range before/after the indexes in meters that is
%       looked at. Default is 100.
%
% OUTPUTS:
% 
%       station1_minus_range_index - starting index after range transformation
%
%       station2_plus_range_index - ending index ater range transformation
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
% REVISION HISTORY:
% 2024_08_02 - Aneesh Batchu
% -- wrote the code originally
% 2024_08_12 - Aleksandr Goncharov
% -- functionalized the code

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

% Does user want to specify the range_of_LiDAR?
range_of_LiDAR = 100;
if 4 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        range_of_LiDAR = temp;
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
station_1 = starting_index; 
station_2 = ending_index;

% Calculate the differences between consecutive points
differences = diff(vehicle_positionsXY);

% Compute the Euclidean distances for each pair of consecutive points
distances = sqrt(sum(differences.^2, 2));

% Compute the cumulative sum of the distances
cumulative_distances = [0; cumsum(distances)];

% Find the cumulative distance of station 1 and station 2
station1_distance = cumulative_distances(station_1);
station2_distance = cumulative_distances(station_2);

% Find the index of the point before station 1 whose distance is
% approximately equal to the range of the LiDAR.
station1_minus_range_index = find(cumulative_distances <= station1_distance - range_of_LiDAR, 1, 'last');
%station1_minus_range_point = vehicle_positionsXY(station1_minus_range_index, :);

% Find the index of the point after station 2 whose distance is
% approximately equal to the range of the LiDAR.
station2_plus_range_index = find(cumulative_distances >= station2_distance + range_of_LiDAR, 1, 'first');
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

% (NOTHING TO PLOT)

end