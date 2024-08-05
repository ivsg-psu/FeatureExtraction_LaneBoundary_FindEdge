% script_test_fcn_Laps_checkZoneType.m
% tests fcn_Laps_checkZoneType.m

% Revision history
%     2022_07_23 - sbrennan@psu.edu
%     -- wrote the code originally, using breakDataIntoLaps as starter

%% Set up the workspace
close all
clc


%% Check assertions for basic path operations and function testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              _   _                 
%      /\                     | | (_)                
%     /  \   ___ ___  ___ _ __| |_ _  ___  _ __  ___ 
%    / /\ \ / __/ __|/ _ \ '__| __| |/ _ \| '_ \/ __|
%   / ____ \\__ \__ \  __/ |  | |_| | (_) | | | \__ \
%  /_/    \_\___/___/\___|_|   \__|_|\___/|_| |_|___/
%                                                    
%                                                    
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Assertions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This is a standard call for point-type
input_start_zone_definition = [2 3 0 0]; % Radius of 2, 3 points, centered at 0 0
[flag_start_is_a_point_type, output_start_zone_definition] = ...
    fcn_Laps_checkZoneType(input_start_zone_definition, 'start_definition');

% Make sure its type is correct
assert(isequal(1,flag_start_is_a_point_type))

% Make sure the output is correct
assert(isequal(output_start_zone_definition,[2 3 0 0]))

%% This is a standard call for segment-type
input_start_zone_definition = [2 3; 0 0]; % starts at 2 3, ends at 0 0
[flag_start_is_a_point_type, output_start_zone_definition] = ...
    fcn_Laps_checkZoneType(input_start_zone_definition, 'start_definition');

% Make sure its type is correct
assert(isequal(0,flag_start_is_a_point_type))

% Make sure the output is correct
assert(isequal(output_start_zone_definition,[2 3; 0 0]))

%% Fail conditions
if 1==0
    %% WARNING for point-type, due to 3D
    input_start_zone_definition = [2 3 0 0 0]; % Radius of 2, 3 points, centered at 0 0 0
    [flag_start_is_a_point_type, output_start_zone_definition] = ...
        fcn_Laps_checkZoneType(input_start_zone_definition, 'start_definition');
    
    % Make sure its type is correct
    assert(isequal(1,flag_start_is_a_point_type))
    
    % Make sure the output is correct
    assert(isequal(output_start_zone_definition,[2 3 0 0]))
    
    %% ERROR for point-type, due to bad array size
    input_start_zone_definition = [2 3]; % Radius of 2, 3 points, centered at ???
    [~, ~] = ...
        fcn_Laps_checkZoneType(input_start_zone_definition, 'start_definition');
    
   
    %% ERROR for point-type, due to bad array size
    input_start_zone_definition = [2 3 4 5 6 7 8]; % Radius of 2, 3 points, centered at ???
    [~, ~] = ...
        fcn_Laps_checkZoneType(input_start_zone_definition, 'start_definition');
    
    
    %% WARNING for segment-type, due to 3D
    input_start_zone_definition = [2 3 0; 0 0 0]; % starts at 2 3 0, ends at 0 0 0
    [flag_start_is_a_point_type, output_start_zone_definition] = ...
        fcn_Laps_checkZoneType(input_start_zone_definition, 'start_definition');
    
    % Make sure its type is correct
    assert(isequal(0,flag_start_is_a_point_type))
    
    % Make sure the output is correct
    assert(isequal(output_start_zone_definition,[2 3; 0 0]))
    
    %% ERROR for segment-type, due to bad array size
    input_start_zone_definition = [2 3 0 4; 0 0 0 0]; % starts at ???, ends at ???
    [~, ~] = ...
        fcn_Laps_checkZoneType(input_start_zone_definition, 'start_definition');
    
    
    %% ERROR for segment-type, due to bad array size
    input_start_zone_definition = [2; 3]; % starts at ????, ends at ????
    [~, ~] = ...
        fcn_Laps_checkZoneType(input_start_zone_definition, 'start_definition');
    

end
