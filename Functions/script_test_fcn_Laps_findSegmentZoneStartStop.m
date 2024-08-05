% script_test_fcn_Laps_findSegmentZoneStartStop.m
% tests fcn_Laps_findSegmentZoneStartStop.m

% Revision history
%     2022_07_12
%     -- first write of the code

close all
clc

% Create some data to plot
full_steps = (-1:0.1:1)';
zero_full_steps = 0*full_steps; %#ok<*NASGU>
ones_full_steps = ones(length(full_steps(:,1)),1);
half_steps = (-1:0.1:0)';
zero_half_steps = 0*half_steps;
ones_half_steps = ones(length(half_steps(:,1)),1); %#ok<*PREALL>


%% Check assertions
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

%% This one returns nothing since there is no portion of the path in the
% criteria, even though the path goes right over the criteria
fig_num = 1;

query_path = ...
    [full_steps 0.4*ones_full_steps];

segment_definition = [0 0; 1 0]; % Starts at [0,0], ends at [1 0]
[zone_start_indices, zone_end_indices] = ...
    fcn_Laps_findSegmentZoneStartStop(...
    query_path,...
    segment_definition,...
    fig_num);

assert(isempty(zone_start_indices));
assert(isempty(zone_end_indices));

%% This one returns one hit
fig_num = 1;

query_path = ...
    [full_steps 0.4*ones_full_steps];

segment_definition = [0 0; 0 1]; % Starts at [0,0], ends at [0 1]
[zone_start_indices, zone_end_indices] = ...
    fcn_Laps_findSegmentZoneStartStop(...
    query_path,...
    segment_definition,...
    fig_num);

assert(isequal(zone_start_indices,10));
assert(isequal(zone_end_indices,11));

%% This one returns one hit, right on start of path
fig_num = 2;

query_path = ...
    [full_steps 0.4*ones_full_steps];

segment_definition = [-1 0; -1 1]; % Starts at [-1,0], ends at [-1 1]
[zone_start_indices, zone_end_indices] = ...
    fcn_Laps_findSegmentZoneStartStop(...
    query_path,...
    segment_definition,...
    fig_num);

assert(isequal(zone_start_indices,1));
assert(isequal(zone_end_indices,2));

%% This one returns no hit, crossed wrong way
fig_num = 3;

query_path = ...
    [flipud(full_steps) 0.4*ones_full_steps];

segment_definition = [0 0; 0 1]; % Starts at [0,0], ends at [0 1]
[zone_start_indices, zone_end_indices] = ...
    fcn_Laps_findSegmentZoneStartStop(...
    query_path,...
    segment_definition,...
    fig_num);

assert(isempty(zone_start_indices));
assert(isempty(zone_end_indices));

%% This one returns no hit, also crossed wrong way
fig_num = 31;

query_path = ...
    [full_steps 0.4*ones_full_steps];

segment_definition = [0 1; 0 0]; % Starts at [0,1], ends at [0 0]
[zone_start_indices, zone_end_indices] = ...
    fcn_Laps_findSegmentZoneStartStop(...
    query_path,...
    segment_definition,...
    fig_num);

assert(isempty(zone_start_indices));
assert(isempty(zone_end_indices));


%% This one returns two hits, even though crossed three times
% One crossing is in the wrong direction!
fig_num = 4;

query_path = ...
    [full_steps 0.4*ones_full_steps; 
    flipud(full_steps) 0.6*ones_full_steps;
    full_steps 0.8*ones_full_steps];

segment_definition = [0 0; 0 1]; % Starts at [0,0], ends at [0 1]
[zone_start_indices, zone_end_indices] = ...
    fcn_Laps_findSegmentZoneStartStop(...
    query_path,...
    segment_definition,...
    fig_num);

assert(isequal(zone_start_indices,[10; 52]));
assert(isequal(zone_end_indices,[11; 53]));

%% This is a hard one
% Multiple crossings
fig_num = 5;

query_path = ...
    [-1 2; 1 2; -1 1; 0 1; -1 0; 0 0; 1 0; 0 -1; 1 -1; -1 -2; 0 -2];

segment_definition = [0 -2; 0 2]; % Starts at [0,0], ends at [0 1]
[zone_start_indices, zone_end_indices] = ...
    fcn_Laps_findSegmentZoneStartStop(...
    query_path,...
    segment_definition,...
    fig_num);

assert(isequal(zone_start_indices,[1; 3; 5; 8; 10]));
assert(isequal(zone_end_indices,[2; 4; 6; 9; 11]));



%% Fail conditions
if 1==0
    
    %% Fails because segment_definition is not correct type
    clc
    
    segment_definition = [0; 1]; % Starts at ?, ends at ?
    [~, ~] = ...
        fcn_Laps_findSegmentZoneStartStop(...
        query_path,...
        segment_definition,...
        []);
    
    %% Warns because segment_definition is 3D
    clc
    
    segment_definition = [0 0 0; 1 0 0]; % Starts at [0 0 0], ends at [1 0 0]
    [~, ~] = ...
        fcn_Laps_findSegmentZoneStartStop(...
        query_path,...
        segment_definition,...
        []);
   
end