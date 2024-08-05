% script_test_fcn_Laps_breakDataIntoLaps.m
% tests fcn_Laps_breakDataIntoLaps.m

% Revision history
%     2022_04_03
%     -- first write of the code
%     2022_04_03
%     -- added external call to zone calculation function
%     2022_07_11 - sbrennan@psu.edu
%     -- corrected calls to zone function to allow number of points,
%     changed format to allow 3d circles

%% Set up the workspace
close all
clc
clear laps_array data single_lap

%% Load some test data and plot it in figure 1 
% Call the function to fill in an array of "path" type
laps_array = fcn_Laps_fillSampleLaps;


% Convert paths to traversals structures. Each traversal instance is a
% "traversal" type, and the array called "data" below is a "traversals"
% type.
for i_Path = 1:length(laps_array)
    traversal = fcn_Path_convertPathToTraversalStructure(laps_array{i_Path});
    data.traversal{i_Path} = traversal;
end

% Plot the last one
fig_num = 1;
example_lap_data.traversal{1} = data.traversal{end};
fcn_Laps_plotLapsXY(example_lap_data,fig_num);

%% Call the function to show it operating, and plot in figure 2
start_definition = [10 3 0 0]; % Radius 10, 3 points must pass near [0,0]
end_definition = [30 3 0 -60]; % Radius 30, 3 points must pass near [0,-60]
excursion_definition = []; % empty
fig_num = 2;
lap_traversals = fcn_Laps_breakDataIntoLaps(...
    example_lap_data.traversal{1},...
    start_definition,...
    end_definition,...
    excursion_definition,...
    fig_num);

% Do we get 3 laps?
assert(isequal(3,length(lap_traversals.traversal)));

% Are the laps different lengths?
assert(isequal(87,length(lap_traversals.traversal{1}.X)));
assert(isequal(98,length(lap_traversals.traversal{2}.X)));
assert(isequal(79,length(lap_traversals.traversal{3}.X)));

% Plot the lap traversals (should have 3)
fig_num = 3;
fcn_Laps_plotLapsXY(lap_traversals,fig_num);


%% Show how the output works if use full argument list
% This returns the input and exit traversals
fig_num = 4;
[lap_traversals, entry_traversal, exit_traversal] = ...
    fcn_Laps_breakDataIntoLaps(...
    example_lap_data.traversal{1},...
    start_definition,...
    end_definition,...
    excursion_definition,...
    fig_num);

% Do we get 3 laps?
assert(isequal(3,length(lap_traversals.traversal)));

% Are the laps different lengths?
assert(isequal(87,length(lap_traversals.traversal{1}.X)));
assert(isequal(98,length(lap_traversals.traversal{2}.X)));
assert(isequal(79,length(lap_traversals.traversal{3}.X)));

% Do we have an entry_traversal?
assert(isequal(2,length(entry_traversal.X)));

% Do we have an exit_traversal?
assert(isequal(28,length(exit_traversal.X)));

% Plot the results
fig_num = 5;
figure(fig_num);
subplot(1,3,1);
fcn_Laps_plotLapsXY(lap_traversals,fig_num);
axis_limits = axis;
title('All laps');

subplot(1,3,2);
single_lap.traversal{1} = entry_traversal;
fcn_Laps_plotLapsXY(single_lap,fig_num);
axis(axis_limits); % Inheret axis limits from main laps plot
title('Entry segment');

subplot(1,3,3);
single_lap.traversal{1} = exit_traversal;
fcn_Laps_plotLapsXY(single_lap,fig_num);
axis(axis_limits); % Inheret axis limits from main laps plot
title('Exit segment');

%% Show how a lap is missed if start zone is not big enough
start_definition = [6 3 0 0]; % Radius 10, 3 points must pass near [0,0]
end_definition = [30 3 0 -60]; % Radius 30, 3 points must pass near [0,-60]
excursion_definition = []; % empty
fig_num = 6;
lap_traversals = fcn_Laps_breakDataIntoLaps(...
    example_lap_data.traversal{1},...
    start_definition,...
    end_definition,...
    excursion_definition,...
    fig_num);

% Do we get 2 laps?
assert(isequal(2,length(lap_traversals.traversal)));

% Are the laps different lengths?
assert(isequal(86,length(lap_traversals.traversal{1}.X)));
assert(isequal(78,length(lap_traversals.traversal{2}.X)));

% Plot the lap traversals (should have 2)
fig_num = 7;
fcn_Laps_plotLapsXY(lap_traversals,fig_num);

%% Show the use of segment definition
start_definition = [10 0; -10 0]; % start at [10 0], end at [-10 0]
end_definition = [30 3 0 -60]; % Radius 30, 3 points must pass near [0,-60]
excursion_definition = []; % empty
fig_num = 6;
lap_traversals = fcn_Laps_breakDataIntoLaps(...
    example_lap_data.traversal{1},...
    start_definition,...
    end_definition,...
    excursion_definition,...
    fig_num);

% Do we get 3 laps?
assert(isequal(3,length(lap_traversals.traversal)));

% Are the laps different lengths?
assert(isequal(86,length(lap_traversals.traversal{1}.X)));
assert(isequal(97,length(lap_traversals.traversal{2}.X)));
assert(isequal(78,length(lap_traversals.traversal{3}.X)));

% Plot the lap traversals (should have 2)
fig_num = 7;
fcn_Laps_plotLapsXY(lap_traversals,fig_num);



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

%% This one returns nothing since there is no portion of the path in the
% criteria
traversal = fcn_Path_convertPathToTraversalStructure([-1 1; 1 1]);
start_definition = [0.2 3 0 0]; % Located at [0,0] with radius 0.2, 3 points

[lap_traversals, entry_traversal,exit_traversal] = fcn_Laps_breakDataIntoLaps(...
    traversal,...
    start_definition);

assert(isempty(lap_traversals));
assert(isequal(entry_traversal,traversal));
assert(isempty(exit_traversal));

%% This one returns nothing since there is one point in criteria
traversal = fcn_Path_convertPathToTraversalStructure([-1 1; 0 0; 1 1]);
start_definition = [0.2 3 0 0]; % Located at [0,0] with radius 0.2, 3 points
[lap_traversals, entry_traversal,exit_traversal] = fcn_Laps_breakDataIntoLaps(...
    traversal,...
    start_definition);

assert(isempty(lap_traversals));
assert(isequal(entry_traversal,traversal));
assert(isempty(exit_traversal));

%% This one returns nothing since there is only two points in criteria
traversal = fcn_Path_convertPathToTraversalStructure([-1 1; 0 0; 0.1 0; 1 1]);
start_definition = [0.2 3 0 0]; % Located at [0,0] with radius 0.2, 3 points
[lap_traversals, entry_traversal,exit_traversal] = fcn_Laps_breakDataIntoLaps(...
    traversal,...
    start_definition);

assert(isempty(lap_traversals));
assert(isequal(entry_traversal,traversal));
assert(isempty(exit_traversal));

%% This one returns nothing since the minimum point is at the start
% and so there is no strong minimum inside the zone
traversal = fcn_Path_convertPathToTraversalStructure([-1 1; 0 0; 0.01 0; 0.02 0; 0.03 0; 1 1]);
start_definition = [0.2 3 0 0]; % Located at [0,0] with radius 0.2, 3 points
[lap_traversals, entry_traversal,exit_traversal] = fcn_Laps_breakDataIntoLaps(...
    traversal,...
    start_definition);

assert(isempty(lap_traversals));
assert(isequal(entry_traversal,traversal));
assert(isempty(exit_traversal));

%% This one returns nothing since the minimum point is at the end
% and so there is no strong minimum inside the zone
traversal = fcn_Path_convertPathToTraversalStructure([-1 1; -0.03 0; -0.02 0; -0.01 0; 0 0; 1 1]);
start_definition = [0.2 3 0 0]; % Located at [0,0] with radius 0.2, 3 points
[lap_traversals, entry_traversal,exit_traversal] = fcn_Laps_breakDataIntoLaps(...
    traversal,...
    start_definition);

assert(isempty(lap_traversals));
assert(isequal(entry_traversal,traversal));
assert(isempty(exit_traversal));

%% This one returns nothing since the path doesn't come back to start
% There is no end after the start
fig_num = 123;
traversal = fcn_Path_convertPathToTraversalStructure([-1 1; -0.03 0; -0.02 0; 0 0; 0.1 0; 1 1]);
start_definition = [0.2 3 0 0]; % Located at [0,0] with radius 0.2, 3 points

[lap_traversals, entry_traversal,exit_traversal] = fcn_Laps_breakDataIntoLaps(...
    traversal,...
    start_definition,...
    [],...
    [],...
    fig_num);

assert(isempty(lap_traversals));
assert(isequal(entry_traversal,traversal));
assert(isempty(exit_traversal));

%% This one returns nothing since the end is incomplete
fig_num = 123;

% Create some data to plot
full_steps = (-1:0.1:1)';
zero_full_steps = 0*full_steps;
half_steps = (-1:0.1:0)';
zero_half_steps = 0*half_steps;

traversal = fcn_Path_convertPathToTraversalStructure(...
    [full_steps zero_full_steps; zero_half_steps half_steps]);

start_definition = [0.2 3 0 0]; % Located at [0,0] with radius 0.2, 3 points

[lap_traversals, entry_traversal,exit_traversal] = fcn_Laps_breakDataIntoLaps(...
    traversal,...
    start_definition,...
    [],...
    [],...
    fig_num);

assert(isempty(lap_traversals));
assert(isequal(entry_traversal,traversal));
assert(isempty(exit_traversal));

%% Show that the start and end points can overlap by their boundaries
fig_num = 1234;

traversal = fcn_Path_convertPathToTraversalStructure(...
    [full_steps zero_full_steps]);
start_definition = [0.5 3 -0.5 0]; % Located at [-0.5,0] with radius 0.5, 3 points
end_definition = [0.5 3 0.5 0]; % Located at [0.5,0] with radius 0.5, 3 points

[lap_traversals, entry_traversal,exit_traversal] = fcn_Laps_breakDataIntoLaps(...
    traversal,...
    start_definition,...
    end_definition,...
    [],...
    fig_num);

% Do we get 1 laps?
assert(isequal(1,length(lap_traversals.traversal)));


%% Show that the start and end points can be at the absolute ends
fig_num = 1234;

traversal = fcn_Path_convertPathToTraversalStructure(...
    [full_steps zero_full_steps]);

start_definition = [0.5 3 -1 0]; % Located at [-1,0] with radius 0.5, 3 points
end_definition = [0.5 3 1 0]; % Located at [1,0] with radius 0.5, 3 points

[lap_traversals, entry_traversal,exit_traversal] = fcn_Laps_breakDataIntoLaps(...
    traversal,...
    start_definition,...
    end_definition,...
    [],...
    fig_num);

% Do we get 1 laps?
assert(isequal(1,length(lap_traversals.traversal)));


%% Fail conditions
if 1==0
    
    %% Fails because start_definition is not correct type
    clc
    start_definition = [1 2];
    [lap_traversals, input_and_exit_traversals] = fcn_Laps_breakDataIntoLaps(...
        single_lap.traversal{1},...
        start_definition); %#ok<*ASGLU>
    
    %% Fails because start_definition is not correct type
    % Radius input is negative
    clc
    start_definition = [-1 2 3 4];
    [lap_traversals, input_and_exit_traversals] = fcn_Laps_breakDataIntoLaps(...
        single_lap.traversal{1},...
        start_definition); 

    %% Fails because start_definition is not correct type
    % Radius input is negative
    clc
    start_definition = [0 2 3 4];
    [lap_traversals, input_and_exit_traversals] = fcn_Laps_breakDataIntoLaps(...
        single_lap.traversal{1},...
        start_definition); 
    
    %% Fails because start_definition is not correct type
    % Num_inputs input is not positive
    clc
    start_definition = [1 0 3 4];
    [lap_traversals, input_and_exit_traversals] = fcn_Laps_breakDataIntoLaps(...
        single_lap.traversal{1},...
        start_definition); 
    
    %% Warning because start_definition is 3D not 2D
    % Start_zone definition is a 3D point [radius num_points X Y Z]
    clc
    start_definition = [1 2 3 4 5];
    [lap_traversals, input_and_exit_traversals] = fcn_Laps_breakDataIntoLaps(...
        single_lap.traversal{1},...
        start_definition); 
    
    %% Warning because start_definition is 3D not 2D
    % Start_zone definition is a 3D point [X Y Z; X Y Z]
    clc
    start_definition = [1 2 3; 4 5 6];
    [lap_traversals, input_and_exit_traversals] = fcn_Laps_breakDataIntoLaps(...
        single_lap.traversal{1},...
        start_definition); 
    
    %% Warning because end_definition is 3D not 2D
    % End_zone definition is a 3D point [radius num_points X Y Z]
    clc
    start_definition = [1 2 3 4];
    end_definition = [1 2 3 4 5];

    [lap_traversals, input_and_exit_traversals] = fcn_Laps_breakDataIntoLaps(...
        single_lap.traversal{1},...
        start_definition,...
        end_definition); 
    
    %% Warning because excursion_definition is 3D not 2D
    % Excursion_zone definition is a 3D point [radius num_points X Y Z]
    clc
    start_definition = [1 2 3 4];
    end_definition = [1 2 3 4];
    excursion_definition = [1 2 3 4 5];

    [lap_traversals, input_and_exit_traversals] = fcn_Laps_breakDataIntoLaps(...
        single_lap.traversal{1},...
        start_definition,...
        end_definition,...
        excursion_definition); 
end