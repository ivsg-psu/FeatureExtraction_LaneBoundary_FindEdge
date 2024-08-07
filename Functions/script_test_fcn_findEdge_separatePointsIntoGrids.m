% script_test_fcn_findEdge_separatePointsIntoGrids.m
%
% This is a script to exercise the function: fcn_findEdge_separatePointsIntoGrids.m
% This function was written on  2024_01_22 by V. Wagh (vbw5054@psu.edu)

% Revision history:
% 2024_01_22  by V. Wagh
% -- first write of the code
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it

close all;

%% Basic Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ____            _        ______                           _
%  |  _ \          (_)      |  ____|                         | |
%  | |_) | __ _ ___ _  ___  | |__  __  ____ _ _ __ ___  _ __ | | ___
%  |  _ < / _` / __| |/ __| |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \
%  | |_) | (_| \__ \ | (__  | |____ >  < (_| | | | | | | |_) | |  __/
%  |____/ \__,_|___/_|\___| |______/_/\_\__,_|_| |_| |_| .__/|_|\___|
%                                                      | |
%                                                      |_|
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§
%% BASIC example 1.01: X with centered data, and one bad data point
fig_num = 101;
figure(fig_num);
clf;

inputPoints = [0; 1; 3; 4; 5];
gridSize = 2;
gridBoundaries = [0 4]; 
[gridIndices,grid_AABBs,gridCenters,nGrids] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));


assert(isequal(length(gridIndices(:,1)),length(inputPoints(:,1))));
assert(isequal(length(grid_AABBs(:,1)),2));
assert(isequal(length(gridCenters(:,1)),2));
assert(isequal(length(nGrids(:,1)),1));

assert(isequal(nGrids,2));



%% BASIC example 1.02: many X values
fig_num = 102;
figure(fig_num);
clf;

inputPoints = 10*rand(100,1);
gridSize = 1;
gridBoundaries = [3 9]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));


assert(isequal(length(gridIndices(:,1)),length(inputPoints(:,1))));
assert(isequal(length(grid_AABBs(:,1)),6));
assert(isequal(length(gridCenters(:,1)),6));

%% BASIC example 1.90: X with out-of-range point
fig_num = 190;
figure(fig_num);
clf;

inputPoints = [1; 3; 4; 5];
gridSize = 2;
gridBoundaries = [0 4]; 
[gridIndices,grid_AABBs,gridCenters,nGrids] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));


assert(isequal(length(gridIndices(:,1)),length(inputPoints(:,1))));
assert(isequal(length(grid_AABBs(:,1)),2));
assert(isequal(length(gridCenters(:,1)),2));
assert(isequal(length(nGrids(:,1)),1));

assert(isequal(nGrids,2));

%% BASIC example 1.91: X with weird grid boundaries
% Shows that the grid boundaries get rounded to the outside values of the
% integer gridSize spacing.
fig_num = 191;
figure(fig_num);
clf;

inputPoints = [1; 2];
gridSize = 2;
gridBoundaries = [0.2 3.8]; % This should round to [0 4] producing 3 grids
[gridIndices,grid_AABBs,gridCenters,nGrids] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));


assert(isequal(length(gridIndices(:,1)),length(inputPoints(:,1))));
assert(isequal(length(grid_AABBs(:,1)),2));
assert(isequal(length(gridCenters(:,1)),2));
assert(isequal(length(nGrids(:,1)),1));

assert(isequal(nGrids,2));

%% BASIC example 2.01: XY with centered data, and one bad data point
fig_num = 201;
figure(fig_num);
clf;

inputPoints = [0.2 0; 1 1; 3 3; 3 1; 3 5; 1 5; 1 3; 4 6; 6 6];
gridSize = 2;
gridBoundaries = [0 4 0 6]; 
[gridIndices,grid_AABBs,gridCenters, nGrids] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));


assert(isequal(length(gridIndices(:,1)),length(inputPoints(:,1))));
assert(isequal(length(grid_AABBs(:,1)),6));
assert(isequal(length(gridCenters(:,1)),6));

assert(isequal(nGrids,[2; 3]));


%% BASIC example 2.02: XY with lots of data
fig_num = 202;
figure(fig_num);
clf;

inputPoints = 10*rand(300,2);
gridSize = 2;
gridBoundaries = [2 6 2 8]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

assert(isequal(length(gridIndices(:,1)),length(inputPoints(:,1))));
assert(isequal(length(grid_AABBs(:,1)),6));
assert(isequal(length(gridCenters(:,1)),6));

%% BASIC example 2.03: XY with lots of data and lots of grids
fig_num = 203;
figure(fig_num);
clf;

N_points = 10000;
inputPoints = ones(N_points,1)*[60 110].*rand(N_points,2) - 5*ones(N_points,2);
gridSize = 10;
gridBoundaries = [0 50 0 100]; 

[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

assert(isequal(length(gridIndices(:,1)),length(inputPoints(:,1))));
assert(isequal(length(grid_AABBs(:,1)),50));
assert(isequal(length(gridCenters(:,1)),50));

%% BASIC example 204: XY with bad boundaries
fig_num = 204;
figure(fig_num);
clf;

inputPoints = 10*rand(300,2);
gridSize = 2;
gridBoundaries = [2 2 2 8]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

assert(isequal(length(gridIndices(:,1)),length(inputPoints(:,1))));
assert(isequal(length(grid_AABBs(:,1)),0));
assert(isequal(length(gridCenters(:,1)),0));

%% BASIC example 2.05: XY with bad boundaries
fig_num = 205;
figure(fig_num);
clf;

inputPoints = 10*rand(300,2);
gridSize = 2;
gridBoundaries = [2 -2 2 8]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

assert(isequal(length(gridIndices(:,1)),length(inputPoints(:,1))));
assert(isequal(length(grid_AABBs(:,1)),0));
assert(isequal(length(gridCenters(:,1)),0));

%% BASIC example 2.99: XY with failed test case
fig_num = 299;
figure(fig_num);
clf;

inputPoints = [
    -2.0000    0.2778
    -1.5556    0.2778
    -1.1111    0.7778
    -0.6667    1.2778
    -0.2222    1.7778
    0.2222    2.2778
    0.6667    2.7778
    1.1111    3.2778
    1.5556    3.7778
    2.0000    4.2778
    -1.3333    0.5000
    -0.8889    1.0000
    -0.4444    1.5000
    0.0000    2.0000
    0.4444    2.5000
    0.8889    3.0000
    1.3333    3.5000
    1.7778    4.0000];
gridSize = 0.25;
gridBoundaries = [ -2     2  0     5   ];
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

assert(isequal(length(gridIndices(:,1)),length(inputPoints(:,1))));
assert(isequal(length(grid_AABBs(:,1)),320));
assert(isequal(length(gridCenters(:,1)),320));
   

%% BASIC example 3.01: XYZ with many points
fig_num = 301;
figure(fig_num);
clf;

N_points = 100;
inputPoints = ones(N_points,1)*[4 6 6].*rand(N_points,3);
gridSize = 2;
gridBoundaries = [0 4 0 6 0 6]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

assert(isequal(length(gridIndices(:,1)),length(inputPoints(:,1))));
assert(isequal(length(grid_AABBs(:,1)),18));
assert(isequal(length(gridCenters(:,1)),18));



%% Speed test
N_points = 100;
inputPoints = ones(N_points,1)*[4 6 6].*rand(N_points,3);
gridSize = 2;
gridBoundaries = [0 4 0 6 0 6]; 

% Perform the calculation in slow mode
fig_num = [];
REPS = 1000; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [gridIndices,grid_AABBs,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));

    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fcn_geometry_fillCircleTestPoints without speed setting (slow) and with speed setting (fast):\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Fail conditions
if 1==0
    % FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end

