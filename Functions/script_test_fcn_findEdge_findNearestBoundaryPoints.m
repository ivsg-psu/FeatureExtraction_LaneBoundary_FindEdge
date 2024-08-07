% script_test_findEdge_findNearestBoundaryPoints
% Exercises the function: findEdge_findNearestBoundaryPoints

% Revision history:
% 2024_07_25 - Jiabao Zhao 
% -- wrote the code
% 2024_07_31 - S. Brennan
% -- rewrote entire code
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it


%% Test 1 simple example with straight path
fig_num = 1; 
gridSize = 1; 
gridBoundaries = [0 8 0 4];
boundaryPointsXY = [1.5 1.5; 1.5 2.5; 1.5 3.5; 2.5 1.5; 3.5 1.5; 3.5 2.5; 3.5 3.5; 2.5 3.5];
drivenPathXY = [1 0.5; 6 0.5];

[isNearest, nearestBorderIndicies, nearestBorderXY]  = findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

% Check sizes
assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
assert(isequal(length(nearestBorderXY(1,:)),2));
assert(isequal(length(nearestBorderXY(:,2)),length(nearestBorderIndicies)));

% Check values
assert(isequal(find(isNearest),[1 4 5]')); % These are the boundary point "rows" that are closest
assert(isequal(nearestBorderIndicies,[10    11    12]')); % These are the indicies in the grid array


%% Test 2 simple example with curving path
fig_num = 2; 
gridSize = 1; 
gridBoundaries = [0 8 0 4];
boundaryPointsXY = [1.5 1.5; 1.5 2.5; 1.5 3.5; 2.5 1.5; 3.5 1.5; 3.5 2.5; 3.5 3.5; 2.5 3.5];
drivenPathXY = [1 0.5; 4 0.5; 6 2.2];

[isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

% Check sizes
assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
assert(isequal(length(nearestBorderXY(1,:)),2));
assert(isequal(length(nearestBorderXY(:,2)),length(nearestBorderIndicies)));

% Check values
assert(isequal(find(isNearest),[1 4 5 6]')); % These are the boundary point "rows" that are closest
assert(isequal(nearestBorderIndicies,[10    11    12    20]')); % These are the indicies in the grid array



%% Test 3.1 - a larger array locked to grid
fig_num = 31;
figure(fig_num);
clf;


% Create some data
N_points = 10;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,5,15);
gridSize = x_range(2) - x_range(1); 

[X,Y] = meshgrid(x_range,y_range);

% Create Z data that is same size
Z = zeros(size(X));
% Make all Z data which has XY data above line y = x + 2 equal to 1
Y_line = X + 2;
flag_larger_than = Y>Y_line;
Z(flag_larger_than) = 1;
if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

x_limits = [min(x_range) max(x_range)];  
y_limits = [min(y_range) max(y_range)]; 

% Calculate boundary points
boundaryPointsXY = findEdge_findBoundaryPoints(X,Y,Z,gridSize,x_limits,y_limits,-1);

gridSize = 0.25; 
gridBoundaries = [-2 2 -2 5];
drivenPathXY = [-1 4; 0 4.5];

% Convert boundary points into on-grid points? (optional)
if 1==1
    [gridIndicesBoundaryPoints,~,gridCenters, ~] = findEdge_separatePointsIntoGrids(boundaryPointsXY, gridSize, gridBoundaries, (-1));

    % Plot the driven path converted to nearest indicies
    boundaryPointsXYOnGrid = gridCenters(gridIndicesBoundaryPoints,:);

    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXYOnGrid, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

    % Check sizes
    assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
    assert(isequal(length(nearestBorderXY(1,:)),2));
    assert(isequal(length(nearestBorderXY(:,2)),length(nearestBorderIndicies)));

    % Check values
    assert(isequal(find(isNearest),[ 3     4     5     9    10    12    13]')); % These are the boundary point "rows" that are closest
    assert(isequal(nearestBorderIndicies,[180   181   214   215   248   383   416]')); % These are the indicies in the grid array

else
    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
        drivenPathXY, gridSize, gridBoundaries, fig_num);

    % Check sizes
    assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
    assert(isequal(length(nearestBorderXY(1,:)),2));
    assert(isequal(length(nearestBorderXY(:,2)),length(nearestBorderIndicies)));

    % Check values
    assert(isequal(find(isNearest),[ 3     4     5     9    10    12    13]')); % These are the boundary point "rows" that are closest
    assert(isequal(nearestBorderIndicies,[180   181   214   215   248   383   416]')); % These are the indicies in the grid array

end


%% Test 3.2 - a larger array NOT locked to grid 
% This shows how, if the points are not locked onto the grid prior to the
% function call, that weird results can occur

fig_num = 32;
figure(fig_num);
clf;


% Create some data
N_points = 10;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,5,15);
gridSize = x_range(2) - x_range(1); 

[X,Y] = meshgrid(x_range,y_range);

% Create Z data that is same size
Z = zeros(size(X));
% Make all Z data which has XY data above line y = x + 2 equal to 1
Y_line = X + 2;
flag_larger_than = Y>Y_line;
Z(flag_larger_than) = 1;
if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

x_limits = [min(x_range) max(x_range)];  
y_limits = [min(y_range) max(y_range)]; 

% Calculate boundary points
boundaryPointsXY = findEdge_findBoundaryPoints(X,Y,Z,gridSize,x_limits,y_limits,-1);

gridSize = 0.25; 
gridBoundaries = [-2 2 -2 5];
drivenPathXY = [-1 4; 0 4.5];

% Convert boundary points into on-grid points? (optional)
if 1==1
    [gridIndicesBoundaryPoints,~,gridCenters, ~] = findEdge_separatePointsIntoGrids(boundaryPointsXY, gridSize, gridBoundaries, (-1));

    % Plot the driven path converted to nearest indicies
    boundaryPointsXYOnGrid = gridCenters(gridIndicesBoundaryPoints,:);

    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXYOnGrid, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

    % Check sizes
    assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
    assert(isequal(length(nearestBorderXY(1,:)),2));
    assert(isequal(length(nearestBorderXY(:,2)),length(nearestBorderIndicies)));

    % Check values
    assert(isequal(find(isNearest),[ 3     4     5     9    10    12    13]')); % These are the boundary point "rows" that are closest
    assert(isequal(nearestBorderIndicies,[180   181   214   215   248   383   416]')); % These are the indicies in the grid array

else
    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
        drivenPathXY, gridSize, gridBoundaries, fig_num);

    % Check sizes
    assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
    assert(isequal(length(nearestBorderXY(1,:)),2));
    assert(isequal(length(nearestBorderXY(:,2)),length(nearestBorderIndicies)));

    % Check values
    assert(isequal(find(isNearest),[ 3     4     5     9    10    12    13]')); % These are the boundary point "rows" that are closest
    assert(isequal(nearestBorderIndicies,[180   181   214   215   248   383   416]')); % These are the indicies in the grid array

end



%% Test 4.1 - a driven path that is partially outside the range
fig_num = 4;


% Create some data
N_points = 20;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,5,N_points);
gridSize = x_range(2) - x_range(1); 

[X,Y] = meshgrid(x_range,y_range);

% Create Z data that is same size
Z = zeros(size(X));

% Make all Z data which has XY data below line y = x + 2 equal to 1
Y_line = X + 2;
flag_larger_than = Y<Y_line;
Z(flag_larger_than) = 1;

if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

x_limits = [min(x_range) max(x_range)];  
y_limits = [min(y_range) max(y_range)]; 

% Calculate boundary points
boundaryPointsXY = findEdge_findBoundaryPoints(X,Y,Z,gridSize,x_limits,y_limits,-1);

gridSize = 0.25; 
gridBoundaries = [-2 2 -2 5];
drivenPathXY = [-1 4; 0 6.5];

% Convert boundary points into on-grid points? (optional)
if 1==1
    [gridIndicesBoundaryPoints,~,gridCenters, ~] = findEdge_separatePointsIntoGrids(boundaryPointsXY, gridSize, gridBoundaries, (-1));

    % Plot the driven path converted to nearest indicies
    boundaryPointsXYOnGrid = gridCenters(gridIndicesBoundaryPoints,:);

    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXYOnGrid, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

else
    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
        drivenPathXY, gridSize, gridBoundaries, fig_num);
end

% Check sizes
assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
assert(isequal(length(nearestBorderXY(1,:)),2));
assert(isequal(length(nearestBorderXY(:,2)),length(nearestBorderIndicies)));

% Check values
% assert(isequal(find(isNearest),[ 5     7     8     9    12    13    14    15    16    17    18]')); % These are the boundary point "rows" that are closest
% assert(isequal(nearestBorderIndicies,[53    87   120   121   154   187   188   221   222   255   256]')); % These are the indicies in the grid array

%% Test 4.1 - a driven path that is entirely outside the range
fig_num = 4;


% Create some data
N_points = 20;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,5,N_points);
gridSize = x_range(2) - x_range(1); 

[X,Y] = meshgrid(x_range,y_range);

% Create Z data that is same size
Z = zeros(size(X));

% Make all Z data which has XY data below line y = x + 2 equal to 1
Y_line = X + 2;
flag_larger_than = Y<Y_line;
Z(flag_larger_than) = 1;

if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

x_limits = [min(x_range) max(x_range)];  
y_limits = [min(y_range) max(y_range)]; 

% Calculate boundary points
boundaryPointsXY = findEdge_findBoundaryPoints(X,Y,Z,gridSize,x_limits,y_limits,-1);

gridSize = 0.25; 
gridBoundaries = [-2 2 -2 5];
drivenPathXY = [-1 7; 0 6.5];

% Convert boundary points into on-grid points? (optional)
if 1==1
    [gridIndicesBoundaryPoints,~,gridCenters, ~] = findEdge_separatePointsIntoGrids(boundaryPointsXY, gridSize, gridBoundaries, (-1));

    % Plot the driven path converted to nearest indicies
    boundaryPointsXYOnGrid = gridCenters(gridIndicesBoundaryPoints,:);

    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXYOnGrid, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

else
    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
        drivenPathXY, gridSize, gridBoundaries, fig_num);
end

% Check sizes
assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
assert(isempty(nearestBorderXY));
assert(isempty(nearestBorderIndicies));


% Check values
assert(all(isNearest==0))

%% EXAMPLE 5 - a circular domain
fig_num = 5;
figure(fig_num);
clf;


% Create some data
N_points = 20;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,5,N_points);
gridSize = x_range(2) - x_range(1); 

[X,Y] = meshgrid(x_range,y_range);

% Create Z data that is same size
Z = zeros(size(X));

% Make all Z data which are near to origin equal to 1
distance_squared = X.^2 + Y.^2;
radius = 2;
radius_squared = radius.^2;

flag_larger_than = distance_squared<radius_squared;
Z(flag_larger_than) = 1;

if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

x_limits = [min(x_range) max(x_range)];  
y_limits = [min(y_range) max(y_range)]; 

% Calculate boundary points
boundaryPointsXY = findEdge_findBoundaryPoints(X,Y,Z,gridSize,x_limits,y_limits,-1);

gridSize = 0.25; 
gridBoundaries = [-2 2 -2 5];
drivenPathXY = [-1 4; 0 4.5];

% Convert boundary points into on-grid points? (optional)
if 1==1
    [gridIndicesBoundaryPoints,~,gridCenters, ~] = findEdge_separatePointsIntoGrids(boundaryPointsXY, gridSize, gridBoundaries, (-1));

    % Plot the driven path converted to nearest indicies
    boundaryPointsXYOnGrid = gridCenters(gridIndicesBoundaryPoints,:);

    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXYOnGrid, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

else
    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
        drivenPathXY, gridSize, gridBoundaries, fig_num);
end

% Check sizes
assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
assert(isequal(length(nearestBorderXY(1,:)),2));
assert(isequal(length(nearestBorderXY(:,2)),length(nearestBorderIndicies)));

% Check values
% assert(isequal(find(isNearest),[ 5     7     8     9    12    13    14    15    16    17    18]')); % These are the boundary point "rows" that are closest
% assert(isequal(nearestBorderIndicies,[53    87   120   121   154   187   188   221   222   255   256]')); % These are the indicies in the grid array



%% EXAMPLE 6 - driven path entirely outside of domain of boundary points
fig_num = 6;
figure(fig_num);
clf;


% Create some data
N_points = 20;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,2,N_points);
gridSize = x_range(2) - x_range(1); 

[X,Y] = meshgrid(x_range,y_range);

% Create Z data that is same size
Z = zeros(size(X));

% Make all Z data which are near to origin equal to 1
distance_squared = X.^2 + Y.^2;
radius = 2.5;
radius_squared = radius.^2;

flag_larger_than = find(distance_squared<radius_squared);
Z(flag_larger_than) = 1;

if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

x_limits = [min(x_range) max(x_range)];  
y_limits = [min(y_range) max(y_range)]; 

% Calculate boundary points
boundaryPointsXY = findEdge_findBoundaryPoints(X,Y,Z,gridSize,x_limits,y_limits,-1);

gridSize = 0.25; 
gridBoundaries = [-2 2 -2 5];
drivenPathXY = [-1 4; 0 4.5];

% Convert boundary points into on-grid points? (optional)
if 1==1
    [gridIndicesBoundaryPoints,~,gridCenters, ~] = findEdge_separatePointsIntoGrids(boundaryPointsXY, gridSize, gridBoundaries, (-1));

    % Plot the driven path converted to nearest indicies
    boundaryPointsXYOnGrid = gridCenters(gridIndicesBoundaryPoints,:);

    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXYOnGrid, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

else
    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
        drivenPathXY, gridSize, gridBoundaries, fig_num);
end

% Check sizes
assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
assert(isempty(nearestBorderXY)); 
assert(isempty(nearestBorderIndicies));

% Check values
% assert(isequal(find(isNearest),[ 5     7     8     9    12    13    14    15    16    17    18]')); % These are the boundary point "rows" that are closest
% assert(isequal(nearestBorderIndicies,[53    87   120   121   154   187   188   221   222   255   256]')); % These are the indicies in the grid array



%% Test2 Real Data
fig_num = 1002;
figure(fig_num);
clf;

drivenPathXY = [-112.3080   55.5933
 -112.3080   56.3933
 -111.5080   56.3933
 -112.3080   57.1933
 -111.5080   57.1933
 -111.5080   57.9933
 -110.7080   57.9933
 -111.5080   58.7933
 -110.7080   58.7933
 -110.7080   59.5933
 -109.9080   59.5933
 -110.7080   60.3933
 -109.9080   60.3933
 -109.9080   61.1933
 -109.1080   61.1933
 -109.1080   61.9933
 -108.3080   61.9933
 -109.1080   62.7933
 -108.3080   62.7933
];


boundaryPointsXY = [-113.1080   56.7933
 -112.3080   58.3933
 -111.5080   59.9933
 -110.7080   61.5933
 -109.9080   63.1933
 -109.1080   55.9933
 -108.3080   57.5933
 -107.5080   59.1933
 -106.7080   59.9933
 -109.5080   54.7933
 -109.5080   55.5933
 -108.7080   56.3933
 -108.7080   57.1933
 -107.9080   57.9933
 -107.9080   58.7933
 -107.1080   59.5933
 -106.3080   60.3933
 -106.3080   61.1933
 -113.5080   56.3933
 -112.7080   57.1933
 -112.7080   57.9933
 -111.9080   58.7933
 -111.9080   59.5933
 -111.1080   60.3933
 -111.1080   61.1933
 -110.3080   61.9933
 -110.3080   62.7933];


gridSize = 0.5; 

allPoints = [drivenPathXY; boundaryPointsXY];
gridBoundaries = [min(allPoints(:,1)) max(allPoints(:,1)) min(allPoints(:,2)) max(allPoints(:,2))];

% gridBoundaries = [-2 2 -2 5];
% drivenPathXY = [-1 4; 0 4.5];

% Convert boundary points into on-grid points? (optional)
if 1==0
    [gridIndicesBoundaryPoints,~,gridCenters, ~] = findEdge_separatePointsIntoGrids(boundaryPointsXY, gridSize, gridBoundaries, (-1));

    % Plot the driven path converted to nearest indicies
    boundaryPointsXYOnGrid = gridCenters(gridIndicesBoundaryPoints,:);

    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXYOnGrid, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

else
    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
        drivenPathXY, gridSize, gridBoundaries, fig_num);
end

% Check sizes
assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
assert(isequal(length(nearestBorderXY(1,:)),2));
assert(isequal(length(nearestBorderXY(:,2)),length(nearestBorderIndicies)));

%% Test 3, real data

fig_num = 1003;
figure(fig_num);
clf;


boundaryPointsXY = [
393.6480  242.3227
  394.6480  241.3227
  395.6480  240.3227
  396.6480  241.3227
  397.6480  240.3227
  398.6480  240.3227
  399.6480  239.3227
  400.6480  239.3227
  401.6480  239.3227
  402.6480  238.3227
  403.6480  238.3227
  404.6480  226.3227
  404.6480  237.3227
  405.6480  237.3227
  406.6480  225.3227
  406.6480  237.3227
  407.6480  236.3227
  408.6480  236.3227
  409.6480  235.3227
  410.6480  223.3227
  410.6480  235.3227
  411.6480  222.3227
  411.6480  234.3227
  412.6480  234.3227
  413.6480  221.3227
  413.6480  233.3227
  414.6480  232.3227
  415.6480  232.3227
  416.6480  231.3227
  417.6480  230.3227
  418.6480  230.3227
  419.6480  229.3227
  420.6480  229.3227
  421.6480  228.3227
  422.6480  228.3227
  423.6480  227.3227
  424.6480  226.3227
  425.6480  225.3227
  404.6480  227.3227
  406.6480  226.3227
  410.6480  221.3227
  410.6480  224.3227
  411.6480  221.3227
  411.6480  223.3227
  412.6480  220.3227
  413.6480  220.3227
  413.6480  222.3227
  414.6480  221.3227
  415.6480  221.3227
  416.6480  221.3227
  417.6480  221.3227
  418.6480  218.3227
  414.1480  220.8227
  413.1480  221.8227
  411.1480  222.8227
  410.1480  223.8227
  406.1480  225.8227
  425.1480  225.8227
  404.1480  226.8227
  424.1480  226.8227
  423.1480  227.8227
  421.1480  228.8227
  419.1480  229.8227
  417.1480  230.8227
  416.1480  231.8227
  414.1480  232.8227
  413.1480  233.8227
  411.1480  234.8227
  409.1480  235.8227
  407.1480  236.8227
  404.1480  237.8227
  402.1480  238.8227
  399.1480  239.8227
  395.1480  240.8227
  397.1480  240.8227
  394.1480  241.8227
  419.1480  217.8227
  418.1480  218.8227
  418.1480  219.8227
  412.1480  220.8227
  418.1480  220.8227
  414.1480  221.8227
  412.1480  222.8227
  411.1480  223.8227
  407.1480  225.8227
  405.1480  226.8227
  396.1480  240.8227    
]; 


drivenPathXY = [

  424.6480  221.8227
  423.6480  222.8227
  424.6480  222.8227
  421.6480  223.8227
  422.6480  223.8227
  423.6480  223.8227
  420.6480  224.8227
  421.6480  224.8227
  418.6480  225.8227
  419.6480  225.8227
  420.6480  225.8227
  417.6480  226.8227
  418.6480  226.8227
  419.6480  226.8227
  416.6480  227.8227
  417.6480  227.8227
  414.6480  228.8227
  415.6480  228.8227
  416.6480  228.8227
  412.6480  229.8227
  413.6480  229.8227
  414.6480  229.8227
  411.6480  230.8227
  412.6480  230.8227
  409.6480  231.8227
  410.6480  231.8227
  411.6480  231.8227
  407.6480  232.8227
  408.6480  232.8227
  409.6480  232.8227
  405.6480  233.8227
  406.6480  233.8227
  407.6480  233.8227
  403.6480  234.8227
  404.6480  234.8227
  405.6480  234.8227
  401.6480  235.8227
  402.6480  235.8227
  403.6480  235.8227
  398.6480  236.8227
  399.6480  236.8227
  400.6480  236.8227
  401.6480  236.8227
  395.6480  237.8227
  396.6480  237.8227
  397.6480  237.8227
  398.6480  237.8227
  399.6480  237.8227
  393.6480  238.8227
  394.6480  238.8227
  395.6480  238.8227
  396.6480  238.8227
  ]; 

gridSize = 0.5; 

allPoints = [drivenPathXY; boundaryPointsXY];
gridBoundaries = [min(allPoints(:,1)) max(allPoints(:,1)) min(allPoints(:,2)) max(allPoints(:,2))];

% gridBoundaries = [-2 2 -2 5];
% drivenPathXY = [-1 4; 0 4.5];

% Convert boundary points into on-grid points? (optional)
if 1==0
    [gridIndicesBoundaryPoints,~,gridCenters, ~] = findEdge_separatePointsIntoGrids(boundaryPointsXY, gridSize, gridBoundaries, (-1));

    % Plot the driven path converted to nearest indicies
    boundaryPointsXYOnGrid = gridCenters(gridIndicesBoundaryPoints,:);

    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXYOnGrid, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

else
    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
        drivenPathXY, gridSize, gridBoundaries, fig_num);
end

% Check sizes
assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
assert(isequal(length(nearestBorderXY(1,:)),2));
assert(isequal(length(nearestBorderXY(:,2)),length(nearestBorderIndicies)));

 %% Test 4 real data
fig_num = 1004;
figure(fig_num);
clf;

boundaryPointsXY = [361.0678  246.3227
  362.0678  247.3227
  363.0678  246.3227
  364.0678  246.3227
  365.0678  246.3227
  366.0678  246.3227
  367.0678  246.3227
  368.0678  246.3227
  369.0678  246.3227
  370.0678  246.3227
  371.0678  246.3227
  372.0678  246.3227
  373.0678  247.3227
  374.0678  246.3227
  375.0678  246.3227
  376.0678  246.3227
  377.0678  246.3227
  378.0678  245.3227
  379.0678  246.3227
  380.0678  245.3227
  381.0678  245.3227
  382.0678  245.3227
  383.0678  244.3227
  384.0678  244.3227
  385.0678  244.3227
  386.0678  243.3227
  387.0678  243.3227
  388.0678  243.3227
  389.0678  243.3227
  390.0678  243.3227
  391.0678  243.3227
  392.0678  243.3227
  393.0678  242.3227
  394.0678  242.3227
  395.0678  241.3227
  396.0678  241.3227
  397.0678  240.3227
  398.0678  240.3227
  399.0678  241.3227
  400.0678  240.3227
  401.0678  239.3227
  402.0678  239.3227
  403.0678  238.3227
  404.0678  238.3227
  405.0678  237.3227
  406.0678  237.3227
  407.0678  238.3227
  408.0678  236.3227
  409.0678  235.3227
  410.0678  235.3227
  411.0678  234.3227
  412.0678  234.3227
  413.0678  233.3227
  414.0678  233.3227
  415.0678  232.3227
  416.0678  231.3227
  417.0678  231.3227
  418.0678  230.3227
  419.0678  230.3227
  420.0678  229.3227
  421.0678  228.3227
  422.0678  228.3227
  423.0678  227.3227
  424.0678  226.3227
  425.0678  226.3227
  426.0678  225.3227
  366.0678  247.3227
  375.0678  233.3227
  413.0678  219.3227
  414.0678  219.3227
  415.0678  221.3227
  416.0678  220.3227
  414.5678  219.8227
  414.5678  220.8227
  425.5678  225.8227
  423.5678  226.8227
  422.5678  227.8227
  420.5678  228.8227
  419.5678  229.8227
  417.5678  230.8227
  415.5678  231.8227
  374.5678  232.8227
  414.5678  232.8227
  412.5678  233.8227
  410.5678  234.8227
  408.5678  235.8227
  407.5678  236.8227
  404.5678  237.8227
  407.5678  237.8227
  402.5678  238.8227
  400.5678  239.8227
  396.5678  240.8227
  399.5678  240.8227
  394.5678  241.8227
  392.5678  242.8227
  385.5678  243.8227
  382.5678  244.8227
  377.5678  245.8227
  379.5678  245.8227
  362.5678  246.8227
  373.5678  246.8227
  366.5678  247.8227
  416.5678  217.8227
  416.5678  218.8227
  416.5678  219.8227
  415.5678  220.8227
  375.5678  232.8227
  406.5678  237.8227
  398.5678  240.8227
  378.5678  245.8227
  361.5678  246.8227
  372.5678  246.8227
  365.5678  247.8227];


drivenPathXY = [ 424.0678  221.8227
  423.0678  222.8227
  424.0678  222.8227
  425.0678  222.8227
  422.0678  223.8227
  423.0678  223.8227
  421.0678  224.8227
  422.0678  224.8227
  419.0678  225.8227
  420.0678  225.8227
  421.0678  225.8227
  418.0678  226.8227
  419.0678  226.8227
  416.0678  227.8227
  417.0678  227.8227
  418.0678  227.8227
  414.0678  228.8227
  415.0678  228.8227
  416.0678  228.8227
  413.0678  229.8227
  414.0678  229.8227
  415.0678  229.8227
  411.0678  230.8227
  412.0678  230.8227
  413.0678  230.8227
  409.0678  231.8227
  410.0678  231.8227
  411.0678  231.8227
  407.0678  232.8227
  408.0678  232.8227
  409.0678  232.8227
  410.0678  232.8227
  405.0678  233.8227
  406.0678  233.8227
  407.0678  233.8227
  408.0678  233.8227
  403.0678  234.8227
  404.0678  234.8227
  405.0678  234.8227
  406.0678  234.8227
  401.0678  235.8227
  402.0678  235.8227
  403.0678  235.8227
  404.0678  235.8227
  399.0678  236.8227
  400.0678  236.8227
  401.0678  236.8227
  396.0678  237.8227
  397.0678  237.8227
  398.0678  237.8227
  399.0678  237.8227
  393.0678  238.8227
  394.0678  238.8227
  395.0678  238.8227
  396.0678  238.8227
  389.0678  239.8227
  390.0678  239.8227
  391.0678  239.8227
  392.0678  239.8227
  393.0678  239.8227
  394.0678  239.8227
  386.0678  240.8227
  387.0678  240.8227
  388.0678  240.8227
  389.0678  240.8227
  390.0678  240.8227
  381.0678  241.8227
  382.0678  241.8227
  383.0678  241.8227
  384.0678  241.8227
  385.0678  241.8227
  386.0678  241.8227
  387.0678  241.8227
  374.0678  242.8227
  375.0678  242.8227
  376.0678  242.8227
  377.0678  242.8227
  378.0678  242.8227
  379.0678  242.8227
  380.0678  242.8227
  381.0678  242.8227
  382.0678  242.8227
  362.0678  243.8227
  363.0678  243.8227
  364.0678  243.8227
  365.0678  243.8227
  366.0678  243.8227
  367.0678  243.8227
  368.0678  243.8227
  369.0678  243.8227
  370.0678  243.8227
  371.0678  243.8227
  372.0678  243.8227
  373.0678  243.8227
  374.0678  243.8227
  375.0678  243.8227
  376.0678  243.8227
  377.0678  243.8227];

gridSize = 0.5; 

allPoints = [drivenPathXY; boundaryPointsXY];
gridBoundaries = [min(allPoints(:,1)) max(allPoints(:,1)) min(allPoints(:,2)) max(allPoints(:,2))];

% gridBoundaries = [-2 2 -2 5];
% drivenPathXY = [-1 4; 0 4.5];

% Convert boundary points into on-grid points? (optional)
if 1==0
    [gridIndicesBoundaryPoints,~,gridCenters, ~] = findEdge_separatePointsIntoGrids(boundaryPointsXY, gridSize, gridBoundaries, (-1));

    % Plot the driven path converted to nearest indicies
    boundaryPointsXYOnGrid = gridCenters(gridIndicesBoundaryPoints,:);

    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXYOnGrid, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

else
    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
        drivenPathXY, gridSize, gridBoundaries, fig_num);
end

% Check sizes
assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
assert(isequal(length(nearestBorderXY(1,:)),2));
assert(isequal(length(nearestBorderXY(:,2)),length(nearestBorderIndicies)));

%% Test 5 real data

fig_num = 1005;
figure(fig_num);
clf;

boundaryPointsXY = [384.9287   56.3596
  358.9287   19.3596
  359.9287   20.3596
  360.9287   21.3596
  361.9287   21.3596
  362.9287   22.3596
  363.9287   23.3596
  364.9287   23.3596
  365.9287   24.3596
  366.9287   25.3596
  367.9287   26.3596
  368.9287   26.3596
  369.9287   28.3596
  370.9287   29.3596
  371.9287   29.3596
  372.9287   30.3596
  373.9287   30.3596
  374.9287   31.3596
  375.9287   32.3596
  376.9287   32.3596
  377.9287   33.3596
  378.9287   34.3596
  379.9287   34.3596
  380.9287   35.3596
  381.9287   36.3596
  382.9287   36.3596
  383.9287   39.3596
  384.9287   39.3596
  385.9287   39.3596
  386.9287   39.3596
  387.9287   40.3596
  388.9287   41.3596
  389.9287   42.3596
  390.9287   43.3596
  391.9287   43.3596
  392.9287   44.3596
  393.9287   46.3596
  394.9287   46.3596
  395.9287   48.3596
  396.9287   48.3596
  397.9287   49.3596
  398.9287   50.3596
  399.9287   50.3596
  400.9287   51.3596
  401.9287   52.3596
  402.9287   52.3596
  358.4287   18.8596
  359.4287   19.8596
  360.4287   20.8596
  362.4287   21.8596
  363.4287   22.8596
  365.4287   23.8596
  366.4287   24.8596
  367.4287   25.8596
  369.4287   26.8596
  369.4287   27.8596
  370.4287   28.8596
  372.4287   29.8596
  374.4287   30.8596
  375.4287   31.8596
  377.4287   32.8596
  378.4287   33.8596
  380.4287   34.8596
  381.4287   35.8596
  383.4287   36.8596
  383.4287   37.8596
  383.4287   38.8596
  387.4287   39.8596
  388.4287   40.8596
  389.4287   41.8596
  390.4287   42.8596
  392.4287   43.8596
  393.4287   44.8596
  393.4287   45.8596
  395.4287   46.8596
  395.4287   47.8596
  397.4287   48.8596
  398.4287   49.8596
  400.4287   50.8596
  401.4287   51.8596
  385.4287   56.8596];


drivenPathXY = [357.9287   19.8596
  357.9287   20.8596
  358.9287   20.8596
  358.9287   21.8596
  359.9287   21.8596
  360.9287   21.8596
  359.9287   22.8596
  360.9287   22.8596
  361.9287   22.8596
  361.9287   23.8596
  362.9287   23.8596
  362.9287   24.8596
  363.9287   24.8596
  364.9287   24.8596
  364.9287   25.8596
  365.9287   25.8596
  365.9287   26.8596
  366.9287   26.8596
  366.9287   27.8596
  367.9287   27.8596
  368.9287   27.8596
  367.9287   28.8596
  368.9287   28.8596
  369.9287   28.8596
  369.9287   29.8596
  370.9287   29.8596
  370.9287   30.8596
  371.9287   30.8596
  372.9287   30.8596
  371.9287   31.8596
  372.9287   31.8596
  373.9287   31.8596
  373.9287   32.8596
  374.9287   32.8596
  374.9287   33.8596
  375.9287   33.8596
  375.9287   34.8596
  376.9287   34.8596
  377.9287   34.8596
  377.9287   35.8596
  378.9287   35.8596
  378.9287   36.8596
  379.9287   36.8596
  380.9287   36.8596
  379.9287   37.8596
  380.9287   37.8596
  381.9287   37.8596
  381.9287   38.8596
  382.9287   38.8596
  382.9287   39.8596
  383.9287   39.8596
  383.9287   40.8596
  384.9287   40.8596
  384.9287   41.8596
  385.9287   41.8596
  386.9287   41.8596
  386.9287   42.8596
  387.9287   42.8596
  387.9287   43.8596
  388.9287   43.8596
  388.9287   44.8596
  389.9287   44.8596
  389.9287   45.8596
  390.9287   45.8596
  391.9287   45.8596
  390.9287   46.8596
  391.9287   46.8596
  392.9287   46.8596
  392.9287   47.8596
  393.9287   47.8596
  393.9287   48.8596
  394.9287   48.8596
  394.9287   49.8596
  395.9287   49.8596
  396.9287   49.8596
  396.9287   50.8596
  397.9287   50.8596
  397.9287   51.8596
  398.9287   51.8596
  398.9287   52.8596
  399.9287   52.8596
  400.9287   52.8596
  399.9287   53.8596];


gridSize = 0.5; 

allPoints = [drivenPathXY; boundaryPointsXY];
gridBoundaries = [min(allPoints(:,1)) max(allPoints(:,1)) min(allPoints(:,2)) max(allPoints(:,2))];

% gridBoundaries = [-2 2 -2 5];
% drivenPathXY = [-1 4; 0 4.5];

% Convert boundary points into on-grid points? (optional)
if 1==0
    [gridIndicesBoundaryPoints,~,gridCenters, ~] = findEdge_separatePointsIntoGrids(boundaryPointsXY, gridSize, gridBoundaries, (-1));

    % Plot the driven path converted to nearest indicies
    boundaryPointsXYOnGrid = gridCenters(gridIndicesBoundaryPoints,:);

    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXYOnGrid, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

else
    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
        drivenPathXY, gridSize, gridBoundaries, fig_num);
end

% Check sizes
assert(isequal(length(isNearest(:,1)),length(boundaryPointsXY(:,1))));
assert(isequal(length(nearestBorderXY(1,:)),2));
assert(isequal(length(nearestBorderXY(:,2)),length(nearestBorderIndicies)));

%% Speed test
boundaryPointsXY = [384.9287   56.3596
  358.9287   19.3596
  359.9287   20.3596
  360.9287   21.3596
  361.9287   21.3596
  362.9287   22.3596
  363.9287   23.3596
  364.9287   23.3596
  365.9287   24.3596
  366.9287   25.3596
  367.9287   26.3596
  368.9287   26.3596
  369.9287   28.3596
  370.9287   29.3596
  371.9287   29.3596
  372.9287   30.3596
  373.9287   30.3596
  374.9287   31.3596
  375.9287   32.3596
  376.9287   32.3596
  377.9287   33.3596
  378.9287   34.3596
  379.9287   34.3596
  380.9287   35.3596
  381.9287   36.3596
  382.9287   36.3596
  383.9287   39.3596
  384.9287   39.3596
  385.9287   39.3596
  386.9287   39.3596
  387.9287   40.3596
  388.9287   41.3596
  389.9287   42.3596
  390.9287   43.3596
  391.9287   43.3596
  392.9287   44.3596
  393.9287   46.3596
  394.9287   46.3596
  395.9287   48.3596
  396.9287   48.3596
  397.9287   49.3596
  398.9287   50.3596
  399.9287   50.3596
  400.9287   51.3596
  401.9287   52.3596
  402.9287   52.3596
  358.4287   18.8596
  359.4287   19.8596
  360.4287   20.8596
  362.4287   21.8596
  363.4287   22.8596
  365.4287   23.8596
  366.4287   24.8596
  367.4287   25.8596
  369.4287   26.8596
  369.4287   27.8596
  370.4287   28.8596
  372.4287   29.8596
  374.4287   30.8596
  375.4287   31.8596
  377.4287   32.8596
  378.4287   33.8596
  380.4287   34.8596
  381.4287   35.8596
  383.4287   36.8596
  383.4287   37.8596
  383.4287   38.8596
  387.4287   39.8596
  388.4287   40.8596
  389.4287   41.8596
  390.4287   42.8596
  392.4287   43.8596
  393.4287   44.8596
  393.4287   45.8596
  395.4287   46.8596
  395.4287   47.8596
  397.4287   48.8596
  398.4287   49.8596
  400.4287   50.8596
  401.4287   51.8596
  385.4287   56.8596];


drivenPathXY = [357.9287   19.8596
  357.9287   20.8596
  358.9287   20.8596
  358.9287   21.8596
  359.9287   21.8596
  360.9287   21.8596
  359.9287   22.8596
  360.9287   22.8596
  361.9287   22.8596
  361.9287   23.8596
  362.9287   23.8596
  362.9287   24.8596
  363.9287   24.8596
  364.9287   24.8596
  364.9287   25.8596
  365.9287   25.8596
  365.9287   26.8596
  366.9287   26.8596
  366.9287   27.8596
  367.9287   27.8596
  368.9287   27.8596
  367.9287   28.8596
  368.9287   28.8596
  369.9287   28.8596
  369.9287   29.8596
  370.9287   29.8596
  370.9287   30.8596
  371.9287   30.8596
  372.9287   30.8596
  371.9287   31.8596
  372.9287   31.8596
  373.9287   31.8596
  373.9287   32.8596
  374.9287   32.8596
  374.9287   33.8596
  375.9287   33.8596
  375.9287   34.8596
  376.9287   34.8596
  377.9287   34.8596
  377.9287   35.8596
  378.9287   35.8596
  378.9287   36.8596
  379.9287   36.8596
  380.9287   36.8596
  379.9287   37.8596
  380.9287   37.8596
  381.9287   37.8596
  381.9287   38.8596
  382.9287   38.8596
  382.9287   39.8596
  383.9287   39.8596
  383.9287   40.8596
  384.9287   40.8596
  384.9287   41.8596
  385.9287   41.8596
  386.9287   41.8596
  386.9287   42.8596
  387.9287   42.8596
  387.9287   43.8596
  388.9287   43.8596
  388.9287   44.8596
  389.9287   44.8596
  389.9287   45.8596
  390.9287   45.8596
  391.9287   45.8596
  390.9287   46.8596
  391.9287   46.8596
  392.9287   46.8596
  392.9287   47.8596
  393.9287   47.8596
  393.9287   48.8596
  394.9287   48.8596
  394.9287   49.8596
  395.9287   49.8596
  396.9287   49.8596
  396.9287   50.8596
  397.9287   50.8596
  397.9287   51.8596
  398.9287   51.8596
  398.9287   52.8596
  399.9287   52.8596
  400.9287   52.8596
  399.9287   53.8596];


gridSize = 0.5; 

allPoints = [drivenPathXY; boundaryPointsXY];
gridBoundaries = [min(allPoints(:,1)) max(allPoints(:,1)) min(allPoints(:,2)) max(allPoints(:,2))];

[gridIndicesBoundaryPoints,~,gridCenters, nGrids] = findEdge_separatePointsIntoGrids(boundaryPointsXY, gridSize, gridBoundaries, (-1));

% Plot the driven path converted to nearest indicies
boundaryPointsXYOnGrid = gridCenters(gridIndicesBoundaryPoints,:);


% Perform the calculation in slow mode
fig_num = 999;
figure(fig_num);
clf;

REPS = 5; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  findEdge_findNearestBoundaryPoints(boundaryPointsXYOnGrid, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);

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
    [isNearest, nearestBorderIndicies, nearestBorderXY]  =  FINDEDGE_findNearestBoundaryPoints(boundaryPointsXYOnGrid, ...
     drivenPathXY, gridSize, gridBoundaries, fig_num);


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






