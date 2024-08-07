%% script_test_fcn_findEdge_findBoundaryPoints
% - this is is a script written to test the function:
% fcn_findEdge_findBoundaryPoints.m
%

% Revision history:
% 2024_06_18 - S. Brennan
% -- wrote these examples for demonstration
% 2024_06_21 - Jiabao Zhao
% -- organized the script
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it

close all;

%% Simple example 1
fig_num = 1;

x = [1,2];
y = [1,2];
grid_size = x(2) - x(1); 
[X,Y] = meshgrid(x,y);

% draw a boundary line between first row and second row
Z = [1 1; 0 0];

if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

x_limits = [min(x) max(x)];  
y_limits = [min(y) max(y)]; 

% Calculate boundary points
boundary_points = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);
plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

% Plot the boundary line
y_line = [1.5,1.5];
plot(x,y_line,'b-','LineWidth',3);

assert(isequal(boundary_points(1,2),1.5));
assert(isequal(boundary_points(2,2),1.5));

%% Simple example 2

x = [1,2,3];
y = [1,2,3];
grid_size = x(2) - x(1); 
[X,Y] = meshgrid(x,y);

% draw a boundary line between third row and second row

Z = [1 1 1; 1 1 1;0 0 0];

if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

x_limits = [min(x) max(x)];  
y_limits = [min(y) max(y)]; 

% Calculate boundary points
fig_num = 123;
boundary_points = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);
plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

% Plot the boundary line
y_line = [2.5,2.5,2.5];
plot(x,y_line,'b-','LineWidth',3);

assert(isequal(boundary_points(1,2),2.5));
assert(isequal(boundary_points(2,2),2.5));
assert(isequal(boundary_points(3,2),2.5));

%% Simple example 3

x = [1,2,3,4,5];
y = [1,2,3,4,5];
grid_size = x(2) - x(1); 
[X,Y] = meshgrid(x,y);

% draw a boundary lines of a rectangle shape with length of 3

Z = [0 0 0 0 0;0 1 1 1 0;0 1 1 1 0;0 1 1 1 0; 0 0 0 0 0];

if 1==0
    % Plot the data in 3D
    figure(1234);
    clf;
    surf(X,Y,Z)
end

x_limits = [min(x) max(x)];  
y_limits = [min(y) max(y)]; 

% Calculate boundary points
fig_num = 123;
boundary_points = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);
plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

% Plot the boundary line
y_line = [4.5,4.5,4.5,4.5,4.5];
y_line1 = [1.5,1.5,1.5,1.5,1.5];
plot(x,y_line,x,y_line1,y_line,y,y_line1,y,'b-','LineWidth',3);

assert(isequal(boundary_points(1,2),4.5));
assert(isequal(boundary_points(2,2),4.5));
assert(isequal(boundary_points(3,2),4.5));
assert(isequal(boundary_points(4,2),1.5));
assert(isequal(boundary_points(5,2),1.5));
assert(isequal(boundary_points(6,2),1.5));
assert(isequal(boundary_points(7,1),4.5));
assert(isequal(boundary_points(8,1),4.5));
assert(isequal(boundary_points(9,1),4.5));
assert(isequal(boundary_points(10,1),1.5));
assert(isequal(boundary_points(11,1),1.5));
assert(isequal(boundary_points(12,1),1.5));

%% Examples for Aneesh
fig_num = 1;
figure(fig_num);
clf;


% Create some data
N_points = 10;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,5,15);
grid_size = x_range(2) - x_range(1); 

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
fig_num = 123;
boundary_points = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);
plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

% Plot the boundary line
y_line = x_range+2;
plot(x_range,y_line,'b-','LineWidth',3);

% Assert check the length of column of the output  
assert(isequal(length(boundary_points(:,1)),length(boundary_points(:,2))));

%% EXAMPLE 2
fig_num = 2;
figure(fig_num);
clf;

% Create some data
N_points = 20;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,5,N_points);
grid_size = x_range(2) - x_range(1); 

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
fig_num = 124;
boundary_points = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);
plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

% Plot the boundary line
y_line = x_range+2;
plot(x_range,y_line,'b-','LineWidth',3);

% Assert check the length of column of the output  
assert(isequal(length(boundary_points(:,1)),length(boundary_points(:,2))));
%% EXAMPLE 3
fig_num = 3;
figure(fig_num);
clf;


% Create some data
N_points = 20;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,5,N_points);
grid_size = x_range(2) - x_range(1); 

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
boundary_points = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

% Plot the boundary circle
% angles = linspace(0,360,100)'*pi/180;
% plot(radius*cos(angles),radius*sin(angles),'b-','LineWidth',3);

% Assert check the length of column of the output 
assert(isequal(length(boundary_points(:,1)),length(boundary_points(:,2))));

%% EXAMPLE 4
fig_num = 4;
figure(fig_num);
clf;


% Create some data
N_points = 20;
x_range = linspace(-2,2,N_points);
y_range = linspace(-2,2,N_points);
grid_size = x_range(2) - x_range(1); 

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
boundary_points = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num);
plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

% % Plot the boundary circle
% angles = linspace(0,360,100)'*pi/180;
% plot(radius*cos(angles),radius*sin(angles),'b-','LineWidth',3);

% Assert check the length of column of the output 
assert(isequal(length(boundary_points(:,1)),length(boundary_points(:,2))));