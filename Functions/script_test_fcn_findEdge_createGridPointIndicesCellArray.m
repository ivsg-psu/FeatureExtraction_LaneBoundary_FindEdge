%% script_test_fcn_findEdge_createGridPointIndicesCellArray
% Exerrcises the function: fcn_findEdge_createGridPointIndicesCellArray.m
% Revision history:
% 2024_6_17
% Jiabao Zhao wrote the code
% 2024_06_19 - Aneesh Batchu
% -- Added a case when "length_indices_array" is more than max(indices) of
% array
% -- Changed the name of the script


% To-do
% 2024_06_20 - Aneesh Batchu
% -- write assertions
% -- write simple examples. 

close all

%% Test 1 
rng(123)
N_points = 1500; 
true_parameters = [ 0 0 3]';
points = ones(N_points,1)+ 9*rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data
z = true_z;
points1 = [x, y, z]; 
% fig_num = 1011;
% figure(fig_num)
% 
% view(3)
% plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
% 
% grid on
% xlabel('x')
% ylabel('y')
% zlabel('z')
fig_num = 10333; 
figure(fig_num); clf;
inputPoints = points1;
gridSize = 1;
gridBoundaries = [1 10 1 10 2.5 3.5]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_findEdge_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));
output = fcn_findEdge_findRepeatedIndices(gridIndices, length(grid_AABBs));
cell_array = fcn_findEdge_createGridPointIndicesCellArray(gridIndices,length(grid_AABBs), fig_num);
assert((length(cell_array{1}))==21);

%% test 2

% Create data

rng(123)

N_points = 100; 

true_parameters = [ 0 0 3]';
points = ones(N_points,1)+ rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

% Add noise
% true_sigma = 0.1;
% z = true_z + true_sigma*randn(Npoints,1);

points1 = [x, y, z]; 

% Create data

rng(123)

N_points = 100; 

true_parameters = [ 0 0 3]';
points = [2 1].*ones(N_points,1)+ rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

% z = true_z;

% Add noise
true_sigma = 0.15;
z = true_z + true_sigma*randn(N_points,1);

points2 = [x, y, z]; 

true_parameters = [ 1 0 1]';
points = [3 1].*ones(N_points,2) + rand(N_points,2);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

points3 = [x, y, z]; 

fig_num = 1; clf;
% Plot the points
figure(fig_num)
plot3(points1(:,1),points1(:,2),points1(:,3),'.','MarkerSize',20,'Color',[0 0 0]);
hold on
plot3(points2(:,1),points2(:,2),points2(:,3),'.','MarkerSize',20,'Color',[1 0 0]);
hold on
plot3(points3(:,1),points3(:,2),points3(:,3),'.','MarkerSize',20,'Color',[1 0 0]);

grid on
xlabel('x')
ylabel('y')
zlabel('z')

% Surface Analysis

fig_num = 334; 
figure(fig_num); clf;
inputPoints = [points1; points2; points3]; 
gridSize = 1;
gridBoundaries = [1 4 1 2 2.5 5.5]; 
[gridIndices,grid_AABBs,gridCenters] = fcn_findEdge_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries, (fig_num));
output = fcn_findEdge_findRepeatedIndices(gridIndices, length(grid_AABBs));
cell_array = fcn_findEdge_createGridPointIndicesCellArray(gridIndices,length(grid_AABBs), fig_num);
% assert((length(cell_array{1}))==21);


