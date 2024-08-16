%% script_test_fcn_findEdge_determineSTDInZError
% - this is is a script written to test the function:
% fcn_findEdge_determineSTDInZError.m
% Revision history:
% 2024_08_16 - Aleksandr Goncharov
% -- Wrote the script, generated data does not fully reflect real data
%% TEST 1 - BASIC GENERATED DATA

fig_1=1;
fig_2=2;
fig_3=3;

LiDAR_allPoints = rand(1000, 4);

gridIndices_cell_array = cell(100, 1);
for i = 1:100
    gridIndices_cell_array{i} = randi([1, 1000], randi([50, 100]), 1);
end

original_qualified_grids = 1:50;
gridCenters_qualified_grids = rand(50, 2); 
gridCenters_driven_path = rand(10, 2);
current_qualified_grids = original_qualified_grids;
grid_AABBs = rand(100, 4); 
grid_size = 1;
gridIndices = randi([1, 100], 1000, 1);
current_grid_numbers_of_driven_path = 1:10; 


[input_points,original_mapped_gridIndices_cell,total_mapped_grids,total_points_in_mapped_grids,standard_deviation_in_z,gridlines_mapped_grids,...
    driven_path_grid_indices_in_current_mapped_grids,std_in_z_driven_path,std_in_z_other_mapped_grids,mean_std_in_z_driven_path,mean_std_in_z_not_driven_path,max_std_in_z_not_driven_path]...
    = fcn_findEdge_determineSTDInZError(LiDAR_allPoints,gridIndices_cell_array,original_qualified_grids,gridCenters_qualified_grids,gridCenters_driven_path,... 
    current_qualified_grids,grid_AABBs,grid_size,gridIndices,current_grid_numbers_of_driven_path,fig_1,fig_2,fig_3);

assert(length(total_points_in_mapped_grids)==length(original_mapped_gridIndices_cell) & length(standard_deviation_in_z)==length(driven_path_grid_indices_in_current_mapped_grids))


%% TEST 2 - BASIC GENERATED DATA - Diff Fig_num - MORE DATA Points

fig_1=5;
fig_2=6;
fig_3=7;

LiDAR_allPoints = rand(10000, 4);

gridIndices_cell_array = cell(10000, 1);

for i = 1:100
    gridIndices_cell_array{i} = randi([1, 10000], randi([50, 100]), 1);
end

original_qualified_grids = 1:500;
gridCenters_qualified_grids = rand(500, 2); 
gridCenters_driven_path = rand(100, 2);
current_qualified_grids = original_qualified_grids;
grid_AABBs = rand(1000, 4); 
grid_size = 1;
gridIndices = randi([1, 100], 10000, 1);
current_grid_numbers_of_driven_path = 1:10; 


[~,original_mapped_gridIndices_cell,~,total_points_in_mapped_grids,standard_deviation_in_z,~,...
    driven_path_grid_indices_in_current_mapped_grids,~,~,~,~,~]...
    = fcn_findEdge_determineSTDInZError(LiDAR_allPoints,gridIndices_cell_array,original_qualified_grids,gridCenters_qualified_grids,gridCenters_driven_path,... 
    current_qualified_grids,grid_AABBs,grid_size,gridIndices,current_grid_numbers_of_driven_path,fig_1,fig_2,fig_3);

assert(length(total_points_in_mapped_grids)==length(original_mapped_gridIndices_cell) & length(standard_deviation_in_z)==length(driven_path_grid_indices_in_current_mapped_grids))

%% SPEED MODE

fig_1=[];
fig_2=[];

LiDAR_allPoints = rand(1000, 4);

gridIndices_cell_array = cell(100, 1);
for i = 1:100
    gridIndices_cell_array{i} = randi([1, 1000], randi([50, 100]), 1);
end

original_qualified_grids = 1:50;
gridCenters_qualified_grids = rand(50, 2); 
gridCenters_driven_path = rand(10, 2);
current_qualified_grids = original_qualified_grids;
grid_AABBs = rand(100, 4); 
grid_size = 1;
gridIndices = randi([1, 100], 1000, 1);
current_grid_numbers_of_driven_path = 1:10; 



% Perform the calculation in slow mode
fig_num = [];
REPS = 5; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;


fcn_findEdge_determineSTDInZError(LiDAR_allPoints,gridIndices_cell_array,original_qualified_grids,gridCenters_qualified_grids,gridCenters_driven_path,... 
    current_qualified_grids,grid_AABBs,grid_size,gridIndices,current_grid_numbers_of_driven_path,fig_1,fig_2,fig_num);

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


fcn_findEdge_determineSTDInZError(LiDAR_allPoints,gridIndices_cell_array,original_qualified_grids,gridCenters_qualified_grids,gridCenters_driven_path,... 
    current_qualified_grids,grid_AABBs,grid_size,gridIndices,current_grid_numbers_of_driven_path,fig_1,fig_2,fig_num);



    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fcn_findEdge_plotVehicleXYZ without speed setting (slow) and with speed setting (fast):\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);