%% SCRIPT OF JUST RUNNING THE FUNCTION (NEEDS WORK)


%% TEST 1 -  RANDOM GENERATED DATA - BASIC INPUTS
num_grids = 100;
total_num_grids = 200;

grids_greater_than_zero_points = randi([1, total_num_grids], num_grids, 1);
total_N_points_in_each_grid = randi([0, 20], total_num_grids, 1);
point_density = 10;
gridCenters = rand(total_num_grids, 2) * 100;
format_low_point_density=[];
format_required_density=[];

fig_num=2222;

[grid_indices_with_required_point_density, ~] = fcn_findEdge_classifyGridsBasedOnDensity(grids_greater_than_zero_points,total_N_points_in_each_grid,point_density,gridCenters,format_low_point_density,format_required_density,fig_num);

assert(length(grid_indices_with_required_point_density)==length(grids_greater_than_zero_points))
%% TEST 2 -  DIFFERENT RANDOM GENERATED DATA - WAY MORE POINTS - BASIC INPUT
num_grids = 1000;
total_num_grids = 2000;

grids_greater_than_zero_points = randi([1, total_num_grids], num_grids, 1);
total_N_points_in_each_grid = randi([0, 20], total_num_grids, 1);
point_density = 10;
gridCenters = rand(total_num_grids, 2) * 100;
format_low_point_density=[];
format_required_density=[];

fig_num=3333;

[grid_indices_with_required_point_density, ~] = fcn_findEdge_classifyGridsBasedOnDensity(grids_greater_than_zero_points,total_N_points_in_each_grid,point_density,gridCenters,format_low_point_density,format_required_density,fig_num);

assert(length(grid_indices_with_required_point_density)==length(grids_greater_than_zero_points))
%% TEST 3 -  RANDOM GENERATED DATA - FORMAT CHANGED
num_grids = 100;
total_num_grids = 200;

grids_greater_than_zero_points = randi([1, total_num_grids], num_grids, 1);
total_N_points_in_each_grid = randi([0, 20], total_num_grids, 1);
point_density = 10;
gridCenters = rand(total_num_grids, 2) * 100;
format_low_point_density=sprintf(' ''.'',''Color'',[1, 0.8, 0],''MarkerSize'', 25,''DisplayName'', "Grids with low point density" ');
format_required_density=sprintf(' ''.'',''Color'',[0.2, 1, 1],''MarkerSize'', 25,''DisplayName'', "Grids with required density" ');

fig_num=4444;

[grid_indices_with_required_point_density, ~] = fcn_findEdge_classifyGridsBasedOnDensity(grids_greater_than_zero_points,total_N_points_in_each_grid,point_density,gridCenters,format_low_point_density,format_required_density,fig_num);
assert(length(grid_indices_with_required_point_density)==length(grids_greater_than_zero_points))

%% TEST 4 - SPEED MODE


num_grids = 100;
total_num_grids = 200;

grids_greater_than_zero_points = randi([1, total_num_grids], num_grids, 1);
total_N_points_in_each_grid = randi([0, 20], total_num_grids, 1);
point_density = 10;
gridCenters = rand(total_num_grids, 2) * 100;
format_low_point_density=[];
format_required_density=[];

% Speed Test Calculation
fig_num=[];
REPS=5; minTimeSlow=Inf;
tic;

%slow mode calculation - code copied from plotVehicleXYZ

for i=1:REPS
    tstart=tic;
    
fcn_findEdge_classifyGridsBasedOnDensity(grids_greater_than_zero_points,total_N_points_in_each_grid,point_density,gridCenters,format_low_point_density,format_required_density,fig_num);
    
    telapsed=toc(tstart);
    minTimeSlow=min(telapsed,minTimeSlow);
end
averageTimeSlow=toc/REPS;

%slow mode END

%Fast Mode Calculation
fig_num = -1;
minTimeFast = Inf;
tic;

for i=1:REPS
    tstart = tic;

fcn_findEdge_classifyGridsBasedOnDensity(grids_greater_than_zero_points,total_N_points_in_each_grid,point_density,gridCenters,format_low_point_density,format_required_density,fig_num);
   
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

%Display Console Comparison
if 1==1
fprintf(1,'\n\nComparison of fcn_findEdge_plotLIDARXYZ without speed setting (slow) and with speed setting (fast):\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);
end

%Assertion on averageTime NOTE: Due to the variance, there is a chance that
%the assertion will fail.
assert(averageTimeFast<averageTimeSlow);