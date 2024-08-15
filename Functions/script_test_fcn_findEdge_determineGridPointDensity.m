% Script for fcn_findEdge_determineGridPointDensity
%
% REVISION HISTORY:
% 2024_08_14 - Aleksandr Goncharov
% -- Wrote the script
%
% FORMAT:
%
%fcn_findEdge_determineGridPointDensity(total_points_in_each_grid_with_points_greater_than_zero ...
% ,total_points_in_each_grid_in_the_driven_path,grid_size,(N_bins_grid_with_points_greater_than_zero)...
% (N_bins_grid_in_the_driven_path),(fig_num))

%RUNS FUNCTION BASED ON DEMO RESULTS (STILL NEED ADDITIONAL TEST CASES)

%[point_density] = fcn_findEdge_determineGridPointDensity(total_points_in_each_grid_with_points_greater_than_zero,total_points_in_each_grid_in_the_driven_path,grid_size,[],[],33);


%% TEST 1 - SIMPLE CASES - GENERATED DATA - DEFAULT INPUTS

total_points_in_each_grid_with_points_greater_than_zero=randi([0 400],1000,1);
total_points_in_each_grid_in_the_driven_path=randi([0 400],100,1);
grid_size=1;
N_bins_grid_with_points_greater_than_zero=[]; %20 default
N_bins_grid_in_the_driven_path=[]; %10 default
fig_num=1111;


[~] = fcn_findEdge_determineGridPointDensity( ...
    total_points_in_each_grid_with_points_greater_than_zero, ...
    total_points_in_each_grid_in_the_driven_path,grid_size, ...
    N_bins_grid_with_points_greater_than_zero,N_bins_grid_in_the_driven_path, fig_num);

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
%% TEST 2 - SIMPLE CASES - GENERATED DATA - N_bins changed

total_points_in_each_grid_with_points_greater_than_zero=randi([0 400],1000,1);
total_points_in_each_grid_in_the_driven_path=randi([0 400],100,1);
grid_size=1;
N_bins_grid_with_points_greater_than_zero=40; %20 default
N_bins_grid_in_the_driven_path=20; %10 default
fig_num=2222;


[~] = fcn_findEdge_determineGridPointDensity( ...
    total_points_in_each_grid_with_points_greater_than_zero, ...
    total_points_in_each_grid_in_the_driven_path,grid_size, ...
    N_bins_grid_with_points_greater_than_zero,N_bins_grid_in_the_driven_path, fig_num);

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% TEST 3 - SIMPLE CASES - GENERATED DATA - N_bins changed - Changed Grid Size

total_points_in_each_grid_with_points_greater_than_zero=randi([0 400],1000,1);
total_points_in_each_grid_in_the_driven_path=randi([0 400],100,1);
grid_size=0.5; % 1 is default
N_bins_grid_with_points_greater_than_zero=40; %20 default
N_bins_grid_in_the_driven_path=20; %10 default
fig_num=3333;


[~] = fcn_findEdge_determineGridPointDensity( ...
    total_points_in_each_grid_with_points_greater_than_zero, ...
    total_points_in_each_grid_in_the_driven_path,grid_size, ...
    N_bins_grid_with_points_greater_than_zero,N_bins_grid_in_the_driven_path, fig_num);

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% TEST 4 - SPEED MODE


total_points_in_each_grid_with_points_greater_than_zero=randi([0 400],1000,1);
total_points_in_each_grid_in_the_driven_path=randi([0 400],100,1);
grid_size=0.5; % 1 is default
N_bins_grid_with_points_greater_than_zero=40; %20 default
N_bins_grid_in_the_driven_path=20; %10 default
%fig_num=4444;

% Speed Test Calculation
fig_num=[];
REPS=5; minTimeSlow=Inf;
tic;

%slow mode calculation - code copied from plotVehicleXYZ

for i=1:REPS
    tstart=tic;
    
[~] = fcn_findEdge_determineGridPointDensity( ...
    total_points_in_each_grid_with_points_greater_than_zero, ...
    total_points_in_each_grid_in_the_driven_path,grid_size, ...
    N_bins_grid_with_points_greater_than_zero,N_bins_grid_in_the_driven_path, fig_num);
    
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

 [~] = fcn_findEdge_determineGridPointDensity( ...
    total_points_in_each_grid_with_points_greater_than_zero, ...
    total_points_in_each_grid_in_the_driven_path,grid_size, ...
    N_bins_grid_with_points_greater_than_zero,N_bins_grid_in_the_driven_path, fig_num);
   
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

