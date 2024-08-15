% FORMAT:
%
%      fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid...
%      ,transverse_span_threshold,grids_greater_than_zero_points, gridCenters, (format_1),(format_2),
%      (fig_num))

%% TEST 1 - Randomly Generated Data

num_grids = 100;
total_num_grids = 200;
grid_size = sqrt(total_num_grids);

grids_greater_than_zero_points = (1:num_grids)';
[LON, LAT] = meshgrid(1:grid_size, 1:grid_size);
gridCenters = [LON(:), LAT(:)] * 10;

transverse_span_threshold = 15;
transverse_span_each_grid = randi([5, 25], num_grids, 1);

format_lesser=[];
format_greater=[];


fig_num=1111;

[grid_indices_with_more_than_transverse_span_threshold, ~] = fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid,transverse_span_threshold,grids_greater_than_zero_points,gridCenters,format_lesser,format_greater,fig_num);

assert(length(grid_indices_with_more_than_transverse_span_threshold)==length(grids_greater_than_zero_points))
%% TEST 2 - Randomly Generated Data - Bigger Data

num_grids = 300;
total_num_grids = 600;
grid_size = sqrt(total_num_grids);

grids_greater_than_zero_points = (1:num_grids)';
[LON, LAT] = meshgrid(1:grid_size, 1:grid_size);
gridCenters = [LON(:), LAT(:)] * 10;

transverse_span_threshold = 15;
transverse_span_each_grid = randi([5, 25], num_grids, 1);

format_lesser=[];
format_greater=[];


fig_num=1111;

[grid_indices_with_more_than_transverse_span_threshold, ~] = fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid,transverse_span_threshold,grids_greater_than_zero_points,gridCenters,format_lesser,format_greater,fig_num);

assert(length(grid_indices_with_more_than_transverse_span_threshold)==length(grids_greater_than_zero_points))

%% TEST 3 - Randomly Generated Data - Bigger Data - Different Format

num_grids = 300;
total_num_grids = 600;
grid_size = sqrt(total_num_grids);

grids_greater_than_zero_points = (1:num_grids)';
[LON, LAT] = meshgrid(1:grid_size, 1:grid_size);
gridCenters = [LON(:), LAT(:)] * 10;

transverse_span_threshold = 15;
transverse_span_each_grid = randi([5, 25], num_grids, 1);

format_lesser=sprintf(' ''.'',''Color'',[0, 1, 0],''MarkerSize'', 25,''DisplayName'', "Grids lesser than minimum transverse span threshold" ');
format_greater=sprintf(' ''.'',''Color'',[1, 0, 0],''MarkerSize'', 25,''DisplayName'', "Grids greater than minimum transverse span threshold" ');


fig_num=3333;

[grid_indices_with_more_than_transverse_span_threshold, ~] = fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid,transverse_span_threshold,grids_greater_than_zero_points,gridCenters,format_lesser,format_greater,fig_num);

assert(length(grid_indices_with_more_than_transverse_span_threshold)==length(grids_greater_than_zero_points))

%% TEST 4 - SPEED MODE


num_grids = 300;
total_num_grids = 600;
grid_size = sqrt(total_num_grids);

grids_greater_than_zero_points = (1:num_grids)';
[LON, LAT] = meshgrid(1:grid_size, 1:grid_size);
gridCenters = [LON(:), LAT(:)] * 10;

transverse_span_threshold = 15;
transverse_span_each_grid = randi([5, 25], num_grids, 1);

format_lesser=[];
format_greater=[];


% Speed Test Calculation
fig_num=[];
REPS=5; minTimeSlow=Inf;
tic;

%slow mode calculation - code copied from plotVehicleXYZ

for i=1:REPS
    tstart=tic;
    
fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid,transverse_span_threshold,grids_greater_than_zero_points,gridCenters,format_lesser,format_greater,fig_num);

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

fcn_findEdge_classifyGridsBasedOnTransverseSpan(transverse_span_each_grid,transverse_span_threshold,grids_greater_than_zero_points,gridCenters,format_lesser,format_greater,fig_num);

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
