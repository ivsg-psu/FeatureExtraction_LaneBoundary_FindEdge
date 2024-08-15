% Script for fcn_findEdge_classifyGridsBasedOnScanLines
%
% REVISION HISTORY:
% 2024_08_14 - Aleksandr Goncharov
% -- Wrote the script
% FORMAT:
%
%      fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points...
%      ,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters, (format_1),(format_2),
%      (fig_num))
%

% fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters,[],[],555)

%% TEST 1 - Randomly Generated Values - Simple Script

num_grids = 100;
total_num_grids = 200;
grid_size = sqrt(total_num_grids);
grids_greater_than_zero_points = (1:num_grids)';
[LONG, LAT] = meshgrid(1:grid_size, 1:grid_size);
gridCenters = [LONG(:), LAT(:)]*20;
format_one=[];
format_greater=[];

fig_num=1111;

total_scan_lines_in_each_grid_with_more_than_zero_points = randi([1, 5], num_grids, 1);

fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters,format_one,format_greater,fig_num);

%% TEST 2 - Different Larger data set - Randomly Generated Values - Simple Script

num_grids = 1000;
total_num_grids = 2000;
grid_size = sqrt(total_num_grids);
grids_greater_than_zero_points = (1:num_grids)';
[LONG, LAT] = meshgrid(1:grid_size, 1:grid_size);
gridCenters = [LONG(:), LAT(:)]*10;
format_one=[];
format_greater=[];

fig_num=2222;

total_scan_lines_in_each_grid_with_more_than_zero_points = randi([1, 5], num_grids, 1);

fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters,format_one,format_greater,fig_num);


%% TEST 3 - Format Changed - Green and Red Points

num_grids = 100;
total_num_grids = 200;
grid_size = sqrt(total_num_grids);
grids_greater_than_zero_points = (1:num_grids)';
[LONG, LAT] = meshgrid(1:grid_size, 1:grid_size);
gridCenters = [LONG(:), LAT(:)]*20;
format_one=sprintf(' ''.'',''Color'',[0, 1, 0],''MarkerSize'', 25,''DisplayName'', "Grids with one scan line" ');
format_greater=sprintf(' ''.'',''Color'',[1, 0, 0],''MarkerSize'', 25,''DisplayName'', "Grids with more than one scan line" ');

fig_num=3333;

total_scan_lines_in_each_grid_with_more_than_zero_points = randi([1, 5], num_grids, 1);

fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters,format_one,format_greater,fig_num);

%% TEST 4 - SPEED MODE


num_grids = 100;
total_num_grids = 200;
grid_size = sqrt(total_num_grids);
grids_greater_than_zero_points = (1:num_grids)';
[LONG, LAT] = meshgrid(1:grid_size, 1:grid_size);
gridCenters = [LONG(:), LAT(:)]*20;
format_one=[];
format_greater=[];

total_scan_lines_in_each_grid_with_more_than_zero_points = randi([1, 5], num_grids, 1);

% Speed Test Calculation
fig_num=[];
REPS=5; minTimeSlow=Inf;
tic;

%slow mode calculation - code copied from plotVehicleXYZ

for i=1:REPS
    tstart=tic;
    
fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters,format_one,format_greater,fig_num);

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

fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters,format_one,format_greater,fig_num);

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
