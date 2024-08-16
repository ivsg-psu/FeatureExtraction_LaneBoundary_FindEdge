%% script for fcn_findEdge_classifyQualifiedGrids
% REVISION HISTORY:
% Written on 2024_08_16 by Aleksandr Goncharov
%
% FORMAT:
%
% [current_qualified_grids,current_unqualified_grids,original_qualified_grids,gridCenters_qualified_grids,gridCenters_unqualified_grids] = ...
%     fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,grid_indices_with_more_than_transverse_span_threshold,...
%     grids_greater_than_zero_points,gridCenters,[],[],27);

%
%
%% TEST CASE 1: SIMPLE CASE - GENERATED DATA
fig_num=1111;

numPoints=1000;
numPointsAll=2000;
grid_indices_with_required_point_density=logical(randi([0 1],1,numPoints)');
grid_indices_with_more_than_one_scan_line=logical(randi([0 1],1,numPoints)');
grid_indices_with_more_than_transverse_span_threshold=logical(randi([0 1],1,numPoints)');
format_unqualified=[];
format_qualified=[];

grids_greater_than_zero_points=linspace(1,numPoints,numPoints);
gridCenters=[linspace(350,450,numPointsAll); linspace(200,250,numPointsAll)]';

[current_qualified_grids,current_unqualified_grids,~,gridCenters_qualified_grids,gridCenters_unqualified_grids] = ...
    fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,grid_indices_with_more_than_transverse_span_threshold,...
    grids_greater_than_zero_points,gridCenters,format_unqualified,format_qualified,fig_num);

assert((length(current_qualified_grids)+length(current_unqualified_grids)) == (length(gridCenters_unqualified_grids)+length(gridCenters_qualified_grids)))

%% TEST CASE 2: SIMPLE CASE - GENERATED DATA - FORMAT CHANGED
fig_num=2222;

numPoints=1000;
numPointsAll=2000;
grid_indices_with_required_point_density=logical(randi([0 1],1,numPoints)');
grid_indices_with_more_than_one_scan_line=logical(randi([0 1],1,numPoints)');
grid_indices_with_more_than_transverse_span_threshold=logical(randi([0 1],1,numPoints)');
format_unqualified=sprintf(' ''.'',''Color'',[1, 0, 0],''MarkerSize'', 25,''DisplayName'',''Unqualified Grids''); (legend ');
format_qualified=sprintf(' ''.'',''Color'',[0, 1, 0],''MarkerSize'', 25,''DisplayName'',''Qualified grids''); (legend ');

grids_greater_than_zero_points=linspace(1,numPoints,numPoints);
gridCenters=[linspace(350,450,numPointsAll); linspace(200,250,numPointsAll)]';

[current_qualified_grids,current_unqualified_grids,~,gridCenters_qualified_grids,gridCenters_unqualified_grids] = ...
    fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,grid_indices_with_more_than_transverse_span_threshold,...
    grids_greater_than_zero_points,gridCenters,format_unqualified,format_qualified,fig_num);

assert((length(current_qualified_grids)+length(current_unqualified_grids)) == (length(gridCenters_unqualified_grids)+length(gridCenters_qualified_grids)))


%% TEST CASE 3 - SPEED TEST

numPoints=1000;
numPointsAll=2000;
grid_indices_with_required_point_density=logical(randi([0 1],1,numPoints)');
grid_indices_with_more_than_one_scan_line=logical(randi([0 1],1,numPoints)');
grid_indices_with_more_than_transverse_span_threshold=logical(randi([0 1],1,numPoints)');
format_unqualified=[];
format_qualified=[];

% Speed Test Calculation
fig_num=[];
REPS=5; minTimeSlow=Inf;
tic;

%slow mode calculation - code copied from plotVehicleXYZ

for i=1:REPS
    tstart=tic;
    
[~,~,~,~,~] = ...
    fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,grid_indices_with_more_than_transverse_span_threshold,...
    grids_greater_than_zero_points,gridCenters,format_unqualified,format_qualified,fig_num);
    
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

[~,~,~,~,~] = ...
    fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,grid_indices_with_more_than_transverse_span_threshold,...
    grids_greater_than_zero_points,gridCenters,format_unqualified,format_qualified,fig_num);
   
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