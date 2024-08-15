% Script for fcn_findEdge_calcNumberOfGridScanLines
%
% REVISION HISTORY:
% 2024_08_15 - Aleksandr Goncharov
% -- Wrote the script
%
% FORMAT:
%
%[total_scan_lines_in_each_grid_with_more_than_zero_points] = fcn_findEdge_calcNumberOfGridScanLines(gridIndices,LIDAR_scanLines,grids_greater_than_zero_points);
%


%% TEST 1 - SIMPLE CASE - GENERATED DATA

gridIndices=randi([400 500],1000,1);
LIDAR_scanLines=linspace(1400,1600,1000)';
grids_greater_than_zero_points=randi([1 100],100,1);

[total_scan_lines_in_each_grid_with_more_than_zero_points] = fcn_findEdge_calcNumberOfGridScanLines(gridIndices,LIDAR_scanLines,grids_greater_than_zero_points);

assert(length(total_scan_lines_in_each_grid_with_more_than_zero_points) == length(grids_greater_than_zero_points))

%% TEST 2 - SIMPLE CASE - GENERATED DATA - VALUES CHANGED

gridIndices=randi([1 500],10005,1);
LIDAR_scanLines=linspace(1200,1400,100)';
grids_greater_than_zero_points=randi([2 412],123,1);

% [total_scan_lines_in_each_grid_with_more_than_zero_points] = fcn_findEdge_calcNumberOfGridScanLines(gridIndices,LIDAR_scanLines,grids_greater_than_zero_points);
% 
% assert(length(total_scan_lines_in_each_grid_with_more_than_zero_points) == length(grids_greater_than_zero_points))