function [total_scan_lines_in_each_grid_with_more_than_zero_points] = fcn_findEdge_calcNumberOfGridScanLines(gridIndices,LIDAR_scanLines,grids_greater_than_zero_points)
%% fcn_findEdge_calcNumberOfGridScanLines
%
% FORMAT:
%
%   [total_scan_lines_in_each_grid_with_more_than_zero_points] = fcn_findEdge_calcNumberOfGridScanLines(gridIndices,LIDAR_scanLines,grids_greater_than_zero_points)
%
% INPUTS:   
%    
%   gridIndices
%
%   LIDAR_scanLines
%
%   grids_greater_than_zero_points
%
% OUTPUTS:
%
%   total_scan_lines_in_each_grid_with_more_than_zero_points
%
% DEPENDENCIES:
%
% EXAMPLES:
%   See the script:
%   
%   script_test_fcn_findEdge_calcNumberOfGridScanLines
%
%   for a full test suite
%
% This function was written by Aneesh Batchu
%
% REVISION HISTORY:
% 2024_07_25 - Aneesh Batchu
% -- wrote the code originally
% 2024_08_13 - Aneesh Batchu
% -- used accumarray instead of a for loop 
% 2024_08_13 - Aleksandr Goncharov
% -- Functionalized code from the FindEdge Demo 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the array to store scan line counts for each grid index
scanLines_count_per_grid = accumarray(gridIndices, LIDAR_scanLines(:,1), [], @(x) length(unique(x)));

% Extract scan line counts for grids greater than zero points
total_scan_lines_in_each_grid_with_more_than_zero_points = scanLines_count_per_grid(grids_greater_than_zero_points);


% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
% for ith_grid = 1:length(grids_greater_than_zero_points)
%      %Get all points in this domain and plot them
%     rows_in_domain = gridIndices==grids_greater_than_zero_points(ith_grid);
% 
%      %Find number of LiDAR scan lines in each grid
%     scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));
% 
%     % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORK     S WITHOUT A FOR LOOP
%     total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid]; %#ok<AGROW> % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP
% end


%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(NO PLOTS)
end