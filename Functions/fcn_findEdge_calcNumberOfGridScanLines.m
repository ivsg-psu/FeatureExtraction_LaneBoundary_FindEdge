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