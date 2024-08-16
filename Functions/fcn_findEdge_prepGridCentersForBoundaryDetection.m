function [X, Y, Z] = fcn_findEdge_prepGridCentersForBoundaryDetection...
    (gridCenters_qualified_grids, gridCenters_unqualified_grids)
%% fcn_findEdge_prepGridCentersForBoundaryDetection
% Prepare the grid centers of qualified and unqualified for boundary detection
% 
% FORMAT:
%
%      [X, Y, Z] = fcn_findEdge_prepGridCentersForBoundaryDetection...
%      (gridCenters_qualified_grids, gridCenters_unqualified_grids, varargin)
%
% INPUTS:     
%       
%      gridCenters_qualified_grids: grids that are qualified, these grids
%      are considered the real boundary points.
%
%      gridCenters_unqualified_grids: grids that are not qualified, these 
%      grids are not considered the real boundary points.
%      
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      X: X meshgrid that contains qualified and unqualified_gridcenters
%
%      Y: Y meshgrid that contains qualified and unqualified_gridcenters
%
%      Z: Z meshgrid that contains qualified and unqualified_gridcenters
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_findEdge_prepGridCentersForBoundaryDetection.m for a full
%       test suite.
%
% This function was written on 2024_08_14 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu
%
% Revision history
% 2024_06_20 - Aneesh Batchu
% -- Wrote the code originally
% 2024_08-14 - Jiabao Zhao
% -- Functionalized this code


%% flag_check_inputs

flag_check_inputs=1;
%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_check_inputs 
    narginchk(2,2); 
end

XYZ_matrix_qualified_grids = [gridCenters_qualified_grids(:,1:2) ones(length(gridCenters_qualified_grids(:,1)),1)]; 

XYZ_matrix_unqualified_grids = [gridCenters_unqualified_grids(:,1:2) zeros(length(gridCenters_unqualified_grids(:,1)),1)]; 

XYZ_matrix_qualified_unqualified_gridcenters = [XYZ_matrix_qualified_grids; XYZ_matrix_unqualified_grids]; 

% Find the unique elements
XYZ_matrix = unique(XYZ_matrix_qualified_unqualified_gridcenters,'rows'); 

% Reverse the matrix
XYZ_matrix_reversed = flipud(XYZ_matrix);

% Find unique rows based on the first two columns (X and Y)
[~,XYZ_matrix_indices] = unique(XYZ_matrix_reversed(:,1:2),'rows'); 

% Create the new matrix with only unique rows, keeping the last occurrence of each duplicate
XYZ_matrix_reversed = XYZ_matrix_reversed(XYZ_matrix_indices,:);  

% Reverse the matrix back to the original order
XYZ_matrix = flipud(XYZ_matrix_reversed);

XYZ_matrix = round(XYZ_matrix,4);

x_range = unique(XYZ_matrix(:,1))'; 
y_range = unique(XYZ_matrix(:,2))';


% Generate X and Y using meshgrid
[X, Y] = meshgrid(x_range, y_range);

% Initialize Z matrix
Z = NaN(size(X));

% Combine X and Y from meshgrid into pairs
XY_pairs = [X(:) Y(:)];

XY_pairs = round(XY_pairs,4); 

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% C = setdiff(XY_pairs, XY_pairs(idx==0,:), 'rows');

% Fill Z matrix with corresponding Z values
Z(idx~=0) = XYZ_matrix(idx(idx~=0), 3);

% Reshape Z to match the dimensions of X and Y
Z = reshape(Z, size(X));

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

% no plot
%%

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§