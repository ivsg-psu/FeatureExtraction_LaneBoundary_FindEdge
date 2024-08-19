function [Xcoord_gridCenters, Ycoord_gridCenters, Zcoord_gridCenters] = fcn_findEdge_prepGridCentersForBoundaryDetection...
    (gridCenters_qualified_grids, gridCenters_unqualified_grids, varargin)
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
%      Xcoord_gridCenters: X meshgrid that contains qualified and unqualified_gridcenters
%
%      Ycoord_gridCenters: Y meshgrid that contains qualified and unqualified_gridcenters
%
%      Zcoord_gridCenters: Z meshgrid that contains qualified and unqualified_gridcenters
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
% 2024_08_14 - Jiabao Zhao
% -- Functionalized this code
% 2024_08_19 - Aneesh Batchu
% -- Added input checking options
% -- Changed the names of variables

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end
%% check input arguments
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

if 0==flag_max_speed
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(2,3);
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if 3<= nargin && 0==flag_max_speed
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

%% Solve for the Maxs and Mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grid centers of qualified grids
XYZ_matrix_qualified_grids = [gridCenters_qualified_grids(:,1:2) ones(length(gridCenters_qualified_grids(:,1)),1)]; 

% Grid centers if unqualified grids
XYZ_matrix_unqualified_grids = [gridCenters_unqualified_grids(:,1:2) zeros(length(gridCenters_unqualified_grids(:,1)),1)]; 

% Combined grid centers
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

% Round the matrx to the fourth decimal
XYZ_matrix = round(XYZ_matrix,4);

x_range = unique(XYZ_matrix(:,1))'; 
y_range = unique(XYZ_matrix(:,2))';


% Generate X and Y using meshgrid
[Xcoord_gridCenters, Ycoord_gridCenters] = meshgrid(x_range, y_range);

% Initialize Z matrix
Zcoord_gridCenters = NaN(size(Xcoord_gridCenters));

% Combine X and Y from meshgrid into pairs
XY_pairs = [Xcoord_gridCenters(:) Ycoord_gridCenters(:)];

XY_pairs = round(XY_pairs,4); 

% Find the indices of XY pairs in the original XYZ matrix
[~, idx] = ismember(XY_pairs, XYZ_matrix(:, 1:2), 'rows');

% C = setdiff(XY_pairs, XY_pairs(idx==0,:), 'rows');

% Fill Z matrix with corresponding Z values
Zcoord_gridCenters(idx~=0) = XYZ_matrix(idx(idx~=0), 3);

% Reshape Z to match the dimensions of X and Y
Zcoord_gridCenters = reshape(Zcoord_gridCenters, size(Xcoord_gridCenters));

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

if flag_do_plots
    fprintf("No plots here")
end % Ends check if plotting

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end
end

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