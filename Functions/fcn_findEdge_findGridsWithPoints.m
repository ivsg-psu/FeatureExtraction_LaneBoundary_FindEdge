function [gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, grids_with_zero_points,...
    grids_greater_than_zero_points, gridCenters_zero_point_density,...
    gridCenters_greater_than_zero_point_density, gridIndices, grid_AABBs] = fcn_findEdge_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,varargin)
%% fcn_findEdge_findGridsWithPoints
% divide a set of input points into a grid and analyze the distribution of points within these grids 
% 
%
% FORMAT:
%
%      [gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, grids_with_zero_point,
%      grids_greater_than_zero_points, gridCenters_zero_point_density,
%      gridCenters_greater_than_zero_point_density] = fcn_findEdge_findGridsWithPoints(input_points,
%      grid_size,grid_boundaries,(fig_num))
%
% INPUTS:     
%      
%      (MANDATORY INPUTS)
%       
%      input_points: Lidar data from the outer edge of the road.
%
%      grid_size: Size of grid that we want to seperate the points into 
%      grid size of each grid in 3 dimensions [length, width, height]. 
%
%      grid_boundaries: The boundaries of data in term of X, Y, Z
%      coordnate. For example, max and min of X, Y, Z coordination.
%      
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      gridIndices_cell_array: Cell array that find the repeated number and
%      its indices of each number, and count the total.
%
%      total_N_points_in_each_grid: The number of points in each grid.
%
%      gridCenters: Centers of the grids.
%
%      grids_with_zero_points: Grids that contain zeros points.
%
%      grids_greater_than_zero_points: Grids that contain more than zero
%      points.
%
%      gridCenters_zero_point_density: Centers of the grids with zero
%      point density.
%
%      gridCenters_greater_than_zero_point_density: Centers of the
%      grids with more than zero point density.


% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_findEdge_findGridsWithPoints.m for a full
%       test suite.
%
% Revision History
% 2024_06_15 - Aneesh Batchu
% -- Wrote the code originally
% 2024_06_25 - Aneesh Batchu
% -- NaNs from gridIndices are removed before they are used in the surface
% analysis
% 2024_07_15 - Aneesh Batchu
% -- Seperated this code from fcn_geometry_surfaceAnalysis
% 2024_07_15 - Jiabao Zhao
% -- Functionalized this code
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS");
    MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG = getenv("MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG); 
        flag_check_inputs  = str2double(MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS);
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
        narginchk(1,5);

        % % Check the points input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points, '2column_of_numbers',[2 3]);
        %
        % % Check the transverse_tolerance input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);
        %
        % % Check the station_tolerance input is a positive single number
        % if ~isempty(station_tolerance)
        %     fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);
        % end
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if 2<= nargin && 0==flag_max_speed
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

% Divides the data into grids
[gridIndices,grid_AABBs,gridCenters] = fcn_findEdge_separatePointsIntoGrids(input_points, grid_size, grid_boundaries, (-1));

% Remove NaNs from gridIndices if there any
gridIndices = gridIndices(~isnan(gridIndices));

% Save all the indices to a cell array. Each cell contains the indices of
% the inputPoints that belong to the grid. 
% The second input length(gridCenters(:,1)) is to assign the size of an array. Eaxh grid has a grid center.  
gridIndices_cell_array = fcn_findEdge_createGridPointIndicesCellArray(gridIndices,length(gridCenters(:,1)),-1);

% Calculate the number of points in each grid
total_N_points_in_each_grid = fcn_findEdge_findRepeatedIndices(gridIndices,length(grid_AABBs(:,1)), -1);

% Find grids with zero points (empty grid numbers)
grids_with_zero_points = find((total_N_points_in_each_grid(:,1) == 0)); 

% Grid Centers of the grids with zero point density 
gridCenters_zero_point_density = gridCenters(grids_with_zero_points,1:2); 

% Find grids with more than zero points
grids_greater_than_zero_points = find((total_N_points_in_each_grid(:,1) > 0)); 

% Grid Centers of the grids with zero point density (Unmapped grid centers)
gridCenters_greater_than_zero_point_density = gridCenters(grids_greater_than_zero_points,1:2);

% Find the point density based on the Histogram analysis
% [point_density,actual_driving_surface_grids_hist,total_grids_hist] = fcn_geometry_findPointDensity(total_points_in_each_grid_in_actual_driving_surface,total_points_in_each_grid,fig_num)

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
    marker_size = 10;
    RGB_triplet = [0.8, 0.8, 0.8];
    legend_option = 1;
    legend_name = 'Grids with zero point density';
    plot_gridCenters_zero_point_density = [gridCenters_zero_point_density, zeros(length(gridCenters_zero_point_density(:,1)),1)]; 
    [~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_zero_point_density,marker_size,RGB_triplet,[],legend_option,legend_name,[],[],[],[],fig_num);

    % plot grid centers
    marker_size = 25;
    RGB_triplet = [0.8, 0.8, 0.8];
    legend_option = 1;
    legend_name = 'Grids greater than zero point density';
    plot_gridCenters_greater_than_zero_point_density = [gridCenters_greater_than_zero_point_density, zeros(length(gridCenters_greater_than_zero_point_density(:,1)),1)]; 
    [~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_greater_than_zero_point_density,marker_size,RGB_triplet,[],legend_option,legend_name,[],[],[],[],fig_num);
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
