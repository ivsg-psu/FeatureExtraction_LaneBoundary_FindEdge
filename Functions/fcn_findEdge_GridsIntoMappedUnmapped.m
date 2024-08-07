function [original_grids_with_low_point_density, original_grids_with_required_point_density, original_grids_with_more_than_one_scan_line, original_grids_with_one_scan_line, ...
    original_mapped_grids, original_unmapped_grids, gridCenters_low_point_density, gridCenters_required_point_density, gridCenters_with_more_than_one_scan_line, ...
gridCenters_with_one_scan_line, gridCenters_mapped_grids, gridCenters_unmapped_grids, current_grids_with_low_point_density, current_grids_with_required_point_density, ...
current_grids_with_more_than_one_scan_line, current_grids_with_one_scan_line, current_mapped_grids, current_unmapped_grids]...
= fcn_findEdge_GridsIntoMappedUnmapped(point_density, total_N_points_in_each_grid, total_scan_lines_in_each_grid_with_more_than_zero_points, ...
grids_greater_than_zero_points, gridCenters, varargin)
%% fcn_findEdge_GridsIntoMappedUnmapped
% This function classifies the grids with more than one point into mapped
% and unmapped.
% 
% FORMAT:
%
%      [original_grids_with_required_point_density, gridCenters_low_point_density, gridCenters_required_point_density,...
%      current_grids_with_low_point_density, current_grids_with_required_point_density]...
%      = fcn_findEdge_GridsIntoMappedUnmapped(point_density, total_N_points_in_each_grid, ...
%      grids_greater_than_zero_points, gridCenters,(fig_num))
%
% INPUTS:     
%
%     (MANDATORY INPUTS)
%
%      point_density: The density of points (this is could be entered
%      manually) 
%
%      total_N_points_in_each_grid: The number of points in each grid.
%
%      grids_greater_than_zero_points: Grids that contain more than zero
%      points.
%
%      gridCenters: Centers of the grids.
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      original_grids_with_required_point_density: The grid numbers belong 
%      to the initial grid indices cell array
%
%      gridCenters_low_point_density: Grid Centers of the grids with low 
%      point density (Unmapped grid centers)
%
%      gridCenters_required_point_density: Grid Centers of the grids with 
%      required point density (Mapped grid centers)
%
%      current_grids_with_low_point_density: Current grid numbers of the 
%      grids with low point density
%
%      current_grids_with_required_point_density: grid numbers of the grids 
%      with respect to grids_greater_than_zero_points
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_findEdge_GridsIntoMappedUnmapped.m
%       test suite.
%
% Revision History
% 2024_06_15 - Aneesh Batchu
% -- Wrote the code originally
% 2024_07_15 - Aneesh Batchu
% -- Seperated this code from fcn_geometry_surfaceAnalysis
% 2024_07_15 - Jiabao Zhao
% -- Functionalize the code
% 2024_07_20 - Aneesah Batchu
% -- Added "original_grids_with_more_than_one_scan_line",
% "original_grids_with_one_scan_line", "original_mapped_grids",
% "original_unmapped_grids" as outputs
% -- Added grid centers as the outputs for above mentioned grids
% -- Current grid numbers were also computed for all classified grids
% -- Added new inputs for "total_scan_lines_in_each_grid_with_more_than_zero_points"
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
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
        narginchk(5,6);

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

% % Does user want to specify fig_num?
% flag_do_plots = 0;
% if 5<= nargin && 0==flag_max_speed
%     temp = varargin{end};
%     if ~isempty(temp)
%         fig_num = temp;
%         flag_do_plots = 1;
%     end
% end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (6<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

%% Find Mapped and Unmapped
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These grids have low point density but not zero point density 
original_grids_with_low_point_density = grids_greater_than_zero_points((total_N_points_in_each_grid(grids_greater_than_zero_points,1) > 0) & (total_N_points_in_each_grid(grids_greater_than_zero_points,1) < point_density)); 

% Find mapped grids 
% Mapped grids. Later these grids are classified into drivable and non
% drivable

% These are grids that contain points more than or equal to point density
original_grids_with_required_point_density = grids_greater_than_zero_points(total_N_points_in_each_grid(grids_greater_than_zero_points,1) >= point_density);
grid_indices_with_required_point_density = (total_N_points_in_each_grid(grids_greater_than_zero_points,1) >= point_density);

% These grids contain more than one scan line
original_grids_with_more_than_one_scan_line = grids_greater_than_zero_points(total_scan_lines_in_each_grid_with_more_than_zero_points>1);
grid_indices_with_more_than_one_scan_line = (total_scan_lines_in_each_grid_with_more_than_zero_points>1);

% These grids contain only one scan line
original_grids_with_one_scan_line = grids_greater_than_zero_points(~grid_indices_with_more_than_one_scan_line);

% --------------------------Mapped grid numbers-----------------------------------
% The grids with more than one scan line and more/equal to point density
% are classified as mapped grids
original_mapped_grids = grids_greater_than_zero_points((grid_indices_with_required_point_density == 1) & (grid_indices_with_more_than_one_scan_line == 1)); 

% --------------------------Unmapped grid numbers-----------------------------------
% Find grids with low point density or grids with only one LiDAR scan
% The grids that are not mapped grids are clasified as unmapped grids
original_unmapped_grids = grids_greater_than_zero_points(~((grid_indices_with_required_point_density == 1) & (grid_indices_with_more_than_one_scan_line == 1))); 

% Grid Centers of the grids with low point density (Unmapped grid centers)
gridCenters_low_point_density = gridCenters(original_grids_with_low_point_density,1:2);  

% Grid Centers of the grids with required point density 
gridCenters_required_point_density = gridCenters(original_grids_with_required_point_density,1:2);  

% Grid centers of the grids with more than zero points and more than one scan line
gridCenters_with_more_than_one_scan_line = gridCenters(original_grids_with_more_than_one_scan_line,1:2);

% Grid centers of the grids with one scan line
gridCenters_with_one_scan_line = gridCenters(original_grids_with_one_scan_line,1:2);

% Grid centers of mapped grids
gridCenters_mapped_grids = gridCenters(original_mapped_grids,1:2); 

% Grid centers of unmapped grids
gridCenters_unmapped_grids = gridCenters(original_unmapped_grids,1:2); 

% Current grid numbers of the grids with low point density
% These are the grid numbers of the grids with respect to grids_greater_than_zero_points
% The numbering starts with "1"

% These are the grid numbers of low point density
current_grids_with_low_point_density = find((total_N_points_in_each_grid(grids_greater_than_zero_points,1) > 0) & (total_N_points_in_each_grid(grids_greater_than_zero_points,1) < point_density)); 

% These are the current grid numbers of required point density
current_grids_with_required_point_density = find(total_N_points_in_each_grid(grids_greater_than_zero_points,1) >= point_density); 

% Current grid numbers of the grids with more than one scan line
current_grids_with_more_than_one_scan_line = find(total_scan_lines_in_each_grid_with_more_than_zero_points>1);

% Current grid numbers of the grids with one scan line
current_grids_with_one_scan_line = find(~(total_scan_lines_in_each_grid_with_more_than_zero_points>1));

% Current mapped grid numbers 
current_mapped_grids = find((grid_indices_with_required_point_density == 1) & (grid_indices_with_more_than_one_scan_line == 1)); 

% Current unmapped grid numbers 
current_unmapped_grids = find(~((grid_indices_with_required_point_density == 1) & (grid_indices_with_more_than_one_scan_line == 1))); 

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
    
    % plot unmapped grid centers
    marker_size = 30;
    RGB_triplet = [0.8, 0.8, 0.8];
    legend_option = 1;
    legend_name = 'Unmapped grids';

    plot_gridCenters_unmapped_grids = [gridCenters_unmapped_grids, zeros(length(gridCenters_unmapped_grids(:,1)),1)];
    [~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_unmapped_grids,marker_size,RGB_triplet,[],legend_option,legend_name,[],[],[],[],fig_num);

    % plot mapped grid centers
    marker_size = 30;
    RGB_triplet = [0.2, 0.2, 0.2];
    legend_option = 1;
    legend_name = 'Mapped grids';
    plot_gridCenters_mapped_grids = [gridCenters_mapped_grids, zeros(length(gridCenters_mapped_grids(:,1)),1)];
    [~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_mapped_grids,marker_size,RGB_triplet,[],legend_option,legend_name,[],[],[],[],fig_num);

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


% function original_mapped_grids = fcn_INTERNAL_findOrthoDistancesWithinGridsMoreThanOneScan(original_grids_with_more_than_one_scan_line)
% 
% 
% 
% end