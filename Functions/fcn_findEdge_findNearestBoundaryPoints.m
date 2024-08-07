function [isNearest, nearestBorderIndicies, nearestBorderXY] = fcn_findEdge_findNearestBoundaryPoints(boundaryPointsXY, ...
    drivenPathXY, gridSize, gridBoundaries, varargin)
% Find the nearest boundary points of a road
%
% FORMAT:
%
%      [true_borders] = fcn_findEdge_findNearestBoundaryPoints(boundaryPointsXY,
%      gridCenters_non_drivable_grids, gridCenters_driven_path, (fig _num))
%
% INPUTS:
%
%      boundaryPointsXY: The boundary points between the non-drivable
%      path and the drivable path, as an Nx2 matrix.
%
%      drivenPathXY: Grid points that are drivable for vehicles as an Mx2
%      matrix.
%
%      gridSize: the width of the grid in X, XY, or XYZ. The width
%      specifies the decimation interval of the grid. Note: the grid rounds
%      to integer increments of the grid size. For example, if
%      gridBoundaries in X start at 1 and end at 3, then a gridsize of 2
%      starts the grid at 0 and goes to 2 and then 4, not 1 to 3.
%
%      gridBoundaries: a 1x2, 1x4, or 1x6 vector containing the low and
%      high values of the grid in the format [xlow xhigh ylow yhigh zlow
%      zhigh]
%
%      (OPTIONAL INPUT)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      isNearest: an Nx1 array of flags listing whether each
%      boundaryPointsXY point is on the "nearest" boundary
%
%      nearestBorderIndicies: indicies of the nearest borders relative to
%      the grid created by the user-entered gridSize and gridBoundaries.
%      The indicies are specified in the grid indicies listing. Returns an
%      empty matrix if there are no nearest border indicies.
%
%      nearestBorderXY: the XY locations of the nearest borders, where they
%      land on the grid. Returns an empty matrix if there are no nearest
%      border indicies.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_findEdge_findNearestBoundaryPoints.m
%       test suite.
%
% This function was written on 2024_07_25 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu

% Revision history:
% 2024_07_25 - Jiabao Zhao
% -- wrote the code
% 2024_07_31 - S. Brennan
% -- rewrote entire code
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin== 5 && isequal(varargin{end},-1))
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
    debug_fig_num = 999978;
else
    debug_fig_num = [];
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
        narginchk(4,5);

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
if 0==flag_max_speed && 5<= nargin
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

Npoints = length(boundaryPointsXY(:,1));

%% Calculate driven path in terms of grids
% First, convert the driven path into many points such that any grids
% between the user-entered points will also be "hit". We do this by
% "oversampling" the user-given grid to be absolutely sure we do not miss
% any point in-between. So if the user gives a grid size of 0.5 meters,
% then the re-sampled distance will be 10 times smaller, or 0.05 meters. Of
% course, this causes some points on the grid to be represented repeatedly,
% so we remove repeats as a final step

fully_sampled_driven_path = []; % Initialize the grid as empty
oversampling_distance = gridSize/10; % Use a grid 10 times smaller than user entry
for ith_segment = 1:(length(drivenPathXY(:,1))-1)
    segment_start_point = drivenPathXY(ith_segment,:);
    segment_end_point = drivenPathXY(ith_segment+1,:);
    segment_length = real(sum((segment_end_point - segment_start_point).^2,2).^0.5);
    interval_distances = (0:oversampling_distance:segment_length)';
    x_values_on_segment = interp1([0 segment_length],[segment_start_point(1,1) segment_end_point(1,1)],interval_distances);
    y_values_on_segment = interp1([0 segment_length],[segment_start_point(1,2) segment_end_point(1,2)],interval_distances);
    fully_sampled_driven_path = [fully_sampled_driven_path; x_values_on_segment y_values_on_segment]; %#ok<AGROW>
end

% the last point may be missing due to how the intervals above are
% constructed - the x:increment:y method does not necessarily include the y
% endpoint. So add this endpoint
fully_sampled_driven_path = [fully_sampled_driven_path; drivenPathXY(end,:)];

% Find the grid indicies for the driven path. Some are repeated typically,
% and so we keep just the unique ones, as well as keep only the ones that
% are not empty (nan).
gridIndicesDrivenPath = fcn_findEdge_separatePointsIntoGrids(fully_sampled_driven_path, gridSize, gridBoundaries, (-1));
unique_gridIndicesDrivenPath = unique(gridIndicesDrivenPath);
uniqueNonempty_gridIndicesDrivenPath = unique_gridIndicesDrivenPath(~isnan(unique_gridIndicesDrivenPath));

%% Find the indices of boundary points
[gridIndicesBoundaryPoints,~,gridCenters, nGrids] = fcn_findEdge_separatePointsIntoGrids(boundaryPointsXY, gridSize, gridBoundaries, (-1));

% Plot everything for debugging
if flag_do_debug
    figure(debug_fig_num);
    clf;
    hold on;
    grid on;
    axis(gridBoundaries);

    % Plot the grid locations as empty grey points
    plot(gridCenters(:,1),gridCenters(:,2),'.','Color',[0.8 0.8 0.8], 'MarkerSize',20)

    % Plot the input boundary points
    plot(boundaryPointsXY(:,1),boundaryPointsXY(:,2),'r.', 'MarkerSize',30);

    % Plot the raw driven path
    plot(drivenPathXY(:,1),drivenPathXY(:,2),'b.-', 'MarkerSize',10,'LineWidth',3);

    % Plot the driven path converted to nearest indicies
    plot(gridCenters(uniqueNonempty_gridIndicesDrivenPath,1),gridCenters(uniqueNonempty_gridIndicesDrivenPath,2),'.','Color',[0 0 0.9], 'MarkerSize',20)

    % Plot the driven path converted to nearest indicies
    plot(gridCenters(gridIndicesBoundaryPoints,1),gridCenters(gridIndicesBoundaryPoints,2),'.','Color',[0.7 0 0], 'MarkerSize',20)


end

% Remove NaNs from gridIndices if there any
goodGridIndicesBoundaryPoints = gridIndicesBoundaryPoints(~isnan(gridIndicesBoundaryPoints));

% make a zeros matrix that is large enough to fit all the data
z = zeros(nGrids(1),nGrids(2));

% Set the boundary points to 1
z(goodGridIndicesBoundaryPoints) = 1;


%% Find the nearest borders
[drive_path_rows, drive_path_columns] = ind2sub(nGrids', uniqueNonempty_gridIndicesDrivenPath);
drive_path_rows_columns = [drive_path_rows, drive_path_columns];
nearestBorders = fcn_INTERNAL_findNearestBorders(z, drive_path_rows_columns, []);

% Remove repeats
uniqueNearestBorders = unique(nearestBorders,'rows','legacy');

% Convert to indicies
if ~isempty(uniqueNearestBorders)
    nearestBorderIndicies = sub2ind(nGrids',uniqueNearestBorders(:,1),uniqueNearestBorders(:,2));

    % Find XY coordinates
    nearestBorderXY = gridCenters(nearestBorderIndicies,:);

    isNearest = ismember(gridIndicesBoundaryPoints,nearestBorderIndicies);
else
    isNearest = zeros(Npoints,1);
    nearestBorderIndicies = [];
    nearestBorderXY = [];
end


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
    figure(fig_num);
    % temp_h = figure(fig_num);
    % flag_rescale_axis = 0;
    % if isempty(get(temp_h,'Children'))
    %     flag_rescale_axis = 1;
    % end

    % Set up the figure
    clf;
    hold on;
    grid on;
    axis(gridBoundaries);

    legend_texts = {};

    % Plot the grid locations as empty grey points
    plot(gridCenters(:,1),gridCenters(:,2),'.','Color',[0.8 0.8 0.8], 'MarkerSize',20);    
    legend_texts{end+1} = 'Gridcenters';

    %%%
    % Plot the inputs
    % Plot the input boundary points
    plot(boundaryPointsXY(:,1),boundaryPointsXY(:,2),'r.', 'MarkerSize',30);
    legend_texts{end+1} = 'BoundaryPointsXY';

    % Plot the raw driven path
    plot(drivenPathXY(:,1),drivenPathXY(:,2),'b.-', 'MarkerSize',20,'LineWidth',3);
    legend_texts{end+1} = 'DrivenPathXY';

    %%%%%
    % Plot the converted inputs
    % Plot the driven path converted to nearest indicies
    if ~isempty(uniqueNonempty_gridIndicesDrivenPath)
        plot(gridCenters(uniqueNonempty_gridIndicesDrivenPath,1),gridCenters(uniqueNonempty_gridIndicesDrivenPath,2),'o','Color',[0 0 1], 'MarkerSize',10);
        legend_texts{end+1} = 'DrivenPathOnGrid';
    end

    % Plot the driven path converted to nearest indicies
    if ~isempty(gridIndicesBoundaryPoints)
        plot(gridCenters(gridIndicesBoundaryPoints,1),gridCenters(gridIndicesBoundaryPoints,2),'.','Color',[0.7 0 0], 'MarkerSize',20);
        legend_texts{end+1} = 'BoundaryPointsOnGrid';
    end

    %%%%%
    % Plot the results
    if ~isempty(nearestBorderXY)
        plot(nearestBorderXY(:,1),nearestBorderXY(:,2),'.','Color',[1 0 1], 'MarkerSize',20);
        plot(boundaryPointsXY(isNearest,1),boundaryPointsXY(isNearest,2),'o','Color',[0 1 1], 'MarkerSize',20,'LineWidth',5);
        legend_texts{end+1} = 'NearestBorderOnGrid';
        legend_texts{end+1} = 'NearestBorderonBoundary';
    end

    legend(legend_texts);

    % % Make axis slightly larger?
    % if flag_rescale_axis
    %     temp = axis;
    %     %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
    %     axis_range_x = temp(2)-temp(1);
    %     axis_range_y = temp(4)-temp(3);
    %     percent_larger = 0.3;
    %     axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    % end
end
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
%% fcn_INTERNAL_findTrueBorders
function nearest_borders = fcn_INTERNAL_findNearestBorders(border_only_test_grid,drive_path_rows_columns, fig_num)
flag_do_debug = 0;
Nrows = length(border_only_test_grid(:,1));
Ncols = length(border_only_test_grid(1,:));
if flag_do_debug && ~isempty(fig_num)
    %%%
    % Plot the inputs
    figure(38383);
    clf;
    hold on;
    % Plot all the empty grids as points. Note: rows are "height" and columns
    % are "width". Need to be VERY careful keeping track of this difference.
    matrix_indicies = reshape(border_only_test_grid,[],1);
    matrix_indicies_to_plot = matrix_indicies==0;

    row_indicies = repmat((1:Nrows)',Ncols,1);
    column_indicies = reshape(repmat((1:Ncols),Nrows,1),[],1);
    plot(column_indicies(matrix_indicies_to_plot),row_indicies(matrix_indicies_to_plot), '.','Color',[0.8 0.8 0.8]);
    xlim([0 Ncols+1])
    ylim([0 Nrows+1])
    xlabel('Columns (X)');
    ylabel('Rows (Y)');
    % Plot the non-empty values as black dots
    matrix_indicies = reshape(border_only_test_grid,[],1);
    matrix_indicies_to_plot = matrix_indicies==1;
    plot(column_indicies(matrix_indicies_to_plot),row_indicies(matrix_indicies_to_plot),'.','Color',[0 0 0],'MarkerSize',20)
    % Plot the results. Again, be careful: columns are X, rows are Y
    plot(drive_path_rows_columns(:,2),drive_path_rows_columns(:,1),'.','Color',[1 0 0],'MarkerSize',20)
end
% Again, need to be very careful here - the "path" is in XY coordinates and
% is plotted as such. So columns are "X" and rows are "Y".
% Loop through all the path points
nearest_borders = [];
if flag_do_debug
    h_test_drive_path = plot(nan,nan,'o','Color',[1 0 0],'MarkerSize',20);
    h_test_expansion  = plot(nan,nan,'.','Color',[0 0 1],'MarkerSize',20);
    %pause_duration = 0.01;
end
for ith_drive_path_point = 1:length(drive_path_rows_columns)
    % For each path point, find row and column of the current point. All
    % the loops will start at this point and search "outward" for
    % true_borders.
    current_point_to_test = drive_path_rows_columns(ith_drive_path_point,:);
    current_row = current_point_to_test(1,1);
    current_column = current_point_to_test(1,2);
    if flag_do_debug
        set(h_test_drive_path,'Xdata',current_column,'Ydata',current_row);
        % pause(pause_duration);
    end
    % Loop up row-wise
    for ith_row = current_row:Nrows
        if flag_do_debug
            set(h_test_expansion,'Xdata',current_column,'Ydata',ith_row);
            % pause(pause_duration);
        end
        if 1==border_only_test_grid(ith_row,current_column)
            nearest_borders = [nearest_borders; ith_row current_column]; %#ok<AGROW>
            break
        end
    end

    if flag_do_debug && ~isempty(nearest_borders)
        plot(nearest_borders(:,2), nearest_borders(:,1),'g.','MarkerSize',10);
    end
    
    % Loop down row-wise
    for ith_row = current_row:(-1):1
        if flag_do_debug
            set(h_test_expansion,'Xdata',current_column,'Ydata',ith_row);
            %pause(pause_duration);
        end
        if 1==border_only_test_grid(ith_row,current_column)
            nearest_borders = [nearest_borders; ith_row current_column]; %#ok<AGROW>
            break
        end
    end
    if flag_do_debug && ~isempty(nearest_borders)
        plot(nearest_borders(:,2), nearest_borders(:,1),'g.','MarkerSize',10);
    end

    % Loop up column-wise
    for jth_column = current_column:Ncols
        if flag_do_debug
            set(h_test_expansion,'Xdata',jth_column,'Ydata',current_row);
            %  pause(pause_duration);
        end
        if 1==border_only_test_grid(current_row, jth_column)
            nearest_borders = [nearest_borders; current_row jth_column]; %#ok<AGROW>
            break
        end
    end

    if flag_do_debug && ~isempty(nearest_borders)
        plot(nearest_borders(:,2), nearest_borders(:,1),'g.','MarkerSize',10);
    end

    % Loop down column-wise
    for jth_column = current_column:(-1):1
        if flag_do_debug
            set(h_test_expansion,'Xdata',jth_column,'Ydata',current_row);
            %pause(pause_duration);
        end
        if 1==border_only_test_grid(current_row, jth_column)
            nearest_borders = [nearest_borders; current_row jth_column]; %#ok<AGROW>
            break
        end
    end

    if flag_do_debug && ~isempty(nearest_borders)
        plot(nearest_borders(:,2), nearest_borders(:,1),'g.','MarkerSize',10);
    end
end
if ~isempty(fig_num) && ~isempty(nearest_borders)
    % Plot the results
    plot(nearest_borders(:,2),nearest_borders(:,1),'o','Color',[0 1 0],'MarkerSize',20)
end
end % Ends fcn_INTERNAL_findTrueBorders
