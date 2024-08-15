function boundary_points = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,varargin)
%% fcn_findEdge_findMaxMinOfXYZ
%
% Finds the boundary points of a matrix Z wherein the boundary points are
% the locations halfway between an occupied element of Z (value of 1) and
% an un-occupied element (value of 0). The function returns the halfway
% points along each row or column. The method is to look at where the
% derivative of Z changes first in rows, then in columns, and defining the
% point just 1/2 grid distance beyond this point as the boundary.
%
% FORMAT:
%
% [boundary_points] = fcn_findEdge_findBoundaryPoints(X,Y,Z,grid_size,(x_limits),(y_limits),(fig_num))
%
% INPUTS:
%
% X: X-coordinates of mapped grid centers in an NxM matrix
%
% Y: Y-coordinates of mapped grid centers in an NxM matrix
%
% Z: 1 or 0 (drivable and non-drivable grids) in an NxM matrix
%
% grid size: the size of the grid, in meters
%
% (OPTIONAL INPUTS)
%
% x_limits: x_limits of the plot (used only for plotting)
%
% y_limits: y_limits of the plot (used only for plotting)
%
% fig_num: a figure number to plot results. If set to -1, skips any
% input checking or debugging, no figures will be generated, and sets
% up code to maximize speed.
%
% OUTPUTS:
%
% boundary_points: the boundary points of the drivable and non-drivable
% grids in 2D
%
%
% DEPENDENCIES:
%
% (none)
%
% EXAMPLES:
%
% See the script:
% script_test_fcn_findEdge_findBoundaryPoints.m for a full
% test suite.
%
% This function was written on 2024_06_21 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu or sbrennan@psu.edu

% Revision History
% 2024_06_18 - S. Brennan
% -- wrote the code originally
% 2024_06_21 - Jiabao Zhao
% -- Functionalized the code
% 2024_06_21 - Aneesh Batchu
% -- changed plotting options
% -- added grid size as one of the inputs
% -- made x_range and y_range to x_limits and y_limits
% -- fixed some comments
% 2024_06_25 - Aneesh Batchu7
% -- wrote a conditional statement to keep "Z_greater_than" array to be a
% column matrix.
% 2024_07_31 - Sean Brennan
% -- fixed some missing comments in the header
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==7 && isequal(varargin{end},-1))
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
        narginchk(4,7);

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

x_limits = [];
if (5<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        x_limits = temp;
    end
end

y_limits = [];
if (6<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        y_limits = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if 7<= nargin && 0==flag_max_speed
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

% Find the X_interval
% x_interval = x_range(2)-x_range(1);
% y_interval = y_range(2)-y_range(1);
x_interval = grid_size;
y_interval = grid_size; 

[boundary_points_falling_y, boundary_points_rising_y] = fcn_INTERNAL_findBoundaryPointsX(X, Y, Z, y_interval);
[boundary_points_falling_x_transpose, boundary_points_rising_x_transpose] = fcn_INTERNAL_findBoundaryPointsX(Y', X', Z', x_interval);
boundary_points_falling_x = fliplr(boundary_points_falling_x_transpose);
boundary_points_rising_x = fliplr(boundary_points_rising_x_transpose);

boundary_points = [boundary_points_falling_y; boundary_points_rising_y; boundary_points_falling_x; boundary_points_rising_x];


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

    % Plot the data in 2D
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    hold on;
    grid on;
    axis equal
    xlabel('X [m]')
    ylabel('Y [m]')

    % Plot the results
    flag_larger_than = Z>0.5;
    plot(X(flag_larger_than),Y(flag_larger_than),'k.','Markersize',50);

    if ~isempty(x_limits) &&  ~isempty(y_limits)
        xlim([min(x_limits) max(x_limits)]);
        ylim([min(y_limits) max(y_limits)]);
    end

    % Make axis slightly bigger than range
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    
    if isempty(boundary_points_falling_y)
        boundary_points_falling_y = zeros(0,2);
    end
    if isempty(boundary_points_rising_y)
        boundary_points_rising_y = zeros(0,2);
    end
    if isempty(boundary_points_falling_x)
        boundary_points_falling_x = zeros(0,2);
    end
    if isempty(boundary_points_rising_x)
        boundary_points_rising_x = zeros(0,2);
    end


    % Plot the results
    % plot(boundary_points_falling_y(:,1),boundary_points_falling_y(:,2),'r.','Markersize',40);
    % plot(boundary_points_rising_y(:,1),boundary_points_rising_y(:,2),'g.','Markersize',40);
    % plot(boundary_points_falling_x(:,1),boundary_points_falling_x(:,2),'c.','Markersize',40);
    % plot(boundary_points_rising_x(:,1),boundary_points_rising_x(:,2),'m.','Markersize',40);
    
    plot(boundary_points(:,1),boundary_points(:,2),'c.','Markersize',40);
    p1 = plot(boundary_points(:,1),boundary_points(:,2),'b.','Markersize',20);

    legend(p1,{'Boundary points'})

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% fcn_INTERNAL_findBoundaryPointsX
function [boundary_points_falling, boundary_points_rising] = fcn_INTERNAL_findBoundaryPointsX(X, Y, Z, y_interval)

flag_do_debug = 0;

Z_greater_than = find(Z>0.5);

if size(Z_greater_than,2) ~= 1
    Z_greater_than = Z_greater_than'; 
end
% Pad the first and last points
N_points = numel(X); % Total number of elements
Z_greater_than_padded = [0; Z_greater_than; N_points+1];

%% Find changes
changes_in_sequence = diff(Z_greater_than_padded);
indicies_falling_edge = changes_in_sequence>1.5; %  Anything greater than 1 is a change. Have to add 1 because we padded indicies above
indicies_rising_edge = find(changes_in_sequence>1.5)+1; % Anything greater than 1 is a change. Have to add 1 because we padded indicies above

indicies_with_rising_x_edge  = Z_greater_than_padded(indicies_rising_edge);
indicies_with_falling_x_edge = Z_greater_than_padded(indicies_falling_edge);

% Clean up any indicies outside of range
indicies_with_rising_x_edge(indicies_with_rising_x_edge<1) = [];
indicies_with_rising_x_edge(indicies_with_rising_x_edge>N_points) = [];
indicies_with_falling_x_edge(indicies_with_falling_x_edge<1) = [];
indicies_with_falling_x_edge(indicies_with_falling_x_edge>N_points) = [];

% Clean up any indicies on borders
border_dimension = length(X(:,1));
border_indices_rising = mod(indicies_with_rising_x_edge,border_dimension)==1;
indicies_with_rising_x_edge(border_indices_rising) = [];
border_indices_falling = mod(indicies_with_falling_x_edge,border_dimension)==0;
indicies_with_falling_x_edge(border_indices_falling) = [];




boundary_points_rising  = [X(indicies_with_rising_x_edge),Y(indicies_with_rising_x_edge)-y_interval/2];
boundary_points_falling = [X(indicies_with_falling_x_edge),Y(indicies_with_falling_x_edge)+y_interval/2];

% Plot the data in 2D?
if 1==flag_do_debug
    figure(1111);
    clf;
    hold on;

    % Plot the inputs
    plot(X(1:N_points),Y(1:N_points),'.','Color',[0.5 0.5 0.5],'Markersize',20);
    plot(X(Z_greater_than),Y(Z_greater_than),'k.','Markersize',50);

    % Number the results (for clarity)
    for ith_label = 1:length(Z_greater_than)
        label_number = Z_greater_than(ith_label);
        current_text = sprintf('%.0d',label_number);
        text(X(label_number),Y(label_number),current_text,'Color',[0.5 0.5 0.5],'HorizontalAlignment','center');
    end
    xlabel('X [m]');
    ylabel('Y [m]');

    % Start by forcing tight axes
    xlim([min(X,[],'all') max(X,[],'all')]);
    ylim([min(Y,[],'all') max(Y,[],'all')]);

    % Make axis slightly bigger than range
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

    plot(boundary_points_falling(:,1),boundary_points_falling(:,2),'r.','Markersize',30);
    plot(boundary_points_rising(:,1),boundary_points_rising(:,2),'g.','Markersize',30);
end
end