function [gridCenters_driven_path, current_grid_numbers_of_driven_path, total_points_in_each_grid_in_the_driven_path, total_points_in_each_grid_with_points_greater_than_zero]...
    = fcn_findEdge_findDrivenPathGrids(gridCenters_greater_than_zero_point_density, boundary_points_driven_path,...
    grids_greater_than_zero_points, total_N_points_in_each_grid, varargin)
%% fcn_findEdge_findDrivenPathGrids   Find the driven path grids within the grids more than zero points
% 
% FORMAT:
%
%      [total_points_in_each_grid_in_the_driven_path, total_points_in_each_grid_with_points_greater_than_zero]...
%      = fcn_findEdge_findDrivenPathGrids(gridCenters_greater_than_zero_point_density, boundary_points_driven_path,...
%      grids_greater_than_zero_points, (fig_num))
%
% INPUTS:     
%       
%      gridCenters_greater_than_zero_point_density: the centers of the grid
%      with greater than zero point density, this means the grid is not empty. 
%
%      boundary_points_driven_path: boundary points that are near the
%      drivable path.
%
%      grids_greater_than_zero_points: grids with greater than zeros point
%      density, this means the gird is not empty.
%
%      total_N_points_in_each_grid: total numbers of points in each grid
%      
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%     total_points_in_each_grid_in_the_driven_path: the total points in
%     each grid of the drivable path.
%
%     total_points_in_each_grid_with_points_greater_than_zero: the total
%     points in each grid with greater than zero. This means the grid is
%     not empty. 
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_findEdge_findDrivenPathGrids.m for a full
%       test suite.
%
% This function was written on 2024_08_13 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu
%
%%% Revision history

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==10 && isequal(varargin{end},-1))
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
        narginchk(4,10);
    end
end

% Check to see if user passed in a string or color style?
if 5 <= nargin
    input = varargin{1};
    if ~isempty(input)
        format = input;
    end
end

% Check to see if user passed in another string or color style?
if 6 <= nargin
    input = varargin{2};
    if ~isempty(input)
        format2 = input;
    end
end


legend_name = 'Grids greater than zero points';
% Check to see if user passed in a legend name? 
if 7 <= nargin
    input = varargin{3};
    if ~isempty(input)
        legend_name = input;
    end
end


legend_name2 = 'Driven path grids';
% Check to see if user passed in another legend name? 
if 8 <= nargin
    input = varargin{4};
    if ~isempty(input)
        legend_name2 = input;
    end
end


flag_do_plots = 0;
if 9<= nargin && 0==flag_max_speed
    temp = varargin{5};
    if ~isempty(temp)
        ENU_XY_fig_num = temp;
        flag_do_plots = 1;
    end
end

% Does user want to specify another fig_num?
if 10<= nargin && 0==flag_max_speed
    temp = varargin{end};
    if ~isempty(temp)
        ENU_XY_fig_num1 = temp;
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

% -----------------------------NOTE------------------------------
% After finding the grids without anypoints, the grids are completely
% removed from the analysis. Only, grids with greater than zero points were
% analyzed from here. 
% -----------------------------NOTE------------------------------

% Plot all the grids greater than zero point density

% "inpolygon" is used to find the grids within the boundary points 
[in,~] = inpolygon(gridCenters_greater_than_zero_point_density(:,1),gridCenters_greater_than_zero_point_density(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));

% Original grid numbers of driven path
original_grid_numbers_of_driven_path = grids_greater_than_zero_points(in); 

% Current grid numbers in driven path 
current_grid_numbers_of_driven_path = find(in); 

% Total points in each grid in the driven path
total_points_in_each_grid_in_the_driven_path = total_N_points_in_each_grid(original_grid_numbers_of_driven_path); 

% Total points in each grid with points greater than zero
total_points_in_each_grid_with_points_greater_than_zero = total_N_points_in_each_grid(grids_greater_than_zero_points); 

% Grid centers of the driven path
gridCenters_driven_path = [gridCenters_greater_than_zero_point_density(in,1),gridCenters_greater_than_zero_point_density(in,2)];

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
    figure(ENU_XY_fig_num);
    fcn_findEdge_plotVehicleXY(gridCenters_greater_than_zero_point_density(:,1:2),(format),(ENU_XY_fig_num));
    fcn_findEdge_plotVehicleXY(gridCenters_driven_path(:,1:2),(format2),(ENU_XY_fig_num));


    % Plot the grids with
    figure(ENU_XY_fig_num1);

    % plot computed boundary points
    marker_size = 10;
    RGB_triplet = [0 0 0];
    legend_option = 1;
    legend_position = [];
    marker_type = [];
    plot_gridCenters_greater_than_zero_point_density = [gridCenters_greater_than_zero_point_density(:,1:2), zeros(length(gridCenters_greater_than_zero_point_density(:,1)),1)];
    [~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_greater_than_zero_point_density,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],ENU_XY_fig_num1);

    % plot driven path
    marker_size = 25;
    RGB_triplet = [0 1 0];
    legend_option = 1;
    legend_position = [];
    marker_type = [];
    plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
    [~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name2,legend_position,[],[],[],ENU_XY_fig_num1);
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

% % "inpolygon" is used to find the grids within the boundary points 
% [in_qg,on_qg] = inpolygon(gridCenters_qualified_grids(:,1),gridCenters_qualified_grids(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));
% 
% % Original grid numbers of driven path
% original_grid_numbers_of_driven_path = original_qualified_grids(in_qg); 
% 
% % Current grid numbers in driven path 
% current_grid_numbers_of_driven_path = current_qualified_grids(in_qg); %find(in); 
% 
% % % Total points in each grid in the driven path
% total_points_in_each_grid_in_the_driven_path = total_N_points_in_each_grid(original_grid_numbers_of_driven_path); 
% % 
% % % Total points in each grid with points greater than zero
% total_points_in_each_grid_with_points_greater_than_zero = total_N_points_in_each_grid(current_qualified_grids); 
% 
% % Grid centers of the driven path
% gridCenters_driven_path = [gridCenters_qualified_grids(in_qg,1),gridCenters_qualified_grids(in_qg,2)];
% 
% 
% fig_num = 517;
% figure(fig_num); clf;
% 
% hold on
% grid on
% xlabel('X[m]')
% ylabel('Y[m]')
% title('Grid centers and boundary points')
% 
% plot(gridCenters_qualified_grids(:,1), gridCenters_qualified_grids(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);
% % plot(boundary_points_driven_path(:,1), boundary_points_driven_path(:,2), '.', 'MarkerSize',30, 'Color',[0 1 0]); 
% 
% % plot the grids in the driven path
% plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',2) % points strictly inside
% 
% % for ith_text = 1:length(current_qualified_grids(:,1))
% %     current_text = sprintf('%.0d',ith_text);
% %     % Place the text on the grid center
% %     text(gridCenters_qualified_grids(ith_text,1), gridCenters_qualified_grids(ith_text,2),current_text,'Color',[1 1 1],'HorizontalAlignment','center','FontSize', 6, 'FontWeight','bold');
% % end
% 
% 
% % Plot the grids with 
% fig_num = 804; 
% figure(fig_num);clf
% 
% % plot computed boundary points
% marker_size = 10;
% RGB_triplet = [0 0 0]; 
% legend_option = 1;
% legend_name = 'Qualified grids';
% legend_position = [];
% marker_type = [];
% plot_gridCenters_qualified_grids = [gridCenters_qualified_grids(:,1:2), zeros(length(gridCenters_qualified_grids(:,1)),1)];
% [~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_qualified_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
% 
% 
% % plot driven path
% marker_size = 25;
% RGB_triplet = [0 1 0]; 
% legend_option = 1;
% legend_name = 'Driven path grids';
% legend_position = [];
% marker_type = [];
% plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
% [~] = fcn_findEdge_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);