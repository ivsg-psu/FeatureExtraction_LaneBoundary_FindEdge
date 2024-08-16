function [angle_btw_unit_normals_and_vertical, mean_angle_btw_unit_normals_and_vertical_driven_path,angle_btw_unit_normals_and_vertical_driven_path] = fcn_findEdge_determineAngleDeviation...
    (LiDAR_allPoints, gridIndices_cell_array, original_qualified_grids,...
    gridCenters_qualified_grids,current_qualified_grids,gridCenters_driven_path, ...
    grid_AABBs, grid_size, gridIndices, current_grid_numbers_of_driven_path, varargin)
%% fcn_findEdge_determineAngleDeviation   Qualified grid conditions - angle deviation
% 
% FORMAT:
%
%      [angle_btw_unit_normals_and_vertical, mean_angle_btw_unit_normals_and_vertical_driven_path] = fcn_findEdge_determineAngleDeviation...
%      (LiDAR_allPoints, gridIndices_cell_array, original_qualified_grids,...
%      gridCenters_qualified_grids,current_qualified_grids,gridCenters_driven_path, 
%      (fig_num),(fig_num2),(fig_num3) )
%
% INPUTS:     
%       
%      LiDAR_allPoints: concatenated LIDAR points in the domain, these are 
%      the schosen scans and cells. 
%
%      gridIndices_cell_array: grid indices of each cell array
%
%      original_qualified_grids: qualified grids. These girds passed tests 
%      and will be passed for future analysis.
%      
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      angle_btw_unit_normals_and_vertical: angle between normal and
%      vertical direction. 
%
%      mean_angle_btw_unit_normals_and_vertical_driven_path
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_findEdge_determineAngleDeviation for a full
%       test suite.
%
% This function was written on 2024_08_13 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu
%
%%% Revision history
%

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==13 && isequal(varargin{end},-1))
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
        narginchk(10,13);
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (10<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end


% Does user want to specify another fig_num?
if (0==flag_max_speed) &&  (11<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        fig_num2 = temp;
    end
end

% Does user want to specify third fig_num?
if (0==flag_max_speed) &&  (12<=nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num3 = temp;
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
input_points = LiDAR_allPoints(:,1:3); 

% The indices of the mapped grids are extracted and concatenated 
original_mapped_gridIndices_cell = gridIndices_cell_array(original_qualified_grids); 

% Total number of mapped grids
total_mapped_grids = length(original_mapped_gridIndices_cell); 

% Unit normal vectors of the plane fits of each mapped grid
unit_normal_vectors = zeros(total_mapped_grids,3); 


for ith_mapped_grid = 1:total_mapped_grids
    [~, ~, ~, unit_normal_vectors(ith_mapped_grid,:), ~, ~] =...
    fcn_findEdge_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),-1);
    % mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
    % z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
end

% STEP 2
% Comparing normal vector with vertical direction
unit_vector_vertical_direction = [0 0 1];

% The dot product is computed to find the angle between the vectors
dot_product = sum(unit_normal_vectors.*unit_vector_vertical_direction,2);

% The angle between unit vertical and the unit_normal_vector is computed to
% determine how close the normal vector is to vertical direction.
angle_btw_unit_normals_and_vertical = acos(dot_product);

% Find driven path indices in current mapped grids
driven_path_grid_indices_in_current_mapped_grids = ismember(current_qualified_grids,current_grid_numbers_of_driven_path);

% Standard deviation in Z of driven path grids
angle_btw_unit_normals_and_vertical_driven_path = angle_btw_unit_normals_and_vertical(driven_path_grid_indices_in_current_mapped_grids);

% Standard deviation in Z of other mapped grids
angle_btw_unit_normals_and_vertical_other_mapped_grids = angle_btw_unit_normals_and_vertical(~driven_path_grid_indices_in_current_mapped_grids); 

% Find mean std in z of driven path
mean_angle_btw_unit_normals_and_vertical_driven_path = mean(angle_btw_unit_normals_and_vertical_driven_path); 

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

hold on
axis on
grid on 
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers mapped grids in ENU')

% % plot(gridCenters_required_point_density(:,1), gridCenters_required_point_density(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);
% plot(gridCenters_drivable_grids(:,1), gridCenters_drivable_grids(:,2), '.','MarkerSize',50,'Color',[0 1 0]);
% plot(gridCenters_non_drivable_grids(:,1), gridCenters_non_drivable_grids(:,2), '.','MarkerSize',50,'Color',[1 0 0]);

plot(gridCenters_qualified_grids(:,1), gridCenters_qualified_grids(:,2), '.','MarkerSize',45,'Color',[0.2 0.2 0.2]);

for ith_text = 1:length(current_qualified_grids(:,1))
    current_text = sprintf('%.0d',current_qualified_grids(ith_text));
     % Place the text on the grid center
    text(gridCenters_qualified_grids(ith_text,1), gridCenters_qualified_grids(ith_text,2),current_text,'Color',[1 1 0],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',20,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

% Plot grid lines and standard deviation

figure(fig_num2)

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Angle deviation of mapped grids')

% Pre-allocation: To find the total number of points in each grid
total_points_in_mapped_grids = zeros(length(original_qualified_grids),1);

% Pre-allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied
gridlines_mapped_grids = zeros(11*length(original_qualified_grids),2); % length(gridlines) = 11

for ith_domain = 1:length(original_qualified_grids)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);

    current_color = [0.5 0.5 0.5];
    % Plot current AABB
    current_AABB = grid_AABBs(original_qualified_grids(ith_domain),1:4);

    % Nudge the current AABB inward
    current_AABB = current_AABB + grid_size/100*[1 -1 1 -1];

    % Calculate the gridlines
    gridlines = [...
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,1) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,2) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,4); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,3); ...
        current_AABB(1,2) current_AABB(1,3); ...
        nan nan;
        current_AABB(1,1) current_AABB(1,4); ...
        current_AABB(1,2) current_AABB(1,4); ...
        ];

    % Get all points in this domain and plot them
    rows_in_domain = gridIndices==original_qualified_grids(ith_domain);
    points_in_domain = input_points(rows_in_domain,:);

    total_points_in_mapped_grids(ith_domain) = length(points_in_domain);

    length_gridlines = length(gridlines);
    % % Plot the mapped points green
    % p1 = plot(points_in_original_mapped_grids(:,1),points_in_original_mapped_grids(:,2),'.','MarkerSize',20,'Color',[0.4660 0.6740 0.1880]);

    % Plot the result
    % p1 = plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    % hold on
    % plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)
    gridlines_mapped_grids(1+(ith_domain-1)*length_gridlines:ith_domain*length_gridlines,:) = gridlines;

end

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',30,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

% [0.4660 0.6740 0.1880] - Green (Not so bright)
angle_btw_unit_normals_and_vertical_round = round((angle_btw_unit_normals_and_vertical*180/pi),3);
% mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1));
for ith_text = 1:length(gridCenters_qualified_grids(:,1))
    current_text = sprintf('%.3f',angle_btw_unit_normals_and_vertical_round(ith_text));
    % Place the text on the grid center
    text(gridCenters_qualified_grids(ith_text,1), gridCenters_qualified_grids(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 10, 'FontWeight','bold','FontSmoothing','on');
end



% Plot grid lines and standard deviation

figure(fig_num3)

hold on
grid on
xlabel('Mapped grid centers')
ylabel('Standard deviation in Z')
title('Mapped grid centers vs standard deviation in Z ')

plot(current_qualified_grids, angle_btw_unit_normals_and_vertical*180/pi,'.','MarkerSize',30,'Color',[0.2 0.2 0.2])
plot(current_grid_numbers_of_driven_path, angle_btw_unit_normals_and_vertical_driven_path*180/pi,'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',1.5)
plot(current_qualified_grids(~driven_path_grid_indices_in_current_mapped_grids), angle_btw_unit_normals_and_vertical_other_mapped_grids*180/pi,'.','MarkerSize',10,'Color',[1 0 0])





% Find mean std in z of not driven path
% max_angle_btw_unit_normals_and_vertical_not_driven_path = max(angle_btw_unit_normals_and_vertical_other_mapped_grids); 

% Theta threshold
% RATIO: Find a ratio
% theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 0.04;
% 
% disp(theta_threshold*(180/pi))
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
