function  [transverse_span_threshold,transverse_span_each_grid]= fcn_findEdge_determineTransverseSpanThreshold...
    (grids_greater_than_zero_points, grid_AABBs, grid_size, gridIndices, input_points, LIDAR_scanLines, gridCenters_greater_than_zero_point_density,varargin)
%% fcn_findEdge_classifyGridsBasedOnTransverseSpan  Finding the orthogonal distances of points in the remaining rings by projecting orthogonally from one ring 
%
% FORMAT:
%
%      [transverse_span_threshold,transverse_span_each_grid]= fcn_findEdge_determineTransverseSpanThreshold...
%      (grids_greater_than_zero_points, grid_AABBs, grid_size, gridIndices, input_points, LIDAR_scanLines, fig_num, fig_num2, fig_num3)
%
% INPUTS:
%      grids_greater_than_zero_points: girds that contains greater than zero
%      points, the girds not empty.
%
%      grid_AABBs: 
%
%      grid_size: the size of grid
%
%      gridIndices: the indices of grids 
%
%     input_points: selected LIDAR points 
%
%     LIDAR_scanLines: the scan lines of LIDAR 
%
% (OPTIONAL INPUTS):
%
%      (fig_num): figure number, can be set to -1 for fast mode
% 
%      (fig_num2): figure number, can be set to -1 for fast mode
%
%      (fig_num3): figure number, can be set to -1 for fast mode
%
% OUTPUTS: 
%
%      transverse_span_threshold: Threshold of transverse span
%
%      transverse_span_each_grid: Threshold of transverse span in each
%      grid.
%
%
% DEPENDENCIES:
%       
%      (none)
%
% EXAMPLES:
%   See the script:
%   
%   script_test_fcn_findEdge_determineTransverseSpanThreshold
%
%   for a full test suite
%
% This function was written by Jiabao Zhao on 2024_08_14
%
% Revision history:
% 2024_07_28 - Aneesh Batchu
% -- wrote the code originally
% 2024_08_14 - Jiabao Zhao
% -- Functionalized this code

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

if flag_max_speed == 0
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(7,10);
    end
end 


% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (7<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

% Does user want to specify another fig_num?
if (0==flag_max_speed) &&  (8<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        fig_num2 = temp;
    end
end

% Does user want to specify third fig_num?
if (0==flag_max_speed) &&  (9<=nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num3 = temp;
    end
end

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
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Gridlines and points of the grids greater than zero point density')


% allocation: these are the coordinates of the corners of each grid
% line. Total no of lines required for each grid including NaNs is 11. Therefore, 11 is multiplied 
gridlines_grids_greater_than_zero = zeros(11*length(grids_greater_than_zero_points),2); % length(gridlines) = 11

concatenate_gridPoints_scanLines_rings = [];
orthogonal_dist_each_grid = []; 
transverse_span_each_grid = [];

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(grids_greater_than_zero_points)
    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.2 0.2 0.2];

    % Plot current AABB
    current_AABB = grid_AABBs(grids_greater_than_zero_points(ith_grid),1:4);

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
    rows_in_domain = gridIndices==grids_greater_than_zero_points(ith_grid);
    
    % XY coordinates of the input points that are in ith_grid
    points_in_domain = input_points(rows_in_domain,:);
 
    % Scan lines and rings
    scanLines_and_rings = LIDAR_scanLines(rows_in_domain,:);

    % Find number of LiDAR scan lines in each grid
    scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));
    
    % Combine points_in_domain and ScanLines and rings
    gridPoints_scanLines_rings_to_add = [ith_grid*(ones(length(points_in_domain(:,1)),1)),points_in_domain,scanLines_and_rings];

    % Sort gridPoints_scanLines_rings_to_add based on scan line
    [~,sorted_gridPoints_scanLines_rings] = sort(gridPoints_scanLines_rings_to_add(:,4));

    % Sorted gridPoints_scanLines_rings_to_add matrix
    gridPoints_scanLines_rings_to_add = gridPoints_scanLines_rings_to_add(sorted_gridPoints_scanLines_rings,:);

    % Count occurrences of each unique number in scan lines
    [uniqueNumbers, ~, uniqueNum_idx] = unique(gridPoints_scanLines_rings_to_add(:,4));
    counts = histc(uniqueNum_idx, 1:numel(uniqueNumbers));

    % Create the new array with the counts
    count_array = counts;

    % Index of the scan line with more than one occurence
    index_of_scanLines = find(count_array>1, 1);
    
    if ~isempty(index_of_scanLines)
        % Indices first scan line of the matrix as a seperate matrix
        indices_gridPoints_scanLines_first_scan = find(gridPoints_scanLines_rings_to_add(:,4) == gridPoints_scanLines_rings_to_add(index_of_scanLines,4));

        % Seperate the scan line with more than one occurence of the matrix as a seperate matrix
        gridPoints_scanLines_first_scan = gridPoints_scanLines_rings_to_add(indices_gridPoints_scanLines_first_scan,:);

        % Count occurrences of each unique number in rings
        [uniqueNumbers, ~, uniqueNum_idx] = unique(gridPoints_scanLines_first_scan(:,5));
        counts = histc(uniqueNum_idx, 1:numel(uniqueNumbers));

        % Create the new array with the counts
        count_array = counts;

        % Index of the scan line with more than one occurence
        index_of_rings = find(count_array>1, 1);

        % if length(gridPoints_scanLines_first_scan(:,1)) == 1
        %
        %     indices_gridPoints_scanLines_first_scan

        if length(gridPoints_scanLines_first_scan(:,1)) > 1 & (gridPoints_scanLines_first_scan(index_of_rings,5) == gridPoints_scanLines_first_scan(index_of_rings+1,5))
            change_in_vector = gridPoints_scanLines_first_scan(2,2:3) - gridPoints_scanLines_first_scan(1,2:3);
            unit_change_in_vector = fcn_INTERNAL_calcUnitVector(change_in_vector);
            orth_unit_change_in_vector = unit_change_in_vector*[0 1; -1 0];

            % The remaining number of grids
            remaining_grids = length(indices_gridPoints_scanLines_first_scan)+1:length(gridPoints_scanLines_rings_to_add(:,1));

            %
            vector_from_base_point_first_scan_to_points_in_otherScans_rings = gridPoints_scanLines_rings_to_add(remaining_grids,2:3) - ...
                gridPoints_scanLines_first_scan(1,2:3).*(ones(length(remaining_grids),2));

            % Unit orthogonal vector
            repeated_orth_unit_change_in_vector = orth_unit_change_in_vector.*(ones(length(remaining_grids),2));

            % Calculate the transverse distance
            transverse_dist_grid_points_other_scanLines = sum(vector_from_base_point_first_scan_to_points_in_otherScans_rings.*repeated_orth_unit_change_in_vector,2);

            % Positive transverse distances
            positive_transverse_dist_grid_points_other_scanLines = transverse_dist_grid_points_other_scanLines(transverse_dist_grid_points_other_scanLines>=0);

            % Negative transverse distances
            negative_transverse_dist_grid_points_other_scanLines = transverse_dist_grid_points_other_scanLines(transverse_dist_grid_points_other_scanLines<0);

            % maximum span distance
            if ~isempty(positive_transverse_dist_grid_points_other_scanLines) && ~isempty(negative_transverse_dist_grid_points_other_scanLines)

                maximum_span_distance = max(positive_transverse_dist_grid_points_other_scanLines) + max(abs(negative_transverse_dist_grid_points_other_scanLines));

            elseif ~isempty(positive_transverse_dist_grid_points_other_scanLines) && isempty(negative_transverse_dist_grid_points_other_scanLines)

                maximum_span_distance = max(positive_transverse_dist_grid_points_other_scanLines);

            elseif isempty(positive_transverse_dist_grid_points_other_scanLines) && ~isempty(negative_transverse_dist_grid_points_other_scanLines)

                maximum_span_distance = max(abs(negative_transverse_dist_grid_points_other_scanLines));

            elseif isempty(positive_transverse_dist_grid_points_other_scanLines) && isempty(negative_transverse_dist_grid_points_other_scanLines)

                maximum_span_distance = 0;
                    
            end
            % Conatenate maximum transverse span
            transverse_span_each_grid = [transverse_span_each_grid; maximum_span_distance]; %#ok<AGROW>

            % Mean of absolute values of transverse distances
            mean_dist = mean(abs(transverse_dist_grid_points_other_scanLines));

            % Concatenate the orthogonal distances
            orthogonal_dist_each_grid = [orthogonal_dist_each_grid; mean_dist]; %#ok<AGROW>
        else

            % Conacatenate orthogonal distance
            orthogonal_dist_each_grid = [orthogonal_dist_each_grid; 0]; %#ok<AGROW>

            % Conatenate maximum transverse span
            transverse_span_each_grid = [transverse_span_each_grid; 0]; %#ok<AGROW>

        end

    else

        % Conacatenate orthogonal distance
        orthogonal_dist_each_grid = [orthogonal_dist_each_grid; 0]; %#ok<AGROW>

        % Conatenate maximum transverse span
        transverse_span_each_grid = [transverse_span_each_grid; 0]; %#ok<AGROW>

    end
    % Concatenate the points in domain, scan lines and rings
    concatenate_gridPoints_scanLines_rings = [concatenate_gridPoints_scanLines_rings; gridPoints_scanLines_rings_to_add]; %#ok<AGROW>

    % Save the grid lines of all the grids greater than zero density in a
    % matrix
    length_gridlines = length(gridlines);
    gridlines_grids_greater_than_zero(1+(ith_grid-1)*length_gridlines:ith_grid*length_gridlines,:) = gridlines;
   
    % Plot the result
    % plot the grid lines of the grids greater than zero
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);
    % plot the points in the grid
    plot(points_in_domain(:,1),points_in_domain(:,2),'.','Color',[0.2,0.2,0.2],'MarkerSize',10)

end

% Write the grid number at the grid center for reference. 

for ith_text = 1:length(grids_greater_than_zero_points(:,1))
    % current_text = sprintf('%.0d',current_mapped_grids(ith_text));
    current_text = sprintf('%.0d',(ith_text));

    % Place the text on the grid center
    text(gridCenters_greater_than_zero_point_density(ith_text,1), gridCenters_greater_than_zero_point_density(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end


figure(fig_num2); 
hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(grids_greater_than_zero_points)

    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.2 0.2 0.2];

    % Plot current AABB
    current_AABB = grid_AABBs(grids_greater_than_zero_points(ith_grid),1:4);

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

    % Plot the result
    % plot the grid lines of the grids greater than zero
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);

end

% Write the grid number at the grid center for reference. 

for ith_text = 1:length(grids_greater_than_zero_points(:,1))
    current_text = sprintf('%.3f',orthogonal_dist_each_grid(ith_text));
    % current_text = sprintf('%.3f',orthogonal_dist_each_grid);

    % Place the text on the grid center
    text(gridCenters_greater_than_zero_point_density(ith_text,1), gridCenters_greater_than_zero_point_density(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end


figure(fig_num3); 
hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Transverse span distances')

% total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(grids_greater_than_zero_points)

    % Get current color
    % current_color = fcn_geometry_fillColorFromNumberOrName(ith_domain);
    current_color = [0.2 0.2 0.2];

    % Plot current AABB
    current_AABB = grid_AABBs(grids_greater_than_zero_points(ith_grid),1:4);

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

    % Plot the result
    % plot the grid lines of the grids greater than zero
    plot(gridlines(:,1),gridlines(:,2),'-','Color',current_color,'LineWidth',3);

end

% Write the grid number at the grid center for reference. 

for ith_text = 1:length(grids_greater_than_zero_points(:,1))
    current_text = sprintf('%.3f',transverse_span_each_grid(ith_text));
    % current_text = sprintf('%.3f',orthogonal_dist_each_grid);

    % Place the text on the grid center
    text(gridCenters_greater_than_zero_point_density(ith_text,1), gridCenters_greater_than_zero_point_density(ith_text,2),current_text,'Color',[0.5, 0, 0.5],'HorizontalAlignment','center','FontSize', 8, 'FontWeight','bold');
end

% Threshold of transverse span
transverse_span_threshold = 0.15;
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

%% fcn_INTERNAL_calcUnitVector
function unit_vectors = fcn_INTERNAL_calcUnitVector(input_vectors)
vector_length = sum(input_vectors.^2,2).^0.5;
unit_vectors = input_vectors./vector_length;
end