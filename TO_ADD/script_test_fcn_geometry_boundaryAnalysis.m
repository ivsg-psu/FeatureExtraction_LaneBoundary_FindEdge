%% script_test_fcn_geometry_boundaryAnalysis
%
% This is the main script that performs boundary analysis
%
% The current script follows to the latest schema. Each step (section)
% corresponds to a different block in the flowchart (schema). The parallel
% steps—tasks that are not dependent on each other and can be run
% independently—are numbered identically.

% Revision History
% 2024_08_04 - Aneesh Batchu
% -- wrote the code originally 


%% STEP 1: Load and study the data (Yet to be functionalized)

% This is done in "script_load_sample_LIDAR_data" to load vehicle pose and
% LIDAR data
% Run script_test_geometry_PointsAtRangeOfLiDARFromStation

% These scripts are in the main directory of Geom Class
% Just "run" the scripts. You don't need to run section by section



%% STEP 2: Run script_plot_sample_data.m to capture enough data

% Then, run script_test_geometry_findPointsInDomain. 

% concatenate_LiDAR_XYZ_points
% concatenate_scanLine_rings

%% STEP 2: Find the boundary points of the driven path

% Run script_test_geometry_boundaryPointsDrivenPath to find the boundary
% points

%% STEP 3: Seperate the data into grids

% Figure number
fig_num = 40; 
figure(fig_num);clf

LiDAR_allPoints = [concatenate_LiDAR_XYZ_points(in_domain,:), concatenate_scanLine_rings(in_domain,:)];

% remove NANs
LiDAR_allPoints = LiDAR_allPoints(~isnan(LiDAR_allPoints(:,1)),:);

% scanLine_rings without NaNs
LIDAR_scanLines = LiDAR_allPoints(:,4:5); 

% Input points for seperating the data into grids. The points are in 2D as
% the analysis is carried out in 2D
input_points = LiDAR_allPoints(:,1:2); 

% grid size of each grid in 3 dimensions. For example, grid_size = 1.25:
% [length, width, height] = [1.25 1.25 1.25]
grid_size = 1; %0.8;%1;%1.25; %1.26

% Find minimum and maximum values of x,y,z of LiDAR data
[Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_geometry_findMaxMinOfXYZ(input_points,-1);

% The grids are plotted only within the boundaries in #D
% grid_boundaries = [Min_x Max_x Min_y Max_y Min_z Max_z]; 

% The 2D grid boundaries required for the current analysis
grid_boundaries = [Min_x Max_x Min_y Max_y]; 


% Seperate the LiDAR data into grids
[gridIndices_cell_array, total_N_points_in_each_grid, gridCenters, grids_with_zero_points,...
    grids_greater_than_zero_points, gridCenters_zero_point_density,...
    gridCenters_greater_than_zero_point_density, gridIndices, grid_AABBs] = fcn_geometry_findGridsWithPoints(input_points,...
    grid_size,grid_boundaries,fig_num);

%% STEP 4: Find the driven path grids within the grids more than zero points

% -----------------------------NOTE------------------------------
% After finding the grids without anypoints, the grids are completely
% removed from the analysis. Only, grids with greater than zero points were
% analyzed from here. 
% -----------------------------NOTE------------------------------

% Plot all the grids greater than zero point density

% "inpolygon" is used to find the grids within the boundary points 
[in,on] = inpolygon(gridCenters_greater_than_zero_point_density(:,1),gridCenters_greater_than_zero_point_density(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));

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

fig_num = 51;
figure(fig_num); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers and boundary points')

plot(gridCenters_greater_than_zero_point_density(:,1), gridCenters_greater_than_zero_point_density(:,2), '.','MarkerSize',30,'Color',[0.8 0.8 0.8]);
plot(gridCenters_driven_path(:,1), gridCenters_driven_path(:,2), '.','MarkerSize',30,'Color',[0 1 0]);


% Plot the grids with 
fig_num = 8765; 
figure(fig_num);clf

% plot computed boundary points
marker_size = 10;
RGB_triplet = [0 0 0]; 
legend_option = 1;
legend_name = 'Grids greater than zero points';
legend_position = [];
marker_type = [];
plot_gridCenters_greater_than_zero_point_density = [gridCenters_greater_than_zero_point_density(:,1:2), zeros(length(gridCenters_greater_than_zero_point_density(:,1)),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_greater_than_zero_point_density,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot driven path
marker_size = 25;
RGB_triplet = [0 1 0]; 
legend_option = 1;
legend_name = 'Driven path grids';
legend_position = [];
marker_type = [];
plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 5: Grid conditions - Point density (Do not need to run after determining point density)

% Figure number of histogram
fig_num = 52; 
figure(fig_num); clf; 
% Add labels and title 
hold on
grid on
xlabel('Points per grid');
ylabel('Frequency');
title('Statistic 1: Determining suitable point density');


% edges = (floor(min(total_points_in_each_grid_with_points_greater_than_zero)/10)*10):10:(ceil(max(total_points_in_each_grid_with_points_greater_than_zero)/10)*10); % Define the bin edges

% Create the histogram
% actual_driving_surface_grids_hist = histogram(total_points_in_each_grid_in_the_driven_path,edges,'Visible','on'); 

% total_grids_hist = histogram(total_points_in_each_grid_with_points_greater_than_zero,edges,'Visible','on'); 
total_grids_greater_than_zero_hist = histogram(total_points_in_each_grid_with_points_greater_than_zero,20,'Visible','on'); 
actual_driven_path_grids_hist = histogram(total_points_in_each_grid_in_the_driven_path,10,'Visible','on'); 

% Extract the counts for both histograms
counts1 = total_grids_greater_than_zero_hist.Values;
[~,index_max_counts2] = max(counts1);
counts2 = actual_driven_path_grids_hist.Values;

binEdges = total_grids_greater_than_zero_hist.BinEdges;

% Calculate the overlapping counts
% % overlapCounts = min(counts2, counts1);

% Find a ratio
 % point_density = sum(binEdges(index_max_counts2:(index_max_counts2+1)))/2; 

% point_density = floor(mean(total_points_in_each_grid_in_the_driven_path) - 7.5*(std(total_points_in_each_grid_in_the_driven_path)));
x_location = floor(mean(total_points_in_each_grid_in_the_driven_path) - 1.5*(std(total_points_in_each_grid_in_the_driven_path)));

% point_density = floor(mean(total_points_in_each_grid_in_the_driven_path));


% Minimum number of points required 
point_density = floor(20*((grid_size^2)/(0.3^2))); 
disp('Chosen point density')
disp(point_density)
% mean(total_points_in_each_grid_in_the_driven_path)/mean(total_points_in_each_grid_with_points_greater_than_zero);

plot(point_density,0, 'k.', 'MarkerSize',20)
current_text = sprintf('point density = %.2d',point_density);
text(x_location, 200,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');

%% STEP 5: Grid conditions - Determining number of scan lines in each grid greater than zero points Point density 

% This was also done in STEP 5, however, it was done using a for loop. Need to
% do it without using a for loop.

total_scan_lines_in_each_grid_with_more_than_zero_points = [];
for ith_grid = 1:length(grids_greater_than_zero_points)
     %Get all points in this domain and plot them
    rows_in_domain = gridIndices==grids_greater_than_zero_points(ith_grid);

     %Find number of LiDAR scan lines in each grid
    scan_lines_ith_grid = length(unique(LIDAR_scanLines(rows_in_domain,1)));

    % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORK     S WITHOUT A FOR LOOP
    total_scan_lines_in_each_grid_with_more_than_zero_points = [total_scan_lines_in_each_grid_with_more_than_zero_points; scan_lines_ith_grid]; %#ok<AGROW> % SHOULD NOT BE IN FOR LOOP, NEED TO WRITE A CODE THAT WORKS WITHOUT A FOR LOOP
end

%% STEP 5: Grid conditions - Finding the orthogonal distances of points in the remaining rings by projecting orthogonally from one ring 

% Figure number
fig_num = 50;
figure(fig_num); clf;

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
            unit_change_in_vector = fcn_geometry_calcUnitVector(change_in_vector);
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

fig_num = 5837; 
figure(fig_num); clf;
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

fig_num = 58309; 
figure(fig_num); clf;
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

%% STEP 6: Grids with required point density and low point density

% These grids have low point density but not zero point density 
original_grids_with_low_point_density = grids_greater_than_zero_points((total_N_points_in_each_grid(grids_greater_than_zero_points,1) > 0) & (total_N_points_in_each_grid(grids_greater_than_zero_points,1) < point_density)); 

% These are grids that contain points more than or equal to point density
original_grids_with_required_point_density = grids_greater_than_zero_points(total_N_points_in_each_grid(grids_greater_than_zero_points,1) >= point_density);
grid_indices_with_required_point_density = (total_N_points_in_each_grid(grids_greater_than_zero_points,1) >= point_density);

% Current grid numbers of the grids with low point density
% These are the grid numbers of the grids with respect to grids_greater_than_zero_points
% The numbering starts with "1"

% These are the grid numbers of low point density
current_grids_with_low_point_density = find((total_N_points_in_each_grid(grids_greater_than_zero_points,1) > 0) & (total_N_points_in_each_grid(grids_greater_than_zero_points,1) < point_density)); 

% These are the current grid numbers of required point density
current_grids_with_required_point_density = find(total_N_points_in_each_grid(grids_greater_than_zero_points,1) >= point_density); 

% Grid Centers of the grids with low point density (Unmapped grid centers)
gridCenters_low_point_density = gridCenters(original_grids_with_low_point_density,1:2);  

% Grid Centers of the grids with required point density 
gridCenters_required_point_density = gridCenters(original_grids_with_required_point_density,1:2);  

% Plot the grids with low point density and required density 
fig_num = 70; 
figure(fig_num);clf

marker_size = 25;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Grids with low point density';
legend_position = [];
marker_type = [];
plot_gridCenters_low_point_density = [gridCenters_low_point_density, zeros(length(gridCenters_low_point_density(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_low_point_density,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot grid centers
marker_size = 25;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Grids with required density';
legend_position = [];
marker_type = [];
plot_gridCenters_required_point_density = [gridCenters_required_point_density, zeros(length(gridCenters_required_point_density(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_required_point_density,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 6: Grids with more than one scan line and grids with one scan line

% These grids contain more than one scan line
original_grids_with_more_than_one_scan_line = grids_greater_than_zero_points(total_scan_lines_in_each_grid_with_more_than_zero_points>1);
grid_indices_with_more_than_one_scan_line = (total_scan_lines_in_each_grid_with_more_than_zero_points>1);

% These grids contain only one scan line
original_grids_with_one_scan_line = grids_greater_than_zero_points(~grid_indices_with_more_than_one_scan_line);

% Current grid numbers of the grids with more than one scan line
% These are the grid numbers of the grids with respect to grids_greater_than_zero_points
% The numbering starts with "1"

% Current grid numbers of the grids with more than one scan line
current_grids_with_more_than_one_scan_line = find(grid_indices_with_more_than_one_scan_line);

% Current grid numbers of the grids with one scan line
current_grids_with_one_scan_line = find(~grid_indices_with_more_than_one_scan_line);

% Grid centers of the grids with more than zero points and more than one scan line
gridCenters_with_more_than_one_scan_line = gridCenters(original_grids_with_more_than_one_scan_line,1:2);

% Grid centers of the grids with one scan line
gridCenters_with_one_scan_line = gridCenters(original_grids_with_one_scan_line,1:2);

% Plot grids with one scan line and more than one scan line
fig_num = 71; 
figure(fig_num);clf

marker_size = 25;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Grids with one scan line';
legend_position = [];
marker_type = [];
plot_gridCenters_with_one_scan_line = [gridCenters_with_one_scan_line, zeros(length(gridCenters_with_one_scan_line(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_with_one_scan_line,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot grid centers
marker_size = 25;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Grids with more than one scan line';
legend_position = [];
marker_type = [];
plot_gridCenters_with_more_than_one_scan_line = [gridCenters_with_more_than_one_scan_line, zeros(length(gridCenters_with_more_than_one_scan_line(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_with_more_than_one_scan_line,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 6: Grids greater than minimum transverse span threshold and lesser than minimum transverse span threshold 

% Grid indices of the original grids greater than transverse span threshold
grid_indices_with_more_than_transverse_span_threshold = (transverse_span_each_grid>transverse_span_threshold);

% The transverse span of these original grids are more than transverse span
% threshold
original_grids_with_more_than_transverse_span_threshold = grids_greater_than_zero_points(grid_indices_with_more_than_transverse_span_threshold);

% The transverse span of these grids are less than/equal to transverse span
% threshold
original_grids_with_less_than_transverse_span_threshold = grids_greater_than_zero_points(~grid_indices_with_more_than_transverse_span_threshold);

% Current grid numbers of the grids with more than one scan line
% These are the grid numbers of the grids with respect to grids_greater_than_zero_points
% The numbering starts with "1"

% Current grid numbers of the grids with more than transverse span
% threshold
current_grids_with_more_than_transverse_span_threshold = find(grid_indices_with_more_than_transverse_span_threshold);

% Current grid numbers of the grids with less than/equal transverse span
% threshold
current_grids_with_less_than_transverse_span_threshold = find(~grid_indices_with_more_than_transverse_span_threshold);

% Grid centers of the grids with more than transerse span threshold
gridCenters_with_more_than_transverse_span_threshold = gridCenters(original_grids_with_more_than_transverse_span_threshold,1:2);

% Grid centers of the grids with less than/equal transverse span threshold
gridCenters_with_less_than_transverse_span_threshold = gridCenters(original_grids_with_less_than_transverse_span_threshold,1:2);

% Plot grids greater than minimum transverse span threshold and lesser than minimum transverse span threshold 

fig_num = 72; 
figure(fig_num);clf

marker_size = 25;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Grids lesser than minimum transverse span threshold';
legend_position = [];
marker_type = [];
plot_gridCenters_with_less_than_transverse_span_threshold = [gridCenters_with_less_than_transverse_span_threshold, zeros(length(gridCenters_with_less_than_transverse_span_threshold(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_with_less_than_transverse_span_threshold,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot grid centers
marker_size = 25;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Grids greater than minimum transverse span threshold';
legend_position = [];
marker_type = [];
plot_gridCenters_with_more_than_transverse_span_threshold = [gridCenters_with_more_than_transverse_span_threshold, zeros(length(gridCenters_with_more_than_transverse_span_threshold(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_with_more_than_transverse_span_threshold,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 6: Qualified and unqualified grids: grids that pass all three conditions above are qualified


% Grid indices of original grids with required point density, more than one
% scan line, and more than transverse span threshold
grid_indices_qualified_grids = (grid_indices_with_required_point_density == 1) & (grid_indices_with_more_than_one_scan_line == 1) & (grid_indices_with_more_than_transverse_span_threshold == 1);

% --------------------------Qualified grid numbers-----------------------------------
% The grids with more than one scan line, more/equal to point density and
% greater than minimum transverse span threshold
original_qualified_grids = grids_greater_than_zero_points(grid_indices_qualified_grids); 

% --------------------------Unqualified grid numbers-----------------------------------
% Find grids with low point density, grids with only one LiDAR scan and
% lesser than minimum transverse span threshold
% The grids that are not qualified grids are clasified as unqualified grids
original_unqualified_grids = grids_greater_than_zero_points(~grid_indices_qualified_grids); 

% Current qualified grid numbers 
current_qualified_grids = find(grid_indices_qualified_grids); 

% Current unqualified grid numbers 
current_unqualified_grids = find(~grid_indices_qualified_grids); 

% Grid centers of qualified grids
gridCenters_qualified_grids = gridCenters(original_qualified_grids,1:2); 

% Grid centers of unqualified grids
gridCenters_unqualified_grids = gridCenters(original_unqualified_grids,1:2); 


% Plot qualified and unqualified grids

fig_num = 27; 
figure(fig_num);clf

marker_size = 25;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Unqualified grids';
legend_position = [];
marker_type = [];
plot_gridCenters_unqualified_grids = [gridCenters_unqualified_grids, zeros(length(gridCenters_unqualified_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_unqualified_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot grid centers
marker_size = 25;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Qualified grids';
legend_position = [];
marker_type = [];
plot_gridCenters_qualified_grids = [gridCenters_qualified_grids, zeros(length(gridCenters_qualified_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_qualified_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 6: Plot circles corresponding to each fail condition
% less than required point density - Red (small circle)
% less than one scan line - Green (medium circle)
% less than transverse span threshold - Blue (Large circle)

% Plot qualified and unqualified grids

fig_num = 987; 
figure(fig_num);clf

marker_size = 10;
RGB_triplet = [0.8, 0.8, 0.8]; 
legend_option = 1;
legend_name = 'Unqualified grids';
legend_position = [];
marker_type = [];
plot_gridCenters_unqualified_grids = [gridCenters_unqualified_grids, zeros(length(gridCenters_unqualified_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_unqualified_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot grid centers
marker_size = 10;
RGB_triplet = [0.2, 0.2, 0.2]; 
legend_option = 1;
legend_name = 'Qualified grids';
legend_position = [];
marker_type = [];
plot_gridCenters_qualified_grids = [gridCenters_qualified_grids, zeros(length(gridCenters_qualified_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_qualified_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% Plot the fail conditions
marker_size = 5;
RGB_triplet = [1, 0, 0]; 
legend_option = 1;
legend_name = 'Grids with low point density';
legend_position = [];
marker_type = 'o';
plot_gridCenters_low_point_density = [gridCenters_low_point_density, zeros(length(gridCenters_low_point_density(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_low_point_density,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


marker_size = 7.5;
RGB_triplet = [0, 1, 0]; 
legend_option = 1;
legend_name = 'Grids with one scan line';
legend_position = [];
marker_type = 'o';
plot_gridCenters_with_one_scan_line = [gridCenters_with_one_scan_line, zeros(length(gridCenters_with_one_scan_line(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_with_one_scan_line,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


marker_size = 10;
RGB_triplet = [0, 0, 1]; 
legend_option = 1;
legend_name = 'transverse span threshold fail';
legend_position = [];
marker_type = 'o';
plot_gridCenters_with_less_than_transverse_span_threshold = [gridCenters_with_less_than_transverse_span_threshold, zeros(length(gridCenters_with_less_than_transverse_span_threshold(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_with_less_than_transverse_span_threshold,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 7: Recalculate the driven path grids among qualified grids

% "inpolygon" is used to find the grids within the boundary points 
[in_qg,on_qg] = inpolygon(gridCenters_qualified_grids(:,1),gridCenters_qualified_grids(:,2),boundary_points_driven_path(:,1),boundary_points_driven_path(:,2));

% Original grid numbers of driven path
original_grid_numbers_of_driven_path = original_qualified_grids(in_qg); 

% Current grid numbers in driven path 
current_grid_numbers_of_driven_path = current_qualified_grids(in_qg); %find(in); 

% % Total points in each grid in the driven path
% total_points_in_each_grid_in_the_driven_path = total_N_points_in_each_grid(original_grid_numbers_of_driven_path); 
% 
% % Total points in each grid with points greater than zero
% total_points_in_each_grid_with_points_greater_than_zero = total_N_points_in_each_grid(current_qualified_grids); 

% Grid centers of the driven path
gridCenters_driven_path = [gridCenters_qualified_grids(in_qg,1),gridCenters_qualified_grids(in_qg,2)];


fig_num = 517;
figure(fig_num); clf;

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Grid centers and boundary points')

plot(gridCenters_qualified_grids(:,1), gridCenters_qualified_grids(:,2), '.','MarkerSize',40,'Color',[0.2 0.2 0.2]);
% plot(boundary_points_driven_path(:,1), boundary_points_driven_path(:,2), '.', 'MarkerSize',30, 'Color',[0 1 0]); 

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',2) % points strictly inside

% for ith_text = 1:length(current_qualified_grids(:,1))
%     current_text = sprintf('%.0d',ith_text);
%     % Place the text on the grid center
%     text(gridCenters_qualified_grids(ith_text,1), gridCenters_qualified_grids(ith_text,2),current_text,'Color',[1 1 1],'HorizontalAlignment','center','FontSize', 6, 'FontWeight','bold');
% end


% Plot the grids with 
fig_num = 804; 
figure(fig_num);clf

% plot computed boundary points
marker_size = 10;
RGB_triplet = [0 0 0]; 
legend_option = 1;
legend_name = 'Qualified grids';
legend_position = [];
marker_type = [];
plot_gridCenters_qualified_grids = [gridCenters_qualified_grids(:,1:2), zeros(length(gridCenters_qualified_grids(:,1)),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_qualified_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot driven path
marker_size = 25;
RGB_triplet = [0 1 0]; 
legend_option = 1;
legend_name = 'Driven path grids';
legend_position = [];
marker_type = [];
plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 8: Qualified grid conditions - Standard deviation in Z (Do not need to run after determining a std_threshold)

input_points =  LiDAR_allPoints(:,1:3);

% The indices of the mapped grids are extracted and concatenated 
original_mapped_gridIndices_cell = gridIndices_cell_array(original_qualified_grids); 

% Total number of mapped grids
total_mapped_grids = length(original_mapped_gridIndices_cell); 

% Standard deviations in orthogonal distances of points in the grid to
% plane
% standard_deviation_in_plane_orthogonals = zeros(total_mapped_grids,1); 
standard_deviation_in_z = zeros(total_mapped_grids,1); 

for ith_mapped_grid = 1:total_mapped_grids
    % points = input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:);
    % points = points(~isnan(points(:,1)),:);
    [~, standard_deviation_in_z(ith_mapped_grid,:), ~, ~, ~, ~] =...
    fcn_geometry_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),-1);
    % mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
    % z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
end

fig_num = 80; 
figure(fig_num);clf;

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

% plot the grids in the driven path
plot(gridCenters_driven_path(:,1),gridCenters_driven_path(:,2),'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',2) % points strictly inside


for ith_text = 1:length(current_qualified_grids(:,1))
    current_text = sprintf('%.0d',(ith_text));
     % Place the text on the grid center
    text(gridCenters_qualified_grids(ith_text,1), gridCenters_qualified_grids(ith_text,2),current_text,'Color',[1 1 1],'HorizontalAlignment','center','FontSize', 6, 'FontWeight','bold');
end



% Plot grid lines and standard deviation
fig_num = 81;
figure(fig_num)

hold on
grid on
xlabel('X[m]')
ylabel('Y[m]')
title('Standard deviation of Z of mapped grids')

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
standard_deviation_in_z_round = round(standard_deviation_in_z,3);
% mapped_grid_numbers = 1:length(gridCenters_required_point_density(:,1));
for ith_text = 1:length(gridCenters_qualified_grids(:,1))
    current_text = sprintf('%.3f',standard_deviation_in_z_round(ith_text));
    % Place the text on the grid center
    text(gridCenters_qualified_grids(ith_text,1), gridCenters_qualified_grids(ith_text,2),current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 10, 'FontWeight','bold','FontSmoothing','on');
end


% Find driven path indices in current mapped grids
driven_path_grid_indices_in_current_mapped_grids = ismember(current_qualified_grids,current_grid_numbers_of_driven_path);

% Standard deviation in Z of driven path grids
std_in_z_driven_path = standard_deviation_in_z(driven_path_grid_indices_in_current_mapped_grids);

% Standard deviation in Z of other mapped grids
std_in_z_other_mapped_grids = standard_deviation_in_z(~driven_path_grid_indices_in_current_mapped_grids); 


% Plot grid lines and standard deviation
fig_num = 82;
figure(fig_num);clf;

hold on
grid on
xlabel('Qualified grid centers')
ylabel('Standard deviation in Z')
title('Qualified grid centers vs standard deviation in Z ')

plot(current_qualified_grids, standard_deviation_in_z,'.','MarkerSize',30,'Color',[0.2 0.2 0.2])
plot(current_grid_numbers_of_driven_path, std_in_z_driven_path,'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',1.5)
plot(current_qualified_grids(~driven_path_grid_indices_in_current_mapped_grids), std_in_z_other_mapped_grids,'.','MarkerSize',10,'Color',[1 0 0])


% Find mean std in z of driven path
mean_std_in_z_driven_path = mean(std_in_z_driven_path); 

% Find mean std in z of not driven path
mean_std_in_z_not_driven_path = mean(std_in_z_other_mapped_grids(~isnan(std_in_z_other_mapped_grids))); 

% Find max std in z of not driven path
max_std_in_z_not_driven_path = max(std_in_z_other_mapped_grids); 

% Std Threshold
% Instead of choosing 6, try to find a ratio
% ratio: mean_std_in_z_driven_path/mean_std_of_all_grids
% std_threshold = mean_std_in_z_driven_path*6;
% 
% disp(std_threshold)

%% STEP 9: Histogram of standard deviation - (Do not need to run after determining a std_threshold)

figure(123);clf;
hold on
grid on
xlabel('Standard deviation in z error after plane fit')
ylabel('Frequency')
title('Statistic 3: Determining suitable standard deviation in z')
% histogram(standard_deviation_in_z)
histogram(standard_deviation_in_z,100)
histogram(std_in_z_driven_path,5)

% std_threshold = mean_std_in_z_driven_path + 6*std(std_in_z_driven_path); 
x_coord = mean_std_in_z_driven_path + 100*std(std_in_z_driven_path); 
std_threshold = 0.1;%0.08; 
% plot(std_threshold,max(std_in_z_driven_path), 'k.', 'MarkerSize',20)
disp('mean of std_threshold of driven path')
disp(mean_std_in_z_driven_path)
disp('chosen std_threshold')
disp(std_threshold)
% std_threshold = 0.05; 
plot(std_threshold,0, 'k.', 'MarkerSize',18)
current_text = sprintf('std threshold = %.4f',std_threshold);
text(x_coord, 80,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');

%% STEP 8: Qualified grid conditions - angle deviation (Do not need to run after determining a theta_threshold)

input_points = LiDAR_allPoints(:,1:3); 

% The indices of the mapped grids are extracted and concatenated 
original_mapped_gridIndices_cell = gridIndices_cell_array(original_qualified_grids); 

% Total number of mapped grids
total_mapped_grids = length(original_mapped_gridIndices_cell); 

% Unit normal vectors of the plane fits of each mapped grid
unit_normal_vectors = zeros(total_mapped_grids,3); 


for ith_mapped_grid = 1:total_mapped_grids
    [~, ~, ~, unit_normal_vectors(ith_mapped_grid,:), ~, ~] =...
    fcn_geometry_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),-1);
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

fig_num_ENU_statistic_four = 83; 
figure(fig_num_ENU_statistic_four);clf;

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
fig_num = 84;
figure(fig_num)

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

% Find driven path indices in current mapped grids
driven_path_grid_indices_in_current_mapped_grids = ismember(current_qualified_grids,current_grid_numbers_of_driven_path);

% Standard deviation in Z of driven path grids
angle_btw_unit_normals_and_vertical_driven_path = angle_btw_unit_normals_and_vertical(driven_path_grid_indices_in_current_mapped_grids);

% Standard deviation in Z of other mapped grids
angle_btw_unit_normals_and_vertical_other_mapped_grids = angle_btw_unit_normals_and_vertical(~driven_path_grid_indices_in_current_mapped_grids); 


% Plot grid lines and standard deviation
fig_num = 85;
figure(fig_num)

hold on
grid on
xlabel('Mapped grid centers')
ylabel('Standard deviation in Z')
title('Mapped grid centers vs standard deviation in Z ')

plot(current_qualified_grids, angle_btw_unit_normals_and_vertical*180/pi,'.','MarkerSize',30,'Color',[0.2 0.2 0.2])
plot(current_grid_numbers_of_driven_path, angle_btw_unit_normals_and_vertical_driven_path*180/pi,'o','MarkerSize',10,'Color',[0 1 0], 'LineWidth',1.5)
plot(current_qualified_grids(~driven_path_grid_indices_in_current_mapped_grids), angle_btw_unit_normals_and_vertical_other_mapped_grids*180/pi,'.','MarkerSize',10,'Color',[1 0 0])


% Find mean std in z of driven path
mean_angle_btw_unit_normals_and_vertical_driven_path = mean(angle_btw_unit_normals_and_vertical_driven_path); 

% Find mean std in z of not driven path
max_angle_btw_unit_normals_and_vertical_not_driven_path = max(angle_btw_unit_normals_and_vertical_other_mapped_grids); 

% Theta threshold
% RATIO: Find a ratio
% theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 0.04;
% 
% disp(theta_threshold*(180/pi))

%% STEP 9: Histogram of angle deviation - (Do not need to run after determining a theta_threshold)

figure(1223);clf;
hold on
grid on
xlabel('Angle deviation in z error after plane fit')
ylabel('Frequency')
title('Statistic 3: Determining suitable angle deviation')
% histogram(standard_deviation_in_z)
histogram(angle_btw_unit_normals_and_vertical,100)
histogram(angle_btw_unit_normals_and_vertical_driven_path,5)

x_coord = mean_angle_btw_unit_normals_and_vertical_driven_path + 20*std(angle_btw_unit_normals_and_vertical_driven_path); 
% theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 3.3*std(angle_btw_unit_normals_and_vertical_driven_path); 
% theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 3*std(angle_btw_unit_normals_and_vertical_driven_path);

% plot(theta_threshold,max(angle_btw_unit_normals_and_vertical), 'k.', 'MarkerSize',20)
theta_threshold = 0.1745; 
% theta_threshold = 0.1504; 
% theta_threshold = 0.3; 

disp('mean of angle_deviation_driven_path')
disp(mean_angle_btw_unit_normals_and_vertical_driven_path*(180/pi))
disp('theta threshold')
disp(theta_threshold*(180/pi))

plot(theta_threshold,0, 'k.', 'MarkerSize',18)
current_text = sprintf('theta threshold = %.1f',theta_threshold*(180/pi));
text(x_coord, 30,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');

%% STEP 10: Voting - Drivable, Non-drivable and Uncertain

% NOTE: The code refers non-drivable grids as failed grids and non-drivable
% as (failed grids and uncertain grids) - NEED to change

input_points = LiDAR_allPoints(:,1:3); 

% Give the thresholds if already found
% std_threshold = 0.05; 
% theta_threshold = 7*pi/180;
% theta_threshold = 30*pi/180;
% gridCenters

fig_num = 6000; 
figure(fig_num);clf

% theta_threshold = [];
% std_threshold = [];


% Classify mapped grids into drivable and drivable
[standard_deviation_in_z, angle_btw_unit_normals_and_vertical, ...
    original_drivable_grids, original_non_drivable_grids, current_drivable_grid_numbers_in_mapped_grids, current_non_drivable_grid_numbers_in_mapped_grids, ...
    current_failed_grid_numbers_in_mapped_grids, current_uncertain_grid_numbers_in_mapped_grids, gridCenters_failed_grids, gridCenters_uncertain_grids, ...
    gridCenters_drivable_grids,gridCenters_non_drivable_grids, concatenate_gridCenters_drivable_non_drivable_grids] = ...
    fcn_geometry_classifyGridsAsDrivable(gridIndices_cell_array, original_qualified_grids, input_points, std_threshold, theta_threshold, gridCenters, fig_num);

%% STEP 10: Plot the drivable, non-drivable and uncertain grid centers 

% (The code refers non-drivable as failed) Non-drivable: Everything (uncertain and failed) other than drivable
 
% Plot qualified and unqualified grids

fig_num = 987; 
figure(fig_num);clf

marker_size = 25;
RGB_triplet = [0 1 0]; 
legend_option = 1;
legend_name = 'Drivable grids';
legend_position = [];
marker_type = [];
plot_gridCenters_drivable_grids = [gridCenters_drivable_grids(:,1:2), zeros(length(gridCenters_drivable_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_drivable_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot grid centers
marker_size = 25;
RGB_triplet = [1 0 0]; 
legend_option = 1;
legend_name = 'Non-drivable grids';
legend_position = [];
marker_type = [];
plot_gridCenters_failed_grids = [gridCenters_failed_grids(:,1:2), zeros(length(gridCenters_failed_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_failed_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot grid centers
marker_size = 25;
RGB_triplet = [0 0 1];%[0.9290 0.6940 0.1250]; 
legend_option = 1;
legend_name = 'Uncertain grids';
legend_position = [];
marker_type = [];
plot_gridCenters_uncertain_grids = [gridCenters_uncertain_grids(:,1:2), zeros(length(gridCenters_uncertain_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_uncertain_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


if ~isempty(std_threshold)
    % Find the indices of std threshold failed indices
    current_std_threshold_failed_indices = current_uncertain_grid_numbers_in_mapped_grids(standard_deviation_in_z(current_uncertain_grid_numbers_in_mapped_grids)>std_threshold);

    % Find the grid centers of the std_threshold failed grid centers
    std_threshold_failed_gridCenters = gridCenters_qualified_grids(current_std_threshold_failed_indices,1:2);
end


marker_size = 10;
RGB_triplet = [1 1 0]; 
legend_option = 1;
legend_name = 'Std threshold failed uncertain grids';
marker_type = 'o';
legend_position = [];
plot_std_threshold_failed_gridCenters = [std_threshold_failed_gridCenters, zeros(length(std_threshold_failed_gridCenters(:,1)),1)]; 
[LLA_data_std_threshold_failed_grids] = fcn_geometry_plotPointsinLLA(plot_std_threshold_failed_gridCenters,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],-1);
% fcn_geometry_plotPointsinLLA(plot_gridCenters_updated_original_unmapped_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% fcn_geometry_plotPointsinLLA(plot_gridCenters_non_drivable_grids,marker_size,RGB_triplet,[],legend_option,legend_name,[],[],[],[],fig_num);

geoplot(LLA_data_std_threshold_failed_grids(:,1), LLA_data_std_threshold_failed_grids(:,2), 'o','MarkerSize',10,'Color',RGB_triplet, 'LineWidth',2, 'DisplayName','Std threshold failed grids') 


if ~isempty(theta_threshold)
    % Find the indices of std threshold failed indices
    current_theta_threshold_failed_indices = current_uncertain_grid_numbers_in_mapped_grids(angle_btw_unit_normals_and_vertical(current_uncertain_grid_numbers_in_mapped_grids)>theta_threshold);

    % Find the grid centers of the theta_threshold failed grid centers
    theta_threshold_failed_gridCenters = gridCenters_qualified_grids(current_theta_threshold_failed_indices,1:2);
end


marker_size = 12.5;
RGB_triplet = [0 1 1]; 
legend_option = 1;
legend_name = 'Theta threshold failed uncertain grids';
marker_type = 'o';
legend_position = [];
plot_theta_threshold_failed_gridCenters = [theta_threshold_failed_gridCenters, zeros(length(theta_threshold_failed_gridCenters(:,1)),1)]; 
[LLA_data_theta_threshold_failed_grids] = fcn_geometry_plotPointsinLLA(plot_theta_threshold_failed_gridCenters,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],-1);

geoplot(LLA_data_theta_threshold_failed_grids(:,1), LLA_data_theta_threshold_failed_grids(:,2), 'o','MarkerSize',marker_size,'Color',RGB_triplet, 'LineWidth',2, 'DisplayName','Theta threshold failed grids') 

%% STEP 11: Find all boundary points: When uncertain grids are assumed as drivable

% Part1 - Find the boundary points of mapped and unmapped grids

% Revision History
% Funtionalized this code
% Added plotting options

% INPUTS - gridCenters_low_point_density,
% gridCenters_required_point_density, figure num
% OUTPUTS - X, Y, Z 

fig_num_mapped_unmapped = 767787; 

XYZ_matrix_mapped_grids = [gridCenters_qualified_grids(:,1:2) ones(length(gridCenters_qualified_grids(:,1)),1)]; 

XYZ_matrix_unmapped_grids = [gridCenters_unqualified_grids(:,1:2) zeros(length(gridCenters_unqualified_grids(:,1)),1)]; 

XYZ_matrix_mapped_unmapped_gridcenters = [XYZ_matrix_mapped_grids; XYZ_matrix_unmapped_grids]; 

% Find the unique elements
XYZ_matrix = unique(XYZ_matrix_mapped_unmapped_gridcenters,'rows'); 

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

%%%%%%%%%%%%%%---------------------------------------------------------------------------
% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
figure(fig_num_mapped_unmapped);
clf;
boundary_points_mapped_unmapped = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num_mapped_unmapped);

% assert(length(drivable_grids)>=1)
% assert(length(non_drivable_grids)>=1)
% % assert(length(unmapped_grids)==zeros(0,1))
% assert(length(gridCenters_mapped_grids(:,1))>=1)
% assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(boundary_points_mapped_unmapped)>=1)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundary points of drivable and non-drivable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 2 - Find boundary points of drivable and non-drivable grids

fig_num_drivable_non_drivable = 98898;

XYZ_matrix_drivable_grids = [gridCenters_drivable_grids(:,1:2) ones(length(gridCenters_drivable_grids(:,1)),1)]; 
XYZ_matrix_uncertain_grids = [gridCenters_uncertain_grids(:,1:2) ones(length(gridCenters_uncertain_grids(:,1)),1)];

XYZ_matrix_failed_grids = [gridCenters_failed_grids(:,1:2) zeros(length(gridCenters_failed_grids(:,1)),1)]; 

XYZ_matrix = [XYZ_matrix_drivable_grids; XYZ_matrix_uncertain_grids; XYZ_matrix_failed_grids]; 

% % XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)];
% XYZ_matrix = concatenate_gridCenters_drivable_non_drivable_grids;

XYZ_matrix = unique(XYZ_matrix,'rows'); 

% Reverse the matrix
XYZ_matrix_reversed = flipud(XYZ_matrix);

% Find unique rows based on the first two columns (X and Y)
[~,XYZ_matrix_indices] = unique(XYZ_matrix_reversed(:,1:2),'rows'); 

% Create the new matrix with only unique rows, keeping the last occurrence of each duplicate
XYZ_matrix_reversed = XYZ_matrix_reversed(XYZ_matrix_indices,:);  

% Reverse the matrix back to the original order
XYZ_matrix = flipud(XYZ_matrix_reversed);

XYZ_matrix = round(XYZ_matrix,4); 
% Given x_range and y_range
% x_range = min(XYZ_matrix(:,1)):grid_size:max(XYZ_matrix(:,1)); % min_x:gridSize:max_x 
% y_range = min(XYZ_matrix(:,2)):grid_size:max(XYZ_matrix(:,2)); % min_y:gridSize:max_y 

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

% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points

figure(fig_num_drivable_non_drivable)
clf;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num_drivable_non_drivable);

% assert(length(drivable_grids)>=1)
% assert(length(non_drivable_grids)>=1)
% % assert(length(unmapped_grids)==zeros(0,1))
% assert(length(gridCenters_mapped_grids(:,1))>=1)
% assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(boundary_points)>=1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finding true boundary points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 3 - Find true boundary points by eliminating the boundary points of
% mapped and unmapped from drivable and non-drivable boubndary points. 

fig_num_bd_pts_ENU = 1000; 

[members, id_x] = ismember(boundary_points,boundary_points_mapped_unmapped,'rows'); 

not_boundary_points = boundary_points(members,:);

true_boundary_points = boundary_points(members==0,:);

figure(fig_num_bd_pts_ENU)
clf;
% figure(10001)
hold on
plot(true_boundary_points(:,1), true_boundary_points(:,2), 'c.', 'MarkerSize',40)
plot(true_boundary_points(:,1), true_boundary_points(:,2), 'b.', 'MarkerSize',30)

fig_num = 452;
figure(fig_num); clf; 
% plot computed boundary points
marker_size = 25;
RGB_triplet = [0 1 1]; 
legend_option = 1;
legend_name = 'Computed boundary points';
legend_position = [];
marker_type = [];
plot_true_boundary_points = [true_boundary_points, zeros(length(true_boundary_points),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot computed boundary points

marker_size = 25;
RGB_triplet = [0 1 0]; 
legend_option = 1;
legend_name = 'Drivable grids';
legend_position = [];
marker_type = [];
plot_gridCenters_drivable_grids = [gridCenters_drivable_grids(:,1:2), zeros(length(gridCenters_drivable_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_drivable_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot grid centers
marker_size = 25;
RGB_triplet = [1 0 0]; 
legend_option = 1;
legend_name = 'Non-drivable grids';
legend_position = [];
marker_type = [];
plot_gridCenters_failed_grids = [gridCenters_failed_grids(:,1:2), zeros(length(gridCenters_failed_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_failed_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot grid centers
marker_size = 25;
RGB_triplet = [0 0 1];%[0.9290 0.6940 0.1250]; 
legend_option = 1;
legend_name = 'Uncertain grids';
legend_position = [];
marker_type = [];
plot_gridCenters_uncertain_grids = [gridCenters_uncertain_grids(:,1:2), zeros(length(gridCenters_uncertain_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_uncertain_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% marker_size = 10;
% RGB_triplet = [0 0 1]; 
% legend_option = 0;
% legend_name = 'Computed boundary points';
% legend_position = [];
% marker_type = [];
% [~] = fcn_geometry_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot driven path
marker_size = 25;
RGB_triplet = [0 0.5 0]; 
legend_option = 0;
legend_name = 'Driven path grids';
legend_position = [];
marker_type = [];
plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot driven path
marker_size = 10;
RGB_triplet = [0 0.8 0]; 
legend_option = 0;
legend_name = 'Driven path grids';
legend_position = [];
marker_type = [];
plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 11: Find boundary points: When uncertain grids are assumed as non drivable

% Part1 - Find the boundary points of mapped and unmapped grids

% Revision History
% Funtionalized this code
% Added plotting options

% INPUTS - gridCenters_low_point_density,
% gridCenters_required_point_density, figure num
% OUTPUTS - X, Y, Z 

fig_num_mapped_unmapped = 767787; 

XYZ_matrix_mapped_grids = [gridCenters_qualified_grids(:,1:2) ones(length(gridCenters_qualified_grids(:,1)),1)]; 

XYZ_matrix_unmapped_grids = [gridCenters_unqualified_grids(:,1:2) zeros(length(gridCenters_unqualified_grids(:,1)),1)]; 

XYZ_matrix_mapped_unmapped_gridcenters = [XYZ_matrix_mapped_grids; XYZ_matrix_unmapped_grids]; 

% Find the unique elements
XYZ_matrix = unique(XYZ_matrix_mapped_unmapped_gridcenters,'rows'); 

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

%%%%%%%%%%%%%%---------------------------------------------------------------------------
% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points
figure(fig_num_mapped_unmapped);
clf;
boundary_points_mapped_unmapped = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num_mapped_unmapped);

% assert(length(drivable_grids)>=1)
% assert(length(non_drivable_grids)>=1)
% % assert(length(unmapped_grids)==zeros(0,1))
% assert(length(gridCenters_mapped_grids(:,1))>=1)
% assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(boundary_points_mapped_unmapped)>=1)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundary points of drivable and non-drivable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 2 - Find boundary points of drivable and non-drivable grids

fig_num_drivable_non_drivable = 98898;

% XYZ_matrix = [gridCenters_mapped_grids(:,1:2) gridCenters_mapped_grids(:,4)];
XYZ_matrix = concatenate_gridCenters_drivable_non_drivable_grids;

XYZ_matrix = unique(XYZ_matrix,'rows'); 

% Reverse the matrix
XYZ_matrix_reversed = flipud(XYZ_matrix);

% Find unique rows based on the first two columns (X and Y)
[~,XYZ_matrix_indices] = unique(XYZ_matrix_reversed(:,1:2),'rows'); 

% Create the new matrix with only unique rows, keeping the last occurrence of each duplicate
XYZ_matrix_reversed = XYZ_matrix_reversed(XYZ_matrix_indices,:);  

% Reverse the matrix back to the original order
XYZ_matrix = flipud(XYZ_matrix_reversed);

XYZ_matrix = round(XYZ_matrix,4); 
% Given x_range and y_range
% x_range = min(XYZ_matrix(:,1)):grid_size:max(XYZ_matrix(:,1)); % min_x:gridSize:max_x 
% y_range = min(XYZ_matrix(:,2)):grid_size:max(XYZ_matrix(:,2)); % min_y:gridSize:max_y 

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

% x_limits = [min(x_range) max(x_range)];  
% y_limits = [min(y_range) max(y_range)]; 

x_limits = [];  
y_limits = []; 
% Calculate boundary points

figure(fig_num_drivable_non_drivable)
clf;
boundary_points = fcn_geometry_findBoundaryPoints(X,Y,Z,grid_size,x_limits,y_limits,fig_num_drivable_non_drivable);

% assert(length(drivable_grids)>=1)
% assert(length(non_drivable_grids)>=1)
% % assert(length(unmapped_grids)==zeros(0,1))
% assert(length(gridCenters_mapped_grids(:,1))>=1)
% assert(length(drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(non_drivable_grid_numbers_in_mapped_grids)>=1)
% assert(length(boundary_points)>=1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finding true boundary points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part 3 - Find true boundary points by eliminating the boundary points of
% mapped and unmapped from drivable and non-drivable boubndary points. 

fig_num_bd_pts_ENU = 1000; 

[members, id_x] = ismember(boundary_points,boundary_points_mapped_unmapped,'rows'); 

not_boundary_points = boundary_points(members,:);

true_boundary_points = boundary_points(members==0,:);

figure(fig_num_bd_pts_ENU)
clf;
% figure(10001)
hold on
plot(true_boundary_points(:,1), true_boundary_points(:,2), 'c.', 'MarkerSize',40)
plot(true_boundary_points(:,1), true_boundary_points(:,2), 'b.', 'MarkerSize',30)

fig_num = 434;
figure(fig_num); clf;

marker_size = 25;
RGB_triplet = [0 1 0]; 
legend_option = 1;
legend_name = 'Drivable grids';
legend_position = [];
marker_type = [];
plot_gridCenters_drivable_grids = [gridCenters_drivable_grids(:,1:2), zeros(length(gridCenters_drivable_grids(:,1)),1)]; 
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_drivable_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot grid centers
marker_size = 25;
RGB_triplet = [1 0 0]; 
legend_option = 1;
legend_name = 'NOT drivable grids';
legend_position = [];
marker_type = [];
plot_gridCenters_non_drivable_grids = [gridCenters_non_drivable_grids(:,1:2), zeros(length(gridCenters_non_drivable_grids(:,1)),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_non_drivable_grids,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


% plot computed boundary points
marker_size = 25;
RGB_triplet = [0 1 1]; 
legend_option = 1;
legend_name = 'Computed boundary points';
legend_position = [];
marker_type = [];
plot_true_boundary_points = [true_boundary_points, zeros(length(true_boundary_points),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% % plot computed boundary points
% marker_size = 10;
% RGB_triplet = [0 0 1]; 
% legend_option = 0;
% legend_name = 'Computed boundary points';
% legend_position = [];
% marker_type = [];
% [~] = fcn_geometry_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
% 

% plot driven path
marker_size = 25;
RGB_triplet = [0 0.5 0]; 
legend_option = 0;
legend_name = 'Driven path grids';
legend_position = [];
marker_type = [];
plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% plot driven path
marker_size = 10;
RGB_triplet = [0 0.8 0]; 
legend_option = 0;
legend_name = 'Driven path grids';
legend_position = [];
marker_type = [];
plot_gridCenters_driven_path = [gridCenters_driven_path, zeros(length(gridCenters_driven_path),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_gridCenters_driven_path,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% plot the boundary points
figure(7639);clf;

% plot computed boundary points
marker_size = 25;
RGB_triplet = [0 1 1]; 
legend_option = 1;
legend_name = 'Computed boundary points';
legend_position = [];
marker_type = [];
plot_true_boundary_points = [true_boundary_points, zeros(length(true_boundary_points),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_true_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


%% STEP 12: Find the nearest boundary points

fig_num = 7676;

% Find the nearest boundaries
[~, nearestBorderIndicies, nearestBorderXY] = fcn_geometry_findNearestBoundaryPoints(true_boundary_points, ...
    gridCenters_driven_path, grid_size, grid_boundaries, fig_num);

% Figure number
fig_num = 1224; 
figure(fig_num); clf;

% plot the nearest boundary points
marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 0;
legend_name = 'Nearest Boundary Points';
legend_position = [];
marker_type = []; 
nearest_boundary_points = [nearestBorderXY(:,1), nearestBorderXY(:,2), zeros(length(nearestBorderXY(:,1)),1)]; 

[~] = fcn_geometry_plotPointsinLLA(nearest_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

marker_size = 10;
RGB_triplet = [1 1 1]; 
legend_option = 0;
legend_name = 'Nearest Boundary Points';
legend_position = [];
marker_type = []; 
nearest_boundary_points = [nearestBorderXY(:,1), nearestBorderXY(:,2), zeros(length(nearestBorderXY(:,1)),1)]; 

[~] = fcn_geometry_plotPointsinLLA(nearest_boundary_points,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

%% STEP 13: Seperate the right and left boundaries from the nearest boundaries 

% boundaryLineNumber_start = scanLineNumber_start;%scanLineNumber_start - 8; 
% boundaryLineNumber_end = scanLineNumber_end;%scanLineNumber_end - 6; 

VehiclePose_current = VehiclePose(scanLineNumber_start:scanLineNumber_end,1:2);

% Get the number of rows 
rows_nearest_boundary_points = size(nearest_boundary_points, 1);
rows_VehiclePose_current= size(VehiclePose_current, 1);

% Calculate how many rows need to be added to VehiclePose_current_shifted
rows_to_add = rows_nearest_boundary_points - rows_VehiclePose_current;

if 0 <= rows_to_add
    
    % Method 1
    % % Update the VehiclePose_current with more rows
    % updated_VehiclePose_current = VehiclePose((scanLineNumber_start - rows_to_add):scanLineNumber_end,1:2);
    % 
    
    % Method 2
    offset_distance = grid_size/2; 
    
    additional_rows = VehiclePose((scanLineNumber_end-rows_to_add):scanLineNumber_end-1, 1:2) - [0, offset_distance];

    updated_VehiclePose_current = [VehiclePose_current; additional_rows];
    
    % Method 3
    % % Create the additional rows by repeating the last row of VehiclePose_current_shifted
    % additional_rows = repmat(VehiclePose_current(rows_VehiclePose_current, :), rows_to_add, 1);
    % 
    % % Concatenate the original VehiclePose_current_shifted and the additional rows
    % updated_VehiclePose_current = [VehiclePose_current; additional_rows];
    
else
    % % Find the total number of rows required 
    current_rows = length(VehiclePose_current(:,1)) + rows_to_add; 

    % Remove some rows from the original VehiclePose_current_shifted

    % Define the total number of indices
    total_indices = length(VehiclePose_current(:,1));

    % Define the number of indices to pick
    num_to_pick = current_rows;

    % Calculate the interval between indices
    interval = total_indices / num_to_pick;

    % Generate the indices based on the interval using vectorized operations
    selected_indices = round(((1:num_to_pick) - 0.5) * interval);

    % Ensure all indices are within the range [1, 50]
    selected_indices = min(max(selected_indices, 1), total_indices);

    % start_rows = (floor(abs(rows_to_add)/2)); 
    updated_VehiclePose_current = VehiclePose_current(selected_indices',:); 

end


vehicle_change_in_pose_XY = diff(updated_VehiclePose_current(:,1:2));

% Repeat the last value again, since diff removes one row. We want the same
% number of vectors as the number of points, and diff removed one point.
vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];

% Convert these to unit vectors
unit_vehicle_change_in_pose_XY = fcn_geometry_calcUnitVector(vehicle_change_in_pose_XY);


% Find orthogonal vetors by rotating by 90 degrees in the CCW direction
unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY*[0 1; -1 0];

% Transverse shift 
transverse_shift = 6*3.6576; 

% right transverse shifted points
right_transverse_points = updated_VehiclePose_current - transverse_shift*unit_ortho_vehicle_vectors_XY; 

% Find right nearest boundary points
right_boundary_points_nearest = [updated_VehiclePose_current;
    flipud(right_transverse_points);
    updated_VehiclePose_current(1,:)]; 

% Use inpolygon to find the indices of the right nearest boundary points
[boundary_points_right_indices,~] = inpolygon(nearest_boundary_points(:,1),nearest_boundary_points(:,2),right_boundary_points_nearest(:,1),right_boundary_points_nearest(:,2));

% nearest boundary points on the right
boundary_points_right = nearest_boundary_points(boundary_points_right_indices,1:2); 


% transverse shifted points
left_transverse_points = updated_VehiclePose_current + transverse_shift*unit_ortho_vehicle_vectors_XY; 


% Find left nearest boundary points
left_boundary_points_nearest = [updated_VehiclePose_current;
    flipud(left_transverse_points);
    updated_VehiclePose_current(1,:)]; 

% Use inpolygon to find the indices of the left nearest boundary points
[boundary_points_left_indices,~] = inpolygon(nearest_boundary_points(:,1),nearest_boundary_points(:,2),left_boundary_points_nearest(:,1),left_boundary_points_nearest(:,2));

% nearest boundary points on the left
boundary_points_left = nearest_boundary_points(boundary_points_left_indices,1:2); 


fig_num = 7876; 
figure(fig_num);clf; 
hold on
axis on
grid on 
xlabel('X[m]')
ylabel('Y[m]')
title('Right Points')

% plot(VehiclePose_current_shifted(:,1), VehiclePose_current_shifted(:,2),'.','Color',[0 0 0],'MarkerSize',30) 
plot(updated_VehiclePose_current(:,1), updated_VehiclePose_current(:,2),'.','Color',[0 0 0],'MarkerSize',30) 

% plot(VehiclePose_current_shifted(:,1), VehiclePose_current_shifted(:,2),'.','Color',[0 0 0],'MarkerSize',30) 
plot(VehiclePose_current(:,1), VehiclePose_current(:,2),'.','Color',[0.5 0.5 0.5],'MarkerSize',10) 

plot(nearest_boundary_points(:,1), nearest_boundary_points(:,2), 'c.', 'MarkerSize',40, 'DisplayName','Boundary points');
plot(nearest_boundary_points(:,1), nearest_boundary_points(:,2), 'b.', 'MarkerSize',30, 'DisplayName','Boundary points');

% plot(boundary_points_right(:,1), boundary_points_right(:,2), 'ro', 'MarkerSize',30, 'DisplayName','Boundary points');
% plot(boundary_points_left(:,1), boundary_points_left(:,2), 'bo', 'MarkerSize',30, 'DisplayName','Boundary points');

% quiver(...
%     updated_VehiclePose_current(:,1),updated_VehiclePose_current(:,2),...
%     unit_vehicle_change_in_pose_XY(1:2,1),unit_vehicle_change_in_pose_XY(:,2),'-','LineWidth',3,'Color',[0 1 0]);

% quiver(...
%     updated_VehiclePose_current(:,1),updated_VehiclePose_current(:,2),...
%     vector_from_vehicle_pose_to_boundary_points(:,1),vector_from_vehicle_pose_to_boundary_points(:,2),'-','LineWidth',3,'Color',[0 1 0]);


quiver(...
    updated_VehiclePose_current(:,1),updated_VehiclePose_current(:,2),...
    unit_ortho_vehicle_vectors_XY(:,1),unit_ortho_vehicle_vectors_XY(:,2),'-','LineWidth',3,'Color',[0 0 1]);

% plot the left shifted points
plot(right_transverse_points(:,1), right_transverse_points(:,2), 'k.', 'MarkerSize',30, 'DisplayName','Boundary points');

% Plot the left nearest points
plot(boundary_points_right(:,1), boundary_points_right(:,2), 'ro', 'MarkerSize',30, 'DisplayName','Boundary points');


% plot the left shifted points
plot(left_transverse_points(:,1), left_transverse_points(:,2), 'k.', 'MarkerSize',30, 'DisplayName','Boundary points');

% Plot the left nearest points
plot(boundary_points_left(:,1), boundary_points_left(:,2), 'bo', 'MarkerSize',30, 'DisplayName','Boundary points');

% % Boundary points on the right
% boundary_points_right_abs = nearest_boundary_points(abs(transverse_dist_boundary_points)<5,:);
% plot(boundary_points_right_abs(:,1), boundary_points_right_abs(:,2), 'go', 'MarkerSize',20, 'DisplayName','Boundary points');


%% Plot right boundary points

fig_num = 1224; 
figure(fig_num); 

marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 0;
legend_name = 'Right Boundary Points';
legend_position = [];
marker_type = []; 

plot_boundary_points_right = [boundary_points_right, zeros(length(boundary_points_right),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

marker_size = 12;
RGB_triplet = [1 0 0]; 
legend_option = 0;
legend_name = 'Right Boundary Points';
legend_position = [];
marker_type = []; 

plot_boundary_points_right = [boundary_points_right, zeros(length(boundary_points_right),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

% Plot right boundary points

fig_num = 1241; 
figure(fig_num); 

marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 0;
legend_name = 'Left Boundary Points';
legend_position = [];
marker_type = []; 

plot_boundary_points_left = [boundary_points_left, zeros(length(boundary_points_left),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_boundary_points_left,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

marker_size = 12;
RGB_triplet = [0 0 1]; 
legend_option = 0;
legend_name = 'Left Boundary Points';
legend_position = [];
marker_type = []; 

plot_boundary_points_left = [boundary_points_left, zeros(length(boundary_points_left),1)];
[~] = fcn_geometry_plotPointsinLLA(plot_boundary_points_left,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);




% %% STEP 13: Seperate the right and left boundaries from the nearest boundaries 
% 
% VehiclePose_current = VehiclePose(scanLineNumber_start:scanLineNumber_end,1:2);
% 
% % Get the number of rows 
% rows_nearest_boundary_points = size(nearest_boundary_points, 1);
% rows_VehiclePose_current= size(VehiclePose_current, 1);
% 
% % Calculate how many rows need to be added to VehiclePose_current_shifted
% rows_to_add = rows_nearest_boundary_points - rows_VehiclePose_current;
% 
% if 0 <= rows_to_add
% 
%     % Method 1
%     % % Update the VehiclePose_current with more rows
%     % updated_VehiclePose_current = VehiclePose((scanLineNumber_start - rows_to_add):scanLineNumber_end,1:2);
%     % 
% 
%     % Method 2
%     offset_distance = grid_size/2; 
% 
%     additional_rows = VehiclePose((scanLineNumber_end-rows_to_add):scanLineNumber_end-1, 1:2) - [0, offset_distance];
% 
%     updated_VehiclePose_current = [VehiclePose_current; additional_rows];
% 
%     % Method 3
%     % % Create the additional rows by repeating the last row of VehiclePose_current_shifted
%     % additional_rows = repmat(VehiclePose_current(rows_VehiclePose_current, :), rows_to_add, 1);
%     % 
%     % % Concatenate the original VehiclePose_current_shifted and the additional rows
%     % updated_VehiclePose_current = [VehiclePose_current; additional_rows];
% 
% else
%     % % Find the total number of rows required 
%     % current_rows = length(VehiclePose_current(:,1)) + rows_to_add; 
% 
%     % Remove some rows from the original VehiclePose_current_shifted 
%     updated_VehiclePose_current = VehiclePose_current(-rows_to_add:end-1,:); 
% 
% end
% 
% 
% vehicle_change_in_pose_XY = diff(updated_VehiclePose_current(:,1:2));
% 
% % Repeat the last value again, since diff removes one row. We want the same
% % number of vectors as the number of points, and diff removed one point.
% vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];
% 
% % Convert these to unit vectors
% unit_vehicle_change_in_pose_XY = fcn_geometry_calcUnitVector(vehicle_change_in_pose_XY);
% 
% shift = 5; 
% % Shift
% shift_distance = unit_vehicle_change_in_pose_XY*shift; 
% 
% % Shift the unit vector 
% unit_vehicle_change_in_pose_XY_shifted = unit_vehicle_change_in_pose_XY; %- shift_distance;
% 
% % Shift the vehicle pose
% updated_VehiclePose_current_shifted = updated_VehiclePose_current; %- shift_distance; 
% 
% % Calculate the vectors
% vector_from_vehicle_pose_to_boundary_points = nearest_boundary_points(:,1:2) - updated_VehiclePose_current_shifted;
% 
% % Find orthogonal vetors by rotating by 90 degrees in the CCW direction
% unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY_shifted*[0 1; -1 0];
% 
% % Calculate the transverse distance
% transverse_dist_boundary_points = sum(vector_from_vehicle_pose_to_boundary_points.*unit_ortho_vehicle_vectors_XY,2);
% 
% % Transverse distances of the right boundaries
% transverse_dist_right_boundary_points = transverse_dist_boundary_points(transverse_dist_boundary_points>0,:);
% 
% % Boundary points on the right
% boundary_points_right = nearest_boundary_points(transverse_dist_boundary_points>0,:);
% 
% figure(31543);clf; 
% hold on
% axis on
% grid on 
% xlabel('X[m]')
% ylabel('Y[m]')
% title('Right Points')
% 
% % plot(VehiclePose_current_shifted(:,1), VehiclePose_current_shifted(:,2),'.','Color',[0 0 0],'MarkerSize',30) 
% plot(updated_VehiclePose_current_shifted(:,1), updated_VehiclePose_current_shifted(:,2),'.','Color',[0 0 0],'MarkerSize',30) 
% 
% plot(nearest_boundary_points(:,1), nearest_boundary_points(:,2), 'c.', 'MarkerSize',40, 'DisplayName','Boundary points');
% plot(nearest_boundary_points(:,1), nearest_boundary_points(:,2), 'b.', 'MarkerSize',30, 'DisplayName','Boundary points');
% 
% plot(boundary_points_right(:,1), boundary_points_right(:,2), 'ro', 'MarkerSize',30, 'DisplayName','Boundary points');
% 
% quiver(...
%     updated_VehiclePose_current(1:2,1),updated_VehiclePose_current(1:2,2),...
%     unit_vehicle_change_in_pose_XY(1:2,1),unit_vehicle_change_in_pose_XY(1:2,2),'-','LineWidth',3,'Color',[0 1 0]);
% 
% quiver(...
%     updated_VehiclePose_current(1:2,1),updated_VehiclePose_current(1:2,2),...
%     unit_ortho_vehicle_vectors_XY(1:2,1),unit_ortho_vehicle_vectors_XY(1:2,2),'-','LineWidth',3,'Color',[0 0 1]);
% 
% 
% % Plot right boundary points
% 
% fig_num = 1224; 
% figure(fig_num); 
% 
% marker_size = 30;
% RGB_triplet = [0 1 1]; 
% legend_option = 0;
% legend_name = 'Right Boundary Points';
% legend_position = [];
% marker_type = []; 
% 
% 
% [~] = fcn_geometry_plotPointsinLLA(boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
% 
% marker_size = 12;
% RGB_triplet = [1 0 0]; 
% legend_option = 0;
% legend_name = 'Right Boundary Points';
% legend_position = [];
% marker_type = []; 
% 
% 
% [~] = fcn_geometry_plotPointsinLLA(boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);
% 
% 
% % %% Right boundary points - Using driven path grids
% % 
% % % VehiclePose_current = VehiclePose(scanLineNumber_start:scanLineNumber_end,1:2);
% % 
% % % Get the number of rows 
% % rows_nearest_boundary_points = size(nearest_boundary_points, 1);
% % rows_drivenPath_grids= size(gridCenters_driven_path, 1);
% % 
% % % Calculate how many rows need to be added to VehiclePose_current_shifted
% % rows_to_add = rows_nearest_boundary_points - rows_drivenPath_grids;
% % 
% % if 0 <= rows_to_add
% % 
% %     offset_distance = grid_size/2; 
% % 
% %     % Create the additional rows by repeating the last row of VehiclePose_current_shifted
% %     additional_rows = gridCenters_driven_path(end-rows_to_add:end-1,1:2) - [0, offset_distance];
% % 
% %     % Concatenate the original VehiclePose_current_shifted and the additional rows
% %     updated_drivenPath_grids = [gridCenters_driven_path; additional_rows];
% % 
% % else
% %     % Find the total number of rows required 
% %     current_rows = length(gridCenters_driven_path(:,1)) + rows_to_add; 
% % 
% %     % Remove some rows from the original VehiclePose_current_shifted 
% %     updated_drivenPath_grids = gridCenters_driven_path(1:current_rows,:); 
% % 
% % end
% % 
% % vehicle_change_in_drivenPath_grids_XY = diff(updated_drivenPath_grids(:,1:2));
% % 
% % % Repeat the last value again, since diff removes one row. We want the same
% % % number of vectors as the number of points, and diff removed one point.
% % vehicle_change_in_drivenPath_grids_XY = [vehicle_change_in_drivenPath_grids_XY; vehicle_change_in_drivenPath_grids_XY(end,:)];
% % 
% % % Convert these to unit vectors
% % unit_vehicle_change_in_drivenPath_grids_XY = fcn_geometry_calcUnitVector(vehicle_change_in_drivenPath_grids_XY);
% % 
% % % Calculate the vectors
% % vector_from_drivenPath_grids_to_boundary_points = nearest_boundary_points(:,1:2) - updated_drivenPath_grids;
% % 
% % % Find orthogonal vetors by rotating by 90 degrees in the CCW direction
% % unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_drivenPath_grids_XY*[0 1; -1 0];
% % 
% % % Calculate the transverse distance
% % transverse_dist_boundary_points = sum(vector_from_drivenPath_grids_to_boundary_points.*unit_ortho_vehicle_vectors_XY,2);
% % 
% % % Transverse distances of the right boundaries
% % transverse_dist_right_boundary_points = transverse_dist_boundary_points(transverse_dist_boundary_points>0,:);
% % 
% % % Boundary points on the right
% % boundary_points_right = nearest_boundary_points(transverse_dist_boundary_points>0,:);
% % 
% % figure(1224);clf; 
% % hold on
% % axis on
% % grid on 
% % xlabel('X[m]')
% % ylabel('Y[m]')
% % title('Right Points')
% % 
% % % plot(VehiclePose_current_shifted(:,1), VehiclePose_current_shifted(:,2),'.','Color',[0 0 0],'MarkerSize',30) 
% % plot(updated_drivenPath_grids(:,1), updated_drivenPath_grids(:,2),'.','Color',[0 0 0],'MarkerSize',30) 
% % 
% % plot(nearest_boundary_points(:,1), nearest_boundary_points(:,2), 'c.', 'MarkerSize',40, 'DisplayName','Boundary points');
% % plot(nearest_boundary_points(:,1), nearest_boundary_points(:,2), 'b.', 'MarkerSize',30, 'DisplayName','Boundary points');
% % 
% % plot(boundary_points_right(:,1), boundary_points_right(:,2), 'ro', 'MarkerSize',30, 'DisplayName','Boundary points');
% % 
% % quiver(...
% %     updated_VehiclePose_current(:,1),updated_VehiclePose_current(:,2),...
% %     unit_vehicle_change_in_pose_XY(:,1),unit_vehicle_change_in_pose_XY(:,2),'-','LineWidth',3,'Color',[0 1 0]);
% % 
% % quiver(...
% %     updated_VehiclePose_current(:,1),updated_VehiclePose_current(:,2),...
% %     unit_ortho_vehicle_vectors_XY(:,1),unit_ortho_vehicle_vectors_XY(:,2),'-','LineWidth',3,'Color',[0 0 1]);
% % 
% % 
% % 
% % 
% % 

%%

% fig_num = 34355;
% k = figure(fig_num); clf
% 
% for ith_mapped_grid = 39%;1:total_mapped_grids
%     % points = input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:);
%     % points = points(~isnan(points(:,1)),:);
%     [~, standard_deviation_in_z(ith_mapped_grid,:), ~, ~, ~, ~] =...
%     fcn_geometry_fitPlaneLinearRegression(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},:),fig_num);
%     % mean_z_of_mapped_grids(ith_mapped_grid,:) = mean(input_points(original_mapped_gridIndices_cell{ith_mapped_grid},3));
%     % z_diff_mapped_grids = abs(diff(mean_z_of_mapped_grids)); 
% end
% 
%    sgtitle('Plane fitting of a failed angle deviation grid')
