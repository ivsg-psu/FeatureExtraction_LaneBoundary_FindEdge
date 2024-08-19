function [boundary_points_left, boundary_points_right] = fcn_findEdge_seperateLeftRightBoundaries...
    (VehiclePose,scanLineNumber_start,scanLineNumber_end,nearest_boundary_points, grid_size, transverse_shift, varargin)
%% fcn_findEdge_seperateLeftRightBoundaries    Seperate the right and left boundaries from the nearest boundaries
% 
% FORMAT:
%
%    [boundary_points_left, boundary_points_right] = fcn_findEdge_seperateLeftRightBoundaries...
%    (VehiclePose,station_1,station_2,nearest_boundary_points, grid_size,
%    transverse_shift, (fig_num)).
%
% INPUTS:
%      VehiclePose: the position of the mapping vehicle during mapping.    
%
%      station_1: the starting scan line of the LIDAR data
%
%      station_2: the end of scan line of the LIDAR data
%
%      nearest_boundary_points: the nearest boundary points close to the
%      driven path of of the mapping vehicle. 
%
%      grid_size: the size of the grid. 
%
%      transverse_shift: Transverse shift.
%
% (OPTIONAL INPUTS):
%
%
%      (fig_num): figure number, can be set to -1 for fast mode
%
% OUTPUTS: 
%
%      boundary_points_left: the boundary points on the left side of the of
%      the mapping vehcile 
%
%      boundary_points_right: the boundary points on the right side of the of
%      the mapping vehcile 
%
%
% DEPENDENCIES:
%
% EXAMPLES:
%   See the script:
%   
%   script_test_fcn_findEdge_seperateLeftRightBoundaries
%
%   for a full test suite
%
% This function was written by Jiabao Zhao on 2024_08_14
%
% Revision history
% 2024_08_31 - Aneesh Batchu
% -- wrote the code originally 
% 2024_08_14 - Jiabao Zhao
% -- Functionalized this code


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

if flag_max_speed == 0
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(6,7);
    end
end 


% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (6<=nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
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

% boundaryLineNumber_start = scanLineNumber_start;%scanLineNumber_start - 8; 
% boundaryLineNumber_end = scanLineNumber_end;%scanLineNumber_end - 6; 
% scanLineNumber_start = 1400; 
% scanLineNumber_end = 1450; 

VehiclePose_current = VehiclePose(scanLineNumber_start:scanLineNumber_start,1:2);

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

[unit_ortho_vehicle_vectors_XY, shifted_vehiclePose] = fcn_INTERNAL_findOrthAndShiftedVehiclePose(updated_VehiclePose_current); 

% Find the boundary points on the right
[right_transverse_points, boundary_points_right] = fcn_INTERNAL_findRightBoundaryPoints(nearest_boundary_points, updated_VehiclePose_current, shifted_vehiclePose, unit_ortho_vehicle_vectors_XY, transverse_shift);

% Find the boundary points on the left
[left_transverse_points, boundary_points_left] = fcn_INTERNAL_findLeftBoundaryPoints(nearest_boundary_points, updated_VehiclePose_current,shifted_vehiclePose, unit_ortho_vehicle_vectors_XY, transverse_shift); 


% % Find the difference of updated vehicle pose
% vehicle_change_in_pose_XY = diff(updated_VehiclePose_current(:,1:2));
% 
% % Repeat the last value again, since diff removes one row. We want the same
% % number of vectors as the number of points, and diff removed one point.
% vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];
% 
% % Convert these to unit vectors
% unit_vehicle_change_in_pose_XY = fcn_INTERNAL_calcUnitVector(vehicle_change_in_pose_XY);
% 
% % Find orthogonal vetors by rotating by 90 degrees in the CCW direction
% unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY*[0 1; -1 0];
% 
% shift = 5; 
% % Shift
% shift_distance = unit_vehicle_change_in_pose_XY*shift; 
% 
% % Shift the updated vehicle pose
% shifted_vehiclePose = updated_VehiclePose_current - shift_distance;
% 
% % right transverse shifted points
% right_transverse_points = updated_VehiclePose_current - transverse_shift*unit_ortho_vehicle_vectors_XY; 
% 
% % Find right nearest boundary points
% right_boundary_points_nearest = [shifted_vehiclePose;
%     flipud(right_transverse_points);
%     shifted_vehiclePose(1,:)]; 
% 
% % Use inpolygon to find the indices of the right nearest boundary points
% [boundary_points_right_indices,~] = inpolygon(nearest_boundary_points(:,1),nearest_boundary_points(:,2),right_boundary_points_nearest(:,1),right_boundary_points_nearest(:,2));
% 
% % nearest boundary points on the right
% boundary_points_right = nearest_boundary_points(boundary_points_right_indices,1:2); 



% % transverse shifted points
% left_transverse_points = updated_VehiclePose_current + transverse_shift*unit_ortho_vehicle_vectors_XY; 
% 
% % Find left nearest boundary points
% left_boundary_points_nearest = [shifted_vehiclePose;
%     flipud(left_transverse_points);
%     shifted_vehiclePose(1,:)]; 
% 
% % Use inpolygon to find the indices of the left nearest boundary points
% [boundary_points_left_indices,~] = inpolygon(nearest_boundary_points(:,1),nearest_boundary_points(:,2),left_boundary_points_nearest(:,1),left_boundary_points_nearest(:,2));
% 
% % nearest boundary points on the left
% boundary_points_left = nearest_boundary_points(boundary_points_left_indices,1:2); 

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
title('Right Points')

% plot(updated_VehiclePose_current(:,1), updated_VehiclePose_current(:,2),'.','Color',[0 0 0],'MarkerSize',30)
plot(shifted_vehiclePose(:,1), shifted_vehiclePose(:,2),'.','Color',[0 0 0],'MarkerSize',30)


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


% quiver(...
%     updated_VehiclePose_current(:,1),updated_VehiclePose_current(:,2),...
%     unit_ortho_vehicle_vectors_XY(:,1),unit_ortho_vehicle_vectors_XY(:,2),'-','LineWidth',3,'Color',[0 0 1]);

quiver(...
    shifted_vehiclePose(:,1),shifted_vehiclePose(:,2),...
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

function [unit_ortho_vehicle_vectors_XY, shifted_vehiclePose] = fcn_INTERNAL_findOrthAndShiftedVehiclePose(updated_VehiclePose_current)
% Find the difference of updated vehicle pose
vehicle_change_in_pose_XY = diff(updated_VehiclePose_current(:,1:2));

% Repeat the last value again, since diff removes one row. We want the same
% number of vectors as the number of points, and diff removed one point.
vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];

% Convert these to unit vectors
unit_vehicle_change_in_pose_XY = fcn_INTERNAL_calcUnitVector(vehicle_change_in_pose_XY);

% Find orthogonal vetors by rotating by 90 degrees in the CCW direction
unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY*[0 1; -1 0];

shift = 5; 
% Shift
shift_distance = unit_vehicle_change_in_pose_XY*shift; 

% Shift the updated vehicle pose
shifted_vehiclePose = updated_VehiclePose_current - shift_distance; 
end

function [right_transverse_points, boundary_points_right] = fcn_INTERNAL_findRightBoundaryPoints(nearest_boundary_points, updated_VehiclePose_current, shifted_vehiclePose, unit_ortho_vehicle_vectors_XY, transverse_shift)

% right transverse shifted points
right_transverse_points = updated_VehiclePose_current - transverse_shift*unit_ortho_vehicle_vectors_XY; 

% Find right nearest boundary points
right_boundary_points_nearest = [shifted_vehiclePose;
    flipud(right_transverse_points);
    shifted_vehiclePose(1,:)]; 

% Use inpolygon to find the indices of the right nearest boundary points
[boundary_points_right_indices,~] = inpolygon(nearest_boundary_points(:,1),nearest_boundary_points(:,2),right_boundary_points_nearest(:,1),right_boundary_points_nearest(:,2));

% nearest boundary points on the right
boundary_points_right = nearest_boundary_points(boundary_points_right_indices,1:2); 

end

function [left_transverse_points, boundary_points_left] = fcn_INTERNAL_findLeftBoundaryPoints(nearest_boundary_points, updated_VehiclePose_current,shifted_vehiclePose, unit_ortho_vehicle_vectors_XY, transverse_shift)
% transverse shifted points
left_transverse_points = updated_VehiclePose_current + transverse_shift*unit_ortho_vehicle_vectors_XY;


% Find left nearest boundary points
left_boundary_points_nearest = [shifted_vehiclePose;
    flipud(left_transverse_points);
    shifted_vehiclePose(1,:)];

% Use inpolygon to find the indices of the left nearest boundary points
[boundary_points_left_indices,~] = inpolygon(nearest_boundary_points(:,1),nearest_boundary_points(:,2),left_boundary_points_nearest(:,1),left_boundary_points_nearest(:,2));

% nearest boundary points on the left
boundary_points_left = nearest_boundary_points(boundary_points_left_indices,1:2);
end