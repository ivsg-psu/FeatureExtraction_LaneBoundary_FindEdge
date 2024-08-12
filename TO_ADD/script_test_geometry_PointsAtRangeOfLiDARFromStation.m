%% script_test_geometry_PointsAtRangeOfLiDARFromStation
%
% This code finds the index of a point that is range_of_LIDAR meters before
% station 1 and the index of a point that is range_of_LiDAR meters after
% station 2. 
%
% Later, use station1_minus_range_index and station2_plus_range_index as
% scanLinscanLineNumber_start and scanLinscanLineNumber_end, respectively
% in script_test_plot_data to generate concatenated matrices of LiDAR XYZ
% points, LiDAR intensities, LiDAR scan lines and rings

% Revision History
% 
% 2024_08_02 - Aneesh Batchu
% -- wrote the code originally

%% Main script

% Define the vehicle position
vehicle_positionsXY = VehiclePose(:,1:2); 

% Range of the LiDAR
range_of_LiDAR = 100; % This is the range of LiDAR we use

% Define the indices of stations S1 and S2
station_1 = 960; % Example index for S1
station_2 = 1010; % Example index for S2

% Calculate the differences between consecutive points
differences = diff(vehicle_positionsXY);

% Compute the Euclidean distances for each pair of consecutive points
distances = sqrt(sum(differences.^2, 2));

% Compute the cumulative sum of the distances
cumulative_distances = [0; cumsum(distances)];

% Find the cumulative distance of station 1 and station 2
station1_distance = cumulative_distances(station_1);
station2_distance = cumulative_distances(station_2);

% Find the index of the point before station 1 whose distance is
% approximately equal to the range of the LiDAR.
station1_minus_range_index = find(cumulative_distances <= station1_distance - range_of_LiDAR, 1, 'last');
station1_minus_range_point = vehicle_positionsXY(station1_minus_range_index, :);

% Find the index of the point after station 2 whose distance is
% approximately equal to the range of the LiDAR.
station2_plus_range_index = find(cumulative_distances >= station2_distance + range_of_LiDAR, 1, 'first');
station2_plus_range_point = vehicle_positionsXY(station2_plus_range_index, :);

% Display the indices and station points
disp(['The index of the scan line that is 100 meters before station 1 is ', num2str(station1_minus_range_index), '.']);
% disp(['The coordinates of the point that is 100 meters before S1 are [', num2str(station1_minus_range_point), '].']);
disp(['The index of the scan line that is 100 meters after station 2 is ', num2str(station2_plus_range_index), '.']);
% disp(['The coordinates of the point that is 100 meters after S2 are [', num2str(station2_plus_range_point), '].']);


