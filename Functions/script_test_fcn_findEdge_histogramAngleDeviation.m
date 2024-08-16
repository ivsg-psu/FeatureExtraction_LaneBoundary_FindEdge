% script_test_fcn_findEdge_histogramAngleDeviation
% Exercises the function: fcn_findEdge_histogramAngleDeviation

% Revision History:
% 2024_08_15 - Jiabao Zhao
% -- wrote the script test originally

%% Test 1 simple example
fig_num = 1;
angle_btw_unit_normals_and_vertical = [0.1, 0.15, 0.2]; % 3 points
angle_btw_unit_normals_and_vertical_driven_path = [0.12, 0.18, 0.22]; % 3 points
mean_angle_btw_unit_normals_and_vertical_driven_path = mean(angle_btw_unit_normals_and_vertical_driven_path);


[theta_threshold] = fcn_findEdge_histogramAngleDeviation(angle_btw_unit_normals_and_vertical, ...
    angle_btw_unit_normals_and_vertical_driven_path, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path,fig_num); 

assert(isequal((theta_threshold),0.1745))
% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% Test 2 no plotting
fig_num = 2;
angle_btw_unit_normals_and_vertical = [0.1, 0.15, 0.2]; % 3 points
angle_btw_unit_normals_and_vertical_driven_path = [0.12, 0.18, 0.22]; % 3 points
mean_angle_btw_unit_normals_and_vertical_driven_path = mean(angle_btw_unit_normals_and_vertical_driven_path);

[theta_threshold] = fcn_findEdge_histogramAngleDeviation(angle_btw_unit_normals_and_vertical, ...
    angle_btw_unit_normals_and_vertical_driven_path, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path,-1); 

assert(isequal((theta_threshold),0.1745))
% Check that the figure plotted
temp_h = figure(fig_num);
assert(isempty(get(temp_h,'Children')))
close(fig_num)

%% Test 3 random data points 
fig_num = 3;
% Example data generation
angle_btw_unit_normals_and_vertical = randn(1000,1) * 0.05; % 1000 points, mean 0, std 0.05 radians
angle_btw_unit_normals_and_vertical_driven_path = randn(100,1) * 0.03 + 0.1; % 100 points, mean 0.1, std 0.03 radians
mean_angle_btw_unit_normals_and_vertical_driven_path = mean(angle_btw_unit_normals_and_vertical_driven_path);

[theta_threshold] = fcn_findEdge_histogramAngleDeviation(angle_btw_unit_normals_and_vertical, ...
    angle_btw_unit_normals_and_vertical_driven_path, ...
    mean_angle_btw_unit_normals_and_vertical_driven_path,fig_num); 

assert(isequal((theta_threshold),0.1745))
% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))
