% script_test_findEdge_fitPlaneLinearRegression
% Exercises the function: findEdge_fitPlaneLinearRegression
% Revision history:
% 2021_05_24
% -- wrote the code
% -- revised from findEdge_fitSlopeInterceptNPoints
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it
close all;


%% Test 1: a basic test

fig_num = 1;
figure(fig_num);
clf;
rng(1823);

true_parameters = [ 0.1 0.2 3]';
points = randn(100,3);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

[parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_findEdge_fitPlaneLinearRegression([x y z],fig_num);


true_ABC = fcn_findEdge_calcUnitVector([-true_parameters(1) -true_parameters(2) 1]);
assert(isequal(round(true_parameters,4),round(parameters,4)));
assert(isequal(round(standard_deviation_in_z,4),round(0.00,4)));
assert(length(z_fit(:,1))==length(points(:,1)));
assert(isequal(round(true_parameters(3)*true_ABC,1),round(base_point,1)));
assert(isequal(round(true_ABC,4),round(unit_vector,4)));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),round(0.00,4)));

%% Test 2: add noise

fig_num = 1;
figure(fig_num);
clf;
rng(1823);

Npoints = 100;
true_parameters = [ 0.1 0.2 3]';
points = randn(Npoints,3);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

true_sigma = 0.1;
z = true_z + true_sigma*randn(Npoints,1);

[fitted_parameters, standard_deviation_in_z] = findEdge_fitPlaneLinearRegression([x y z],fig_num);

assert(isequal(round(true_parameters,1),round(fitted_parameters,1)));
assert(isequal(round(true_sigma,1),round(standard_deviation_in_z,1)));

%% Test of fast implementation mode 

rng(1823);

Npoints = 100;
true_parameters = [ 0.1 0.2 3]';
points = randn(Npoints,3);

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [root_point, unit_vector] = findEdge_fitPlaneLinearRegression(points, (fig_num));
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [root_point, unit_vector] = findEdge_fitPlaneLinearRegression(points, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of findEdge_fitPlaneLinearRegression:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);
%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [root_point, unit_vector] = findEdge_fitPlaneLinearRegression(points,fig_num);
    fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));
end
