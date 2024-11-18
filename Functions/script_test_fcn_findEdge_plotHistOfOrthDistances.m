% script_test_fcn_findEdge_plotHistOfOrthDistances
% Exercises the function: fcn_findEdge_plotHistOfOrthDistances
%
% Revision History:
%
% 2024_10_24 - Aneesh Batchu
% - wrote the code originally

clc
close all

%% Example 1: Plot the histogram of computed boundary points of Sample Run 1 (Initial Run)

% Load the data
savefile = fullfile(pwd, 'Data', 'one_ComputedBoundaryPoints_Sample1.mat');
load(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')

% Figure number of histogram
fig_num1 = 11;
figure(fig_num1)
clf; 
% Figure number of LLA plot
fig_num2 = 12;
figure(fig_num2)
clf; 
title('Sample Run 1: Boundary points error plotted in color ')

% Figure number of ENU plot
fig_num = 13;
figure(fig_num)
clf; 

% Set by user - this determines the color scale
max_transverse_error = 2; 

% Number of points to plot (must be less than or equal to size(boundary_points_test_track_right,1))
nPoints = size(boundary_points_test_track_right,1); 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, max_transverse_error, fig_num1, fig_num2, nPoints, fig_num); 

assert(length(transverse_distances)>1); 

figure(fig_num1)

sigma_std = round(std(transverse_distances),2);
mu_mean = round(mean(transverse_distances),2);

text(mu_mean - 3.5*sigma_std, 8, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 3.5*sigma_std, 7, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
title('Sample Run 1: Histogram of orthogonal projections for 1 meter grid')
xlim([-8, 8])
ylim([0 20])

figure(fig_num2)
title('Sample Run 1: Boundary point errors plotted in color (1 meter grid)')

% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/one_Sample1.jpg')

%% Example 2: Plot the histogram of computed boundary points of Sample Run 2 (Run 1)

% Load the data
savefile = fullfile(pwd, 'Data', 'one_ComputedBoundaryPoints_Sample2.mat');
load(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')

% Figure number of histogram
fig_num1 = 111;
figure(fig_num1)
clf;

% Figure number of LLA plot
fig_num2 = 112;
figure(fig_num2)
clf;

% Figure number of ENU plot
fig_num = 113; 
figure(fig_num)
clf;

% Set by user - this determines the color scale
max_transverse_error = 2; 

% Number of points to plot (must be less than or equal to size(boundary_points_test_track_right,1))
nPoints = size(boundary_points_test_track_right,1); 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, max_transverse_error, fig_num1, fig_num2, nPoints, fig_num); 

assert(length(transverse_distances)>1); 

figure(fig_num1)

sigma_std = round(std(transverse_distances),2);
mu_mean = round(mean(transverse_distances),2);

text(mu_mean - 3.5*sigma_std, 8, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 3.5*sigma_std, 7, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
title('Sample Run 2: Histogram of orthogonal projections for 1 meter grid')
xlim([-8, 8])
ylim([0 20])

figure(fig_num2)
title('Sample Run 2: Boundary point errors plotted in color ((1 meter grid))')

% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/one_Sample2.jpg')

%% Example 3: Plot the histogram of computed boundary points of Sample Run 3 (Run 2)

% Load the data
savefile = fullfile(pwd, 'Data', 'one_ComputedBoundaryPoints_Sample3.mat');
load(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')

% % NOTE: the hand_labeled_path is backwards (CW). Need to fix it
% hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 211;

% Figure number of LLA plot
fig_num2 = 212;

% Figure number of ENU plot
fig_num = 213; 

% Set by user - this determines the color scale
max_transverse_error = 2; 

% Number of points to plot (must be less than or equal to size(boundary_points_test_track_right,1))
nPoints = size(boundary_points_test_track_right,1); 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, max_transverse_error, fig_num1, fig_num2, nPoints, fig_num); 

assert(length(transverse_distances)>1); 

figure(fig_num1)

sigma_std = round(std(transverse_distances),2);
mu_mean = round(mean(transverse_distances),2);

text(mu_mean - 3.5*sigma_std, 8, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 3.5*sigma_std, 7, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
title('Sample Run 3: Histogram of orthogonal projections for 1 meter grid')
xlim([-8, 8])
ylim([0 20])

figure(fig_num2)
title('Sample Run 3: Boundary point errors plotted in color (1 meter grid)')

% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/one_Sample3.jpg')
% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/one_Sample3.jpg')

%% Example 4: Plot the histogram of computed boundary points of Sample Run 4 (Run 3)

% Load the data
savefile = fullfile(pwd, 'Data', 'one_ComputedBoundaryPoints_Sample4.mat');
load(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')

% % NOTE: the hand_labeled_path is backwards (CW). Need to fix it
% hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 311;

% Figure number of LLA plot
fig_num2 = 312;

% Figure number of ENU plot
fig_num = 313; 

% Set by user - this determines the color scale
max_transverse_error = 2; 

% Number of points to plot (must be less than or equal to size(boundary_points_test_track_right,1))
nPoints = size(boundary_points_test_track_right,1); 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, max_transverse_error, fig_num1, fig_num2, nPoints, fig_num); 

assert(length(transverse_distances)>1);

figure(fig_num1)

sigma_std = round(std(transverse_distances),2);
mu_mean = round(mean(transverse_distances),2);

text(mu_mean - 3.5*sigma_std, 8, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 3.5*sigma_std, 7, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
title('Sample Run 4: Histogram of orthogonal projections for 1 meter grid')
xlim([-8, 8])
ylim([0 20])

figure(fig_num2)
title('Sample Run 4: Boundary point errors plotted in color (1 meter grid)')

% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/one_Sample4.jpg')

%% Example 5: Plot the histogram of computed boundary points of Sample Run 5 (Run 4)
% Load the data
savefile = fullfile(pwd, 'Data', 'one_ComputedBoundaryPoints_Sample5.mat');
load(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')

% % NOTE: the hand_labeled_path is backwards (CW). Need to fix it
% hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 411;

% Figure number of LLA plot
fig_num2 = 412;

% Figure number of ENU plot
fig_num = 413; 

% Set by user - this determines the color scale
max_transverse_error = 2; 

% Number of points to plot (must be less than or equal to size(boundary_points_test_track_right,1))
nPoints = size(boundary_points_test_track_right,1); 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, max_transverse_error, fig_num1, fig_num2, nPoints, fig_num); 

assert(length(transverse_distances)>1);

figure(fig_num1)

sigma_std = round(std(transverse_distances),2);
mu_mean = round(mean(transverse_distances),2);

text(mu_mean - 3.5*sigma_std, 8, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 3.5*sigma_std, 7, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
title('Sample Run 5: Histogram of orthogonal projections for 1 meter grid')
xlim([-8, 8])
ylim([0 20])

figure(fig_num2)
title('Sample Run 5: Boundary point errors plotted in color (1 meter grid)')

% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/one_Sample5.jpg')

%% Example 6: Plot the histogram of computed boundary points of Sample Run 6 (Run 5)

% Load the data
savefile = fullfile(pwd, 'Data', 'one_ComputedBoundaryPoints_Sample6.mat');
load(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')

% % NOTE: the hand_labeled_path is backwards (CW). Need to fix it
% hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 511;

% Figure number of LLA plot
fig_num2 = 512;

% Figure number of ENU plot
fig_num = 513; 

% Set by user - this determines the color scale
max_transverse_error = 2; 

% Number of points to plot (must be less than or equal to size(boundary_points_test_track_right,1))
nPoints = size(boundary_points_test_track_right,1); 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, max_transverse_error, fig_num1, fig_num2, nPoints, fig_num); 

assert(length(transverse_distances)>1);

figure(fig_num1)

sigma_std = round(std(transverse_distances),2);
mu_mean = round(mean(transverse_distances),2);

text(mu_mean - 3.5*sigma_std, 8, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 3.5*sigma_std, 7, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
title('Sample Run 6: Histogram of orthogonal projections for 1 meter grid')
xlim([-8, 8])
ylim([0 20])

figure(fig_num2)
title('Sample Run 6: Boundary point errors plotted in color (1 meter grid)')

% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/one_Sample6.jpg')

%% Example 7: Plot the histogram of computed boundary points of Sample Run 7 (Run 6)

% Load the data
savefile = fullfile(pwd, 'Data', 'one_ComputedBoundaryPoints_Sample7.mat');
load(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')

% % NOTE: the hand_labeled_path is backwards (CW). Need to fix it
% hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 611;

% Figure number of LLA plot
fig_num2 = 612;

% Figure number of ENU plot
fig_num = 613; 

% Set by user - this determines the color scale
max_transverse_error = 2; 

% Number of points to plot (must be less than or equal to size(boundary_points_test_track_right,1))
nPoints = size(boundary_points_test_track_right,1); 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, max_transverse_error, fig_num1, fig_num2, nPoints, fig_num); 

assert(length(transverse_distances)>1);

figure(fig_num1)

sigma_std = round(std(transverse_distances),2);
mu_mean = round(mean(transverse_distances),2);

text(mu_mean - 3*sigma_std, 8, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 3*sigma_std, 7, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
title('Sample Run 7: Histogram of orthogonal projections for 1 meter grid')
xlim([-8, 8])
ylim([0 20])

figure(fig_num2)
title('Sample Run 7: Boundary point errors plotted in color (1 meter grid)')

% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/one_Sample7.jpg')

%% Example 7 - (0.3 grid size)Plot the histogram of computed boundary points of Sample Run 1 (Initial Run)

% Load the data
savefile = fullfile(pwd, 'Data', 'point3_ComputedBoundaryPoints_Sample1.mat');
load(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')
% 
% % NOTE: the hand_labeled_path is backwards (CW). Need to fix it
% hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 711;

% Figure number of LLA plot
fig_num2 = 712;

% Figure number of ENU plot
fig_num = 713; 

% Set by user - this determines the color scale
max_transverse_error = 2; 

% Number of points to plot (must be less than or equal to size(boundary_points_test_track_right,1))
nPoints = 1397; 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, max_transverse_error, fig_num1, fig_num2, nPoints, fig_num); 

figure(fig_num1);

% axis equal
title( 'Histogram of all orthogonal distance projections for grid size = 0.3 m');
xlabel('Orthogonal Distance (m)')

assert(length(transverse_distances)>1);

figure(fig_num1)

sigma_std = round(std(transverse_distances),2);
mu_mean = round(mean(transverse_distances),2);

text(mu_mean - 4*sigma_std, 8, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 4*sigma_std, 7, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
title('Sample Run 1: Histogram of orthogonal projections for 0.3 meter grid')
xlim([-8, 8])
ylim([0 20])

figure(fig_num2)
title('Sample Run 1: Boundary point errors plotted in color (0.3 meter grid)')

% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/point3_Sample1.jpg')

%% Example 8 - (0.3 grid size)Plot the histogram of computed boundary points of Sample Run 2 (Run 1)

% Load the data
savefile = fullfile(pwd, 'Data', 'point3_ComputedBoundaryPoints_Sample2.mat');
load(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')
% 
% % NOTE: the hand_labeled_path is backwards (CW). Need to fix it
% hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 811;

% Figure number of LLA plot
fig_num2 = 812;

% Figure number of ENU plot
fig_num = 813; 

% Set by user - this determines the color scale
max_transverse_error = 2; 

% Number of points to plot (must be less than or equal to size(boundary_points_test_track_right,1))
nPoints = 831; 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, max_transverse_error, fig_num1, fig_num2, nPoints, fig_num); 

figure(fig_num1);

% axis equal
title( 'Histogram of all orthogonal distance projections for grid size = 0.3 m');
xlabel('Orthogonal Distance (m)')

assert(length(transverse_distances)>1);

figure(fig_num1)

sigma_std = round(std(transverse_distances),2);
mu_mean = round(mean(transverse_distances),2);

text(mu_mean - 4*sigma_std, 8, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 4*sigma_std, 7, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
title('Sample Run 2: Histogram of orthogonal projections for 0.3 meter grid')
xlim([-8, 8])
ylim([0 20])

figure(fig_num2)
title('Sample Run 2: Boundary point errors plotted in color (0.3 meter grid)')

% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/point3_Sample2.jpg')

%% Example 9 - (0.3 grid size)Plot the histogram of computed boundary points of Sample Run 3 (Run 2)

% Load the data
savefile = fullfile(pwd, 'Data', 'point3_ComputedBoundaryPoints_Sample3.mat');
load(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')
% 
% % NOTE: the hand_labeled_path is backwards (CW). Need to fix it
% hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 911;

% Figure number of LLA plot
fig_num2 = 912;

% Figure number of ENU plot
fig_num = 913; 

% Set by user - this determines the color scale
max_transverse_error = 2; 

% Number of points to plot (must be less than or equal to size(boundary_points_test_track_right,1))
nPoints = 1009; 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, max_transverse_error, fig_num1, fig_num2, nPoints, fig_num); 

figure(fig_num1);

% axis equal
title( 'Histogram of all orthogonal distance projections for grid size = 0.3 m');
xlabel('Orthogonal Distance (m)')

assert(length(transverse_distances)>1);

figure(fig_num1)

sigma_std = round(std(transverse_distances),2);
mu_mean = round(mean(transverse_distances),2);

text(mu_mean - 4*sigma_std, 8, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 4*sigma_std, 7, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
title('Sample Run 3: Histogram of orthogonal projections for 0.3 meter grid')
xlim([-8, 8])
ylim([0 20])

figure(fig_num2)
title('Sample Run 3: Boundary point errors plotted in color (0.3 meter grid)')

% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/point3_Sample3.jpg')

%% Example 10 - (0.3 grid size)Plot the histogram of computed boundary points of Sample Run 4 (Run 3)

% Load the data
savefile = fullfile(pwd, 'Data', 'point3_ComputedBoundaryPoints_Sample4.mat');
load(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')
% 
% % NOTE: the hand_labeled_path is backwards (CW). Need to fix it
% hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 1011;

% Figure number of LLA plot
fig_num2 = 1012;

% Figure number of ENU plot
fig_num = 1013; 

% Set by user - this determines the color scale
max_transverse_error = 2; 

% Number of points to plot (must be less than or equal to size(boundary_points_test_track_right,1))
nPoints = 1103; 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, max_transverse_error, fig_num1, fig_num2, nPoints, fig_num); 

figure(fig_num1);

% axis equal
title( 'Histogram of all orthogonal distance projections for grid size = 0.3 m');
xlabel('Orthogonal Distance (m)')

% assert(length(transverse_distances)>1);

figure(fig_num1)

sigma_std = round(std(transverse_distances),2);
mu_mean = round(mean(transverse_distances),2);

text(mu_mean - 4*sigma_std, 8, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 4*sigma_std, 7, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
title('Sample Run 4: Histogram of orthogonal projections for 0.3 meter grid')
xlim([-8, 8])
ylim([0 20])

figure(fig_num2)
title('Sample Run 4: Boundary point errors plotted in color (0.3 meter grid)')

% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/point3_Sample4.jpg')


%% Example 10 - (0.3 grid size)Plot the histogram of computed boundary points of Sample Run 5 (Run 3)

% Load the data
savefile = fullfile(pwd, 'Data', 'point3_ComputedBoundaryPoints_Sample5.mat');
load(savefile, 'boundary_points_test_track_right', 'hand_labeled_path')
% 
% % NOTE: the hand_labeled_path is backwards (CW). Need to fix it
% hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 1111;

% Figure number of LLA plot
fig_num2 = 1112;

% Figure number of ENU plot
fig_num = 1113; 

% Set by user - this determines the color scale
max_transverse_error = 2; 

% Number of points to plot (must be less than or equal to size(boundary_points_test_track_right,1))
nPoints = 992; 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, max_transverse_error, fig_num1, fig_num2, nPoints, fig_num); 

figure(fig_num1);

% axis equal
title( 'Histogram of all orthogonal distance projections for grid size = 0.3 m');
xlabel('Orthogonal Distance (m)')

% assert(length(transverse_distances)>1);

figure(fig_num1)

sigma_std = round(std(transverse_distances),2);
mu_mean = round(mean(transverse_distances),2);

text(mu_mean - 4*sigma_std, 8, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 4*sigma_std, 7, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
title('Sample Run 5: Histogram of orthogonal projections for 0.3 meter grid')
xlim([-8, 8])
ylim([0 20])

figure(fig_num2)
title('Sample Run 5: Boundary point errors plotted in color (0.3 meter grid)')

% saveas(figure(fig_num2), '/Users/aneeshbatchu/Documents/PennState_Semesters/05_Fall 2024/FA24_IVSG/Point3_one_CombinedHistPlots/point3_Sample5.jpg')

