% script_test_fcn_findEdge_plotHistOfOrthDistances
% Exercises the function: fcn_findEdge_plotHistOfOrthDistances
%
% Revision History:
%
% 2024_10_24 - Aneesh Batchu
% - wrote the code originally

clc
close all

%% Example 1: Plot the histogram of computed boundary points of Initial Run

% Load the data
savefile = fullfile(pwd, 'Data', 'ExampleData_findOrthoError_1.mat');
load(savefile, 'computed_boundary_path', 'hand_labeled_path')

% NOTE: the hand_labeled_path is backwards (CW). Need to fix it
hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 111;

% Figure number of LLA plot
fig_num2 = 112;

% Figure number of ENU plot
fig_num = 113; 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, computed_boundary_path, fig_num1, fig_num2, fig_num); 

assert(length(transverse_distances)>1); 

%% Example 2: Plot the histogram of computed boundary points of Run 1

% Load the data
savefile = fullfile(pwd, 'Data', 'ExampleData_findOrthoError_2.mat');
load(savefile, 'computed_boundary_path', 'hand_labeled_path')

% NOTE: the hand_labeled_path is backwards (CW). Need to fix it
hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 211;

% Figure number of LLA plot
fig_num2 = 212;

% Figure number of ENU plot
fig_num = 213; 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, computed_boundary_path, fig_num1, fig_num2, fig_num); 

assert(length(transverse_distances)>1); 

%% Example 3: Plot the histogram of computed boundary points of Run 2

% Load the data
savefile = fullfile(pwd, 'Data', 'ExampleData_findOrthoError_3.mat');
load(savefile, 'computed_boundary_path', 'hand_labeled_path')

% NOTE: the hand_labeled_path is backwards (CW). Need to fix it
hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 311;

% Figure number of LLA plot
fig_num2 = 312;

% Figure number of ENU plot
fig_num = 313; 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, computed_boundary_path, fig_num1, fig_num2, fig_num); 

assert(length(transverse_distances)>1);

%% Example 4: Plot the histogram of computed boundary points of Run 3

% Load the data
savefile = fullfile(pwd, 'Data', 'ExampleData_findOrthoError_4.mat');
load(savefile, 'computed_boundary_path', 'hand_labeled_path')

% NOTE: the hand_labeled_path is backwards (CW). Need to fix it
hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 411;

% Figure number of LLA plot
fig_num2 = 412;

% Figure number of ENU plot
fig_num = 413; 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, computed_boundary_path, fig_num1, fig_num2, fig_num); 

assert(length(transverse_distances)>1);

%% Example 5: Plot the histogram of computed boundary points of Run 4

% Load the data
savefile = fullfile(pwd, 'Data', 'ExampleData_findOrthoError_5.mat');
load(savefile, 'computed_boundary_path', 'hand_labeled_path')

% NOTE: the hand_labeled_path is backwards (CW). Need to fix it
hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 511;

% Figure number of LLA plot
fig_num2 = 512;

% Figure number of ENU plot
fig_num = 513; 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, computed_boundary_path, fig_num1, fig_num2, fig_num); 

assert(length(transverse_distances)>1);

%% Example 6: Plot the histogram of computed boundary points of Run 5

% Load the data
savefile = fullfile(pwd, 'Data', 'ExampleData_findOrthoError_6.mat');
load(savefile, 'computed_boundary_path', 'hand_labeled_path')

% NOTE: the hand_labeled_path is backwards (CW). Need to fix it
hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 611;

% Figure number of LLA plot
fig_num2 = 612;

% Figure number of ENU plot
fig_num = 613; 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, computed_boundary_path, fig_num1, fig_num2, fig_num); 

assert(length(transverse_distances)>1);

%% Point 3 grid size

% % Load the data
% savefile = fullfile(pwd, 'Data', 'ExampleData_findOrthoError_7.mat');
% load(savefile, 'computed_boundary_path', 'hand_labeled_path')

% NOTE: the hand_labeled_path is backwards (CW). Need to fix it
hand_labeled_path = flipud(hand_labeled_path);

% Figure number of histogram
fig_num1 = 111;

% Figure number of LLA plot
fig_num2 = 112;

% Figure number of ENU plot
fig_num = 113; 

transverse_distances = fcn_findEdge_plotHistOfOrthDistances(hand_labeled_path, boundary_points_test_track_right, fig_num1, fig_num2, fig_num); 

figure(fig_num1);

% axis equal
title( 'Histogram of all orthogonal distance projections for grid size = 0.3 m');
xlabel('Orthogonal Distance (m)')

assert(length(transverse_distances)>1); 






