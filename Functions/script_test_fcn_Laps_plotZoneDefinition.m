% script_test_fcn_Laps_plotZoneDefinition.m.m
% Tests fcn_Laps_plotZoneDefinition.m.m
       
% Revision history:
%      2022_07_23
%      -- first write of the code, based on fcn_Laps_plotSegmentZoneDefinition.m
close all
clc

%% Fill in some working data
segment_zone_definition = [1 2; 3 4]; % Starts at 1,2 and ends at 3,4

zone_center = [-1 2];
num_points = 3;
zone_radius = 2;
point_zone_definition = [zone_radius num_points zone_center];


%% Show default segment zone

fcn_Laps_plotZoneDefinition(segment_zone_definition);
axis([-5 5 -5 5]);

%% Show default point zone 
fcn_Laps_plotZoneDefinition(point_zone_definition);
axis([-5 5 -5 5]);

%% Show that can set the color segment zone
fcn_Laps_plotZoneDefinition(segment_zone_definition,'r');
axis([-5 5 -5 5]);

%% Show that can set the color point zone
fcn_Laps_plotZoneDefinition(point_zone_definition,'r');
axis([-5 5 -5 5]);


%% Specify zone, style, and figure
fig_num = 11121;
plot_style = 'g-';
fcn_Laps_plotZoneDefinition(point_zone_definition,plot_style,fig_num);
fcn_Laps_plotZoneDefinition(segment_zone_definition,plot_style,fig_num);
axis([-5 5 -5 5]);


%% Show that can set the figure without specifying plot color
% leave the plot style empty
fig_num = 222;
fcn_Laps_plotZoneDefinition(segment_zone_definition,[],fig_num);
fcn_Laps_plotZoneDefinition(point_zone_definition,[],fig_num);
axis([-5 5 -5 5]);

