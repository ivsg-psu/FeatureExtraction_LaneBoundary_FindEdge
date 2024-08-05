% script_test_fcn_Laps_plotPointZoneDefinition.m.m
% Tests fcn_Laps_plotPointZoneDefinition.m
       
% Revision history:
%      2022_04_10
%      -- first write of the code
%      2022_07_23 S. Brennan, sbrennan@psu.edu
%      -- more examples

close all
clc

%% Call the plot command to show how it works. 
% Demonstrate that defaults work
zone_center = [1 2];
num_points = 3;
zone_radius = 5;
zone_definition = [zone_radius num_points zone_center];
fcn_Laps_plotPointZoneDefinition(zone_definition);


%% Show specification of plot style. 
% Demonstrate that defaults work
zone_center = [1 2];
num_points = 3;
zone_radius = 5;
zone_definition = [zone_radius num_points zone_center];
fcn_Laps_plotPointZoneDefinition(zone_definition,'r-');

%% Show specification of figure without plot style. 
% Demonstrate that defaults work
zone_center = [1 2];
num_points = 3;
zone_radius = 5;
zone_definition = [zone_radius num_points zone_center];
fig_num = 3;
fcn_Laps_plotPointZoneDefinition(zone_definition,[],fig_num);

%% Show zone, plot style, and fig number together
figure(11);
plot_style = 'g-';
zone_center = [1 2];
num_points = 3;
zone_radius = 5;
zone_definition = [zone_radius num_points zone_center];
fig_num = 4;
fcn_Laps_plotPointZoneDefinition(zone_definition,plot_style, fig_num);

%% Show that the fig_num option works
fig_num = 11121;
plot_style = 'g-.';
fcn_Laps_plotPointZoneDefinition(zone_definition,plot_style,fig_num);

