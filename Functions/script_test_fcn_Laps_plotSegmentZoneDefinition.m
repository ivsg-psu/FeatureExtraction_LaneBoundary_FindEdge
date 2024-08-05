% script_test_fcn_Laps_plotSegmentZoneDefinition.m
% Tests fcn_Laps_plotSegmentZoneDefinition.m
       
% Revision history:
%      2022_04_10
%      -- first write of the code
%      2022_07_23
%      -- more examples

close all
clc


%% Call the plot command to show how it works. 
% For it to use defaults
fcn_Laps_plotSegmentZoneDefinition([1 2; 3 4]);
axis([-5 5 -5 5]);

%% Show that can set the color
% First, put it into our figure
% to show that it will auto-label the axes and create a new figure (NOT
% figure 11 here) to plot the data.
figure(11);
plot_style = 'b.-';
fcn_Laps_plotSegmentZoneDefinition([1 2; 3 4],plot_style);
axis([-5 5 -5 5]);

%% Show that can set the figure without specifying plot color
% leave the plot style empyt
fig_num = 222;
fcn_Laps_plotSegmentZoneDefinition([1 2; 3 4],[],fig_num);
axis([-5 5 -5 5]);

%% Call the plot command to show how it works. First, put it into our figure
% to show that it will auto-label the axes and create a new figure (NOT
% figure 11 here) to plot the data.
fig_num = 11121;
plot_style = 'g.-.';
fcn_Laps_plotSegmentZoneDefinition([1 2; 3 4],plot_style,fig_num);

%% Call to show coloring
fig_num = 11122;
plot_style = 'r.-';
fcn_Laps_plotSegmentZoneDefinition([1 2; 3 4],plot_style,fig_num);