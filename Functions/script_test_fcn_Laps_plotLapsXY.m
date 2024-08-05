% script_test_fcn_Laps_plotTraversalsXY.m
% Tests fcn_Laps_plotLapsXY
       
% Revision history:
%      2022_04_02
%      -- first write of the code

close all
clc


% Fill in some dummy data
laps = fcn_Laps_fillSampleLaps;
 

% Convert laps into traversals
for i_traveral = 1:length(laps)
    traversal = fcn_Path_convertPathToTraversalStructure(laps{i_traveral});
    data.traversal{i_traveral} = traversal;
end


%% Call the plot command to show how it works. First, put it into our figure
% to show that it will auto-label the axes and create a new figure (NOT
% figure 11 here) to plot the data.
figure(11);
fcn_Laps_plotLapsXY(data);


%% Next, specify the figure number to show that it will NOT auto-label the
% axes if figure is already given and it puts the plots into this figure.
fig_num = 12;
fcn_Laps_plotLapsXY(data,fig_num);
