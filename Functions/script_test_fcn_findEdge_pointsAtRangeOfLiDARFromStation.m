%% script_test_fcn_findEdge_pointsAtRangeOfLiDARFromStation
% This is a script for the function fcn_findEdge_pointsAtRangeOfLiDARFromStation
%
% FORMAT:
% [station1_minus_range_index, station2_plus_range_index]... 
% = fcn_findEdge_pointsAtRangeOfLiDARFromStation(VehiclePose,...
% starting_index,ending_index,(range))
% 
% Revision history:
%       2024_08_12 - Aleksandr Goncharov
%       -- Wrote the script for the function

%% TEST 1: Simple Case - Default Values

Test_Vehicle_Pose = [(-100:.1:100);(-100:.1:100)].';
Test_Starting_Index= 1000;
Test_Ending_Index=1200;
Range=[];

[a1,a2] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(Test_Vehicle_Pose,Test_Starting_Index,Test_Ending_Index,Range);

%% TEST 2: Simple Case - Range changed to 50, about half of default

Test_Vehicle_Pose = [(-100:.1:100);(-100:.1:100)].';
Test_Starting_Index= 1000;
Test_Ending_Index=1200;
Range=50;

[a1,a2] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(Test_Vehicle_Pose,Test_Starting_Index,Test_Ending_Index,Range);

%% TEST 3: Simple Case - Range changed to 50, about half of default. And Index moved to 1200-1400

Test_Vehicle_Pose = [(-100:.1:100);(-100:.1:100)].';
Test_Starting_Index= 1200;
Test_Ending_Index=1400;
Range=50;

[a1,a2] = fcn_findEdge_pointsAtRangeOfLiDARFromStation(Test_Vehicle_Pose,Test_Starting_Index,Test_Ending_Index,Range);

%%