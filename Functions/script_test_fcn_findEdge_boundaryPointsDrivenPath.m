%% script_test_fcn_findEdge_boundaryPointsDrivenPath
% Exercises the function: fcn_findEdge_boundaryPointsDrivenPath

% Revision History:
% 2024_08_06 - Aleksandr Goncharov
% -- Wrote the testing script

%% TEST 1: Simple Input

VehiclePose=[linspace(10,100,1000).',linspace(10,100,1000).',linspace(-11,0,1000).'];
scanLineRange=[300 400];
fig_num=1000;
shift=[];
figure(fig_num);
clf;

hold on
%Main Function

[left_boundary_points,right_boundary_points,boundaryLineNumber_start,boundaryLineNumber_end]= fcn_findEdge_boundaryPointsDrivenPath(VehiclePose,scanLineRange,shift,(fig_num));


%Additional Plotting for Step 2
format1=sprintf(' ''.'',''Color'',[0 0 0],''MarkerSize'',30,''LineWidth'',3');
format2=sprintf(' ''.'',''Color'',[1 1 0],''MarkerSize'',30,''LineWidth'',3');

format3=sprintf(' ''.'',''Color'',[0 1 0],''MarkerSize'',30,''LineWidth'',3');
format4=sprintf(' ''.'',''Color'',[0 0 1],''MarkerSize'',30,''LineWidth'',3');

fcn_findEdge_plotLIDARXYZ(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format1,fig_num)
fcn_findEdge_plotLIDARXYZ(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format2,fig_num)

fcn_findEdge_plotLIDARXYZ(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format3,fig_num)
fcn_findEdge_plotLIDARXYZ(right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format4,fig_num)
%
xlabel('East position [m]');
ylabel('North position [m]');
zlabel('Up position [m]');
view(3)
hold off

%% TEST 2: Simple Input - Changed Format (All Green)

VehiclePose=[linspace(10,100,1000).',linspace(10,100,1000).',linspace(-11,-10,1000).'];
scanLineRange=[300 400];
fig_num=2222;
shift=[];
figure(fig_num);
clf;

hold on
%Main Function

[left_boundary_points,right_boundary_points,boundaryLineNumber_start,boundaryLineNumber_end]= fcn_findEdge_boundaryPointsDrivenPath(VehiclePose,scanLineRange,shift,fig_num);

%Additional Plotting for Step 2
format1=sprintf(' ''.'',''Color'',[0 1 0],''MarkerSize'',30,''LineWidth'',3');
format2=sprintf(' ''.'',''Color'',[0 1 0],''MarkerSize'',30,''LineWidth'',3');

format3=sprintf(' ''.'',''Color'',[0 1 0],''MarkerSize'',30,''LineWidth'',3');
format4=sprintf(' ''.'',''Color'',[0 1 0],''MarkerSize'',30,''LineWidth'',3');

fcn_findEdge_plotLIDARXYZ(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format1,fig_num)
fcn_findEdge_plotLIDARXYZ(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format2,fig_num)

fcn_findEdge_plotLIDARXYZ(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format3,fig_num)
fcn_findEdge_plotLIDARXYZ(right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format4,fig_num)

xlabel('East position [m]');
ylabel('North position [m]');
zlabel('Up position [m]');
view(3)
hold off

%% TEST 3: Simple Inputs - All Lines

VehiclePose=[linspace(10,100,1000).',linspace(10,100,1000).',linspace(-11,0,1000).'];
scanLineRange=[1 length(VehiclePose)];
fig_num=3333;
shift=[];
figure(fig_num);
clf;

hold on
%Main Function

[left_boundary_points,right_boundary_points,boundaryLineNumber_start,boundaryLineNumber_end]= fcn_findEdge_boundaryPointsDrivenPath(VehiclePose,scanLineRange,shift,(fig_num));

%Additional Plotting for Step 2
format1=sprintf(' ''.'',''Color'',[0 0 0],''MarkerSize'',30,''LineWidth'',3');
format2=sprintf(' ''.'',''Color'',[1 1 0],''MarkerSize'',30,''LineWidth'',3');

format3=sprintf(' ''.'',''Color'',[0 1 0],''MarkerSize'',30,''LineWidth'',3');
format4=sprintf(' ''.'',''Color'',[0 0 1],''MarkerSize'',30,''LineWidth'',3');

fcn_findEdge_plotLIDARXYZ(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format1,fig_num)
fcn_findEdge_plotLIDARXYZ(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format2,fig_num)

fcn_findEdge_plotLIDARXYZ(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format3,fig_num)
fcn_findEdge_plotLIDARXYZ(right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format4,fig_num)

xlabel('East position [m]');
ylabel('North position [m]');
zlabel('Up position [m]');
view(3)
hold off

%% TEST 4: Simple Inputs - All Lines - Shift 10

VehiclePose=[linspace(10,100,1000).',linspace(10,100,1000).',linspace(-11,0,1000).'];
scanLineRange=[1 length(VehiclePose)];
fig_num=4444;
shift=10;
figure(fig_num);
clf;

hold on
%Main Function

[left_boundary_points,right_boundary_points,boundaryLineNumber_start,boundaryLineNumber_end]= fcn_findEdge_boundaryPointsDrivenPath(VehiclePose,scanLineRange,shift,(fig_num));

%Additional Plotting for Step 2
format1=sprintf(' ''.'',''Color'',[0 0 0],''MarkerSize'',30,''LineWidth'',3');
format2=sprintf(' ''.'',''Color'',[1 1 0],''MarkerSize'',30,''LineWidth'',3');

format3=sprintf(' ''.'',''Color'',[0 1 0],''MarkerSize'',30,''LineWidth'',3');
format4=sprintf(' ''.'',''Color'',[0 0 1],''MarkerSize'',30,''LineWidth'',3');

fcn_findEdge_plotLIDARXYZ(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format1,fig_num)
fcn_findEdge_plotLIDARXYZ(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format2,fig_num)

fcn_findEdge_plotLIDARXYZ(left_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format3,fig_num)
fcn_findEdge_plotLIDARXYZ(right_boundary_points(boundaryLineNumber_start:boundaryLineNumber_end,1:3),[],[],[],format4,fig_num)

xlabel('East position [m]');
ylabel('North position [m]');
zlabel('Up position [m]');
view(3)
hold off

