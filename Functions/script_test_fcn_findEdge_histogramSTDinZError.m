% Script for fcn_findEdge_histogramSTDinZError
%
% FORMAT:
%   fcn_findEdge_histogramSTDinZError(standard_deviation_in_z,N_bins_stdz,std_in_z_driven_path,N_bins_std_drivenpath,mean_std_in_z_driven_path,std_threshold,(fig_num))
 
%
% REVISION HISTORY:
% 2024_08_15 - Aleksandr Goncharov
% -- Wrote the script
%%  TEST 1 - SIMPLE CASES - GENERATED DATA

fig_num=1111;
N_bins_stdz=20;
N_bins_std_drivenpath=10;
standard_deviation_in_z=rand(200,1);
std_in_z_driven_path=rand(50,1);
mean_std_in_z_driven_path=mean(std_in_z_driven_path);
std_threshold=0.1;

fcn_findEdge_histogramSTDinZError(standard_deviation_in_z,N_bins_stdz,std_in_z_driven_path,N_bins_std_drivenpath,mean_std_in_z_driven_path,std_threshold,(fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%%  TEST 2 - SIMPLE CASES - GENERATED DATA - MORE BINS


fig_num=2222;
N_bins_stdz=30;
N_bins_std_drivenpath=20;
standard_deviation_in_z=rand(200,1);
std_in_z_driven_path=rand(50,1);
mean_std_in_z_driven_path=mean(std_in_z_driven_path);
std_threshold=0.1;

fcn_findEdge_histogramSTDinZError(standard_deviation_in_z,N_bins_stdz,std_in_z_driven_path,N_bins_std_drivenpath,mean_std_in_z_driven_path,std_threshold,(fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%%  TEST 3 - SIMPLE CASES - GENERATED DATA - MORE BINS - BIGGER STD_THRESHOLD


fig_num=3333;
N_bins_stdz=30;
N_bins_std_drivenpath=20;
standard_deviation_in_z=rand(200,1);
std_in_z_driven_path=rand(50,1);
mean_std_in_z_driven_path=mean(std_in_z_driven_path);
std_threshold=0.5;

fcn_findEdge_histogramSTDinZError(standard_deviation_in_z,N_bins_stdz,std_in_z_driven_path,N_bins_std_drivenpath,mean_std_in_z_driven_path,std_threshold,(fig_num))

% Check that the figure plotted
temp_h = figure(fig_num);
assert(~isempty(get(temp_h,'Children')))

%% TEST 4 - SPEED MODE


N_bins_stdz=20;
N_bins_std_drivenpath=10;
standard_deviation_in_z=rand(200,1);
std_in_z_driven_path=rand(50,1);
mean_std_in_z_driven_path=mean(std_in_z_driven_path);
std_threshold=0.1;

% Speed Test Calculation
fig_num=[];
REPS=5; minTimeSlow=Inf;
tic;

%slow mode calculation - code copied from plotVehicleXYZ

for i=1:REPS
    tstart=tic;
    
fcn_findEdge_histogramSTDinZError(standard_deviation_in_z,N_bins_stdz,std_in_z_driven_path,N_bins_std_drivenpath,mean_std_in_z_driven_path,std_threshold,(fig_num))

    telapsed=toc(tstart);
    minTimeSlow=min(telapsed,minTimeSlow);
end
averageTimeSlow=toc/REPS;

%slow mode END

%Fast Mode Calculation
fig_num = -1;
minTimeFast = Inf;
tic;

for i=1:REPS
    tstart = tic;

fcn_findEdge_histogramSTDinZError(standard_deviation_in_z,N_bins_stdz,std_in_z_driven_path,N_bins_std_drivenpath,mean_std_in_z_driven_path,std_threshold,(fig_num))

    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

%Display Console Comparison
if 1==1
fprintf(1,'\n\nComparison of fcn_findEdge_plotLIDARXYZ without speed setting (slow) and with speed setting (fast):\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);
end

%Assertion on averageTime NOTE: Due to the variance, there is a chance that
%the assertion will fail.
assert(averageTimeFast<averageTimeSlow);

