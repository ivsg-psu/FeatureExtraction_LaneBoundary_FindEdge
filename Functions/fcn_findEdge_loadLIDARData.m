function [VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData(varargin)
%% fcn_findEdge_loadLIDARData
% Loads LIDAR data from a specified file. Peforms loading of vehicle pose
% and LIDAR data with many conditional checks to confirm the file is
% present, the variables were not already loaded, the files are up-to-date,
% etc.
% 
% FORMAT:
%
%      [VehiclePose, LiDAR_Scan_ENU_Entire_Loop] = fcn_findEdge_loadLIDARData((test_date_string),(vehicle_pose_string), (LIDAR_file_string), (flag_load_all_data), (fig_num))
%
% INPUTS:     
%      
%      (OPTIONAL INPUTS)
%       
%      test_date_string = '2024_06_28'; % The date of testing. This defines
%      the folder where the data should be found within LargeData main
%      folder on the IVSG GitHubMirror wherein large data are stored.
%
%      vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file
%      containing VehiclePose
%
%      LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of
%      the file containing the LIDAR data
%
%      flag_load_all_data: set flat to 0 (default) to check before loading
%      anything, or 1 to force loading of data no matter what
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      VehiclePose: the position of the mapping vehicle during mapping
%
%      LiDAR_Scan_ENU_Entire_Loop: the points from the LIDAR scan 
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
% 
%       script_test_fcn_findEdge_loadLIDARData.m 
%  
%       for a full test suite.
%
% This function was written on 2024_08_05 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu

% Revision history
% 2024_07_17 - S Brennan
% -- wrote the code
% 2024_08_05 - Sean Brennan
% -- Created function by copying out of load script in Geometry library
% 2024_08_05 - Jiabao Zhao
% -- Functionalized this code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS");
    MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG = getenv("MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG); 
        flag_check_inputs  = str2double(MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end


%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0==flag_max_speed
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(0,5);

        % % Check the points input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points, '2column_of_numbers',[2 3]);
        %
        % % Check the transverse_tolerance input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);
        %
        % % Check the station_tolerance input is a positive single number
        % if ~isempty(station_tolerance)
        %     fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);
        % end
    end
end


% Does user want to specify station_tolerance?
test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
if (1<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        test_date_string = temp;
    end
end

% Does user want to specify station_tolerance?
vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
if (2<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        vehicle_pose_string = temp;
    end
end

% Does user want to specify station_tolerance?
LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
if (3<=nargin)
    temp = varargin{3};
    if ~isempty(temp)
        LIDAR_file_string = temp;
    end
end

% Does user want to specify flag_load_all_data?
flag_load_all_data = 0; % FORCE LOAD? Set this manually to 1 to FORCE load
if (4<=nargin)
    temp = varargin{4};
    if ~isempty(temp)
        flag_load_all_data = temp;
    end
end


% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (5<=nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end


%% Solve for the Maxs and Mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the file names
mat_loadFilename_vehiclePose = fullfile(cd,'LargeData',test_date_string, vehicle_pose_string);
mat_loadFilename_LIDARdata   = fullfile(cd,'LargeData',test_date_string, LIDAR_file_string);

persistent permanent_file_date PermanentVehiclePose PermanentLiDAR_Scan_ENU_Entire_Loop

% If the variable does not yet exist, load it
if isempty(PermanentVehiclePose)
    flag_load_all_data = 1;
end
if isempty(PermanentLiDAR_Scan_ENU_Entire_Loop)
    flag_load_all_data = 1;
end

% Does the data match in terms of dates?
if exist('permanent_file_date','var') || ~isempty(permanent_file_date)

    % Do the names match? if not, don't use the permanet ones - need to
    % reload!

    % Check the file's date of creation - if it doesn't match, the file has
    % been edited and needs to be reloaded

    % % BELOW is an example, commented out, of how to check date of this
    % % function:
    % st = dbstack;
    % this_function = st(1).file;
    % file_info = dir(which(this_function));

    file_info = dir(which(mat_loadFilename_vehiclePose));
    file_date = file_info.date;

    % Do we need to fill the file date? 
    if isempty(permanent_file_date)
        permanent_file_date = file_date;
    end

    % Does the file date from the first time we load match the date now? If
    % not, tell the user and plan to reload the data.
    if ~strcmp(file_date,permanent_file_date)
        fprintf(1,'\n\nComparing the current data file''s creation date to the date of the last load, they did not match.\n');
        fprintf(1,'\tCurrent file''s creation date: %s\n',file_date);
        fprintf(1,'\tCreation date of the file that was last loaded: %s\n',permanent_file_date)
        flag_load_all_data = 1;
    end

else
    flag_load_all_data = 1;
end

%% Load the data?
% If something indiactes that the data is wrong, the flag_load_all_data
% will have been set to 1 in the above steps

if 1==flag_load_all_data

    % Does the file exist?
    if exist(mat_loadFilename_vehiclePose,'file')
        fprintf(1,'Loading vehicle pose data from file: %s\n',mat_loadFilename_vehiclePose);
        load(mat_loadFilename_vehiclePose,'VehiclePose');
    else
        % File does not exist - need to warn the user
        error('Unable to find file: %s',mat_loadFilename_vehiclePose);
    end


    % Does the file exist?
    if exist(mat_loadFilename_LIDARdata,'file')
        fprintf(1,'Loading LIDAR scan data from file: %s\n',mat_loadFilename_LIDARdata);
        load(mat_loadFilename_LIDARdata,'LiDAR_Scan_ENU_Entire_Loop');
    else
        % File does not exist - need to load it
        error('Unable to find file: %s',mat_loadFilename_LIDARdata);
    end

    % Make sure both have same length
    if 1==flag_do_plots
        fprintf(1,'\n\nComparing the number of data in the vehicle pose data file to the number of data in the LIDAR scan file:\n');
        fprintf(1,'\tNumber of data in the vehicle pose file: %.0f\n',length(VehiclePose(:,1)));
        fprintf(1,'\tNumber of data in the LIDAR scan file: %.0f\n',length(LiDAR_Scan_ENU_Entire_Loop))
    end

    if length(VehiclePose(:,1)) ~= length(LiDAR_Scan_ENU_Entire_Loop)
        error('The number of data in the vehicle pose data file does not match the number of data in the LIDAR scan file.');
    end

    % Save the results into memory
    PermanentVehiclePose = VehiclePose;
    PermanentLiDAR_Scan_ENU_Entire_Loop = LiDAR_Scan_ENU_Entire_Loop;
else
    VehiclePose = PermanentVehiclePose;
    LiDAR_Scan_ENU_Entire_Loop = PermanentLiDAR_Scan_ENU_Entire_Loop;
    
    if 1==flag_do_plots
        fprintf(1,'Vehice pose and LIDAR data appears unchanged. Using existing variables within workspace.\n')
    end
end


%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plots
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end        

    clf;
    hold on;
    grid on;
    axis equal

    % plot of vehicle pose
    boundaryLineNumber_start = 1;
    boundaryLineNumber_end = length(VehiclePose);
    plot3(VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,1),VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,2),VehiclePose(boundaryLineNumber_start:boundaryLineNumber_end,3),'.','Color',[0 0 0],'MarkerSize',10,'LineWidth',3);
   

    title('ENU plot of vehicle trajectory');
    xlabel('East position [m]');
    ylabel('North position [m]');
    zlabel('Up position [m]');
    view(3)

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

    
end % Ends check if plotting

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§



