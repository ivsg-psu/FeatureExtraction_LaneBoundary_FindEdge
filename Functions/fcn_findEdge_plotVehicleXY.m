function fcn_findEdge_plotVehicleXY(VehiclePose, varargin)
%% fcn_findEdge_plotVehicleXY
% Plots the vehicle's XY position.
% 
% FORMAT:
%
%      fcn_findEdge_loadLIDARData(VehiclePose, (plot_str), (fig_num))
%
% INPUTS:  
%
%      VehiclePose: the position of the mapping vehicle during mapping
%      
%      (OPTIONAL INPUTS)
%
%      format: A format string, e.g. 'b-', that dictates the plot style or
%      a color vector, e.g. [1 0 0.23], that dictates the line color.
%       
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
% 
%       script_test_fcn_findEdge_plotVehicleXY.m 
%  
%       for a full test suite.
%
% This function was written on 2024_08_05 by Sean Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history
% 2024_08_05 - Sean Brennan
% -- Created function by copying out of load script in Geometry library
% 2024_08_09 - Jiabao Zhao
% -- Added format string as a optional input


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
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
        narginchk(1,3);

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

% 
% % Does user want to specify station_tolerance?
% test_date_string = '2024_06_28'; % The date of testing. This defines the folder where the data should be found within LargeData main folder
% if (1<=nargin)
%     temp = varargin{1};
%     if ~isempty(temp)
%         test_date_string = temp;
%     end
% end
% 
% % Does user want to specify station_tolerance?
% vehicle_pose_string = 'VehiclePose_ENU.mat'; % The name of the file containing VehiclePose
% if (2<=nargin)
%     temp = varargin{2};
%     if ~isempty(temp)
%         vehicle_pose_string = temp;
%     end
% end
% 
% % Does user want to specify station_tolerance?
% LIDAR_file_string   = 'Velodyne_LiDAR_Scan_ENU.mat'; % The name of the file containing the LIDAR data
% if (3<=nargin)
%     temp = varargin{3};
%     if ~isempty(temp)
%         LIDAR_file_string = temp;
%     end
% end
% 
% % Does user want to specify flag_load_all_data?
% flag_load_all_data = 0; % FORCE LOAD? Set this manually to 1 to FORCE load
% if (4<=nargin)
%     temp = varargin{4};
%     if ~isempty(temp)
%         flag_load_all_data = temp;
%     end
% end

% Set plotting defaults
plot_str = 'k';
plot_type = 1;  % Plot type refers to 1: a string is given or 2: a color is given - default is 1

% Check to see if user passed in a string or color style?
if 2 <= nargin
    input = varargin{1};
    if ~isempty(input)
        plot_str = input;
        if isnumeric(plot_str)  % Numbers are a color style
            plot_type = 2;
        end
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (2<=nargin)
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


% Nothing to do here - just plotting

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

    % Plot start and end points
    plot(VehiclePose(1,1),VehiclePose(1,2),'.','Color',[0 1 0],'MarkerSize',20);
    plot(VehiclePose(end,1),VehiclePose(end,2),'o','Color',[1 0 0],'MarkerSize',10);


    title('ENU XY plot of vehicle trajectory');
    xlabel('East position [m]');
    ylabel('North position [m]');
    zlabel('Up position [m]');

    % make plots
    if plot_type==1
        if length(plot_str)>3
            eval_string = sprintf('plot(VehiclePose(:,1),VehiclePose(:,2),%s)',plot_str);
            eval(eval_string);
        else
            plot(VehiclePose(:,1),VehiclePose(:,2),plot_str);
        end
    elseif plot_type==2
        plot(VehiclePose(:,1),VehiclePose(:,2),'Color',plot_str);
    end

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



