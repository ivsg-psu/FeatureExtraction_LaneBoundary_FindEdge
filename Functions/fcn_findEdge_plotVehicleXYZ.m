function fcn_findEdge_plotVehicleXYZ(VehiclePose, varargin)
%% fcn_findEdge_plotVehicleXYZ
% Plots the vehicle's XYZ position.
% 
% FORMAT:
%
%      fcn_findEdge_loadLIDARData(VehiclePose, (scanLineRange), (fig_num))
%
% INPUTS:  
%
%      VehiclePose: the position of the mapping vehicle during mapping
%      
%      (OPTIONAL INPUTS)
%       
%      scanLineRange: the [min_number max_number] of the scanlines to
%      extract. The default is to use all scan lines.
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
%       script_test_fcn_findEdge_plotVehicleXYZ.m 
%  
%       for a full test suite.
%
% This function was written on 2024_08_05 by Sean Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history
% 2024_08_05 - Sean Brennan
% -- Created function by copying out of load script in Geometry library


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

% Does user want to specify scanLineRange?
scanLineRange = []; 
if (2<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        scanLineRange = temp;
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (3<=nargin)
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

N_scanLines = length(VehiclePose(:,1));

%% Calculate the vehicle orientation
vehicle_change_in_pose_XY = diff(VehiclePose(:,1:2));

% Repeat the last value again, since diff removes one row. We want the same
% number of vectors as the number of points, and diff removed one point.
vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];

% Convert these to unit vectors
unit_vehicle_change_in_pose_XY = fcn_INTERNAL_calcUnitVector(vehicle_change_in_pose_XY);

% Find orthogonal vetors by rotating by 90 degrees in the CCW direction by
% multiplying by the rotatio matrix
unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY*[0 1; -1 0];


%% Get the scanLineRange
if ~isempty(scanLineRange)
    scanLineNumber_start = scanLineRange(1);
    scanLineNumber_end   = scanLineRange(2);
    if scanLineNumber_start<1
        warning('on','backtrace');
        warning('Incorrect scan line detected');
        error('Starting scan line must be greater than or equal to 1');
    end
    if scanLineNumber_start>N_scanLines
        warning('on','backtrace');
        warning('Incorrect scan line detected');
        error('Starting scan line must be less than or equal to number of scans: %.0f',N_scanLines);
    end

else
    scanLineNumber_start = 1;
    scanLineNumber_end   = N_scanLines;
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

    hold on;
    grid on;
    axis equal

    if 1==0
        % Plot the vehicle's trajectory
        plot3(VehiclePose(:,1),VehiclePose(:,2),VehiclePose(:,3),'-','Color',[1 0 1],'MarkerSize',10,'LineWidth',3);

        % Plot start and end points of trajectory
        plot3(VehiclePose(1,1),VehiclePose(1,2),VehiclePose(1,3),'.','Color',[1 0 0],'MarkerSize',10,'LineWidth',3);
        plot3(VehiclePose(end,1),VehiclePose(end,2),VehiclePose(end,3),'o','Color',[0 1 0],'MarkerSize',10,'LineWidth',3);

    else
        plot3(VehiclePose(scanLineNumber_start:scanLineNumber_end,1),VehiclePose(scanLineNumber_start:scanLineNumber_end,2),VehiclePose(scanLineNumber_start:scanLineNumber_end,3),'.','Color',[0 0 0],'MarkerSize',30,'LineWidth',3);
        plot3(VehiclePose(scanLineNumber_start:scanLineNumber_end,1),VehiclePose(scanLineNumber_start:scanLineNumber_end,2),VehiclePose(scanLineNumber_start:scanLineNumber_end,3),'.','Color',[1 1 0],'MarkerSize',10,'LineWidth',3);

        % Show the orthogonal arrows showing vehicle motion directions. Green
        % is forward, bLue is Left
        quiver3(...
            VehiclePose(scanLineNumber_start,1),VehiclePose(scanLineNumber_start,2),VehiclePose(scanLineNumber_start,3), ...
            unit_vehicle_change_in_pose_XY(scanLineNumber_start,1),unit_vehicle_change_in_pose_XY(scanLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 1 0]);
        quiver3(...
            VehiclePose(scanLineNumber_start,1),VehiclePose(scanLineNumber_start,2),VehiclePose(scanLineNumber_start,3), ...
            unit_ortho_vehicle_vectors_XY(scanLineNumber_start,1),unit_ortho_vehicle_vectors_XY(scanLineNumber_start,2),0,0,'-','LineWidth',3,'Color',[0 0 1]);
    end

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

    hold off
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


%% fcn_INTERNAL_calcUnitVector
function unit_vectors = fcn_INTERNAL_calcUnitVector(input_vectors)
vector_length = sum(input_vectors.^2,2).^0.5;
unit_vectors = input_vectors./vector_length;
end
