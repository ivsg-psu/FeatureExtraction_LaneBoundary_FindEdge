function LLA_VehiclePose = fcn_findEdge_plotVehicleLLA(VehiclePose, varargin)
%fcn_findEdge_plotVehicleLLA    plot the trace of the mapping vehicle during the mapping
% 
% FORMAT:
%
%      LLA_VehiclePose = fcn_findEdge_plotVehicleLLA(VehiclePose, (reference_LLA), (zoom_in_location), (fig_num))
%
% INPUTS:  
%
%      VehiclePose: the position of the mapping vehicle during mapping in
%      ENU coordinates.
%       
%      (OPTIONAL INPUTS)
%
%      reference_LLA: the [reference latitude reference_longitude
%      reference_altitude] of the mapping vehicle during mapping. This
%      have to be in LLA coordinates.
%
%      zoom_in_location: zoom into the location we want. This variable
%      usually contains latitude and longitude. 
%      
%      zoomLevel: The elevation level we want to set for the figure to
%      view.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.     
%
% OUTPUTS:
%
%      LLA_VehiclePose = the position of the mapping vehicle during mapping
%      in LLA coordinates.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
% 
%       script_test_fcn_findEdge_plotVehicleLLA
%  
%       for a full test suite.
%
% This function was written on 2024_08_06 by Jiabao Zhao using code from
% Aneesh Batchu. Edits by S. Brennan.
% Questions or comments? jpz5469@psu.edu or sbrennan@psu.edu

% Revision history
% 2024_08_06 - Jiabao Zhao
% -- Functionalize this code from the "script_demo_Find_Edge". 
% 2024_08_07 - S. Brennan
% -- reordered and simplified the inputs, allowing variable input arguments
% -- minor clean-up of comments. Function overall looks good.
% 2024_08_07 - Jiabao Zhao
% -- Add "zoomLevel" as one of the input arguments.


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
        narginchk(1,5);

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

% Does user want to specify reference_LLA?
reference_LLA = [];
if (2<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        reference_LLA = temp;
    end
end

% Does user want to specify zoom_in_location?
zoom_in_location = [];
if (3<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        zoom_in_location = temp;
    end
end

zoomLevel = [];
if (4<=nargin)
    temp = varargin{3};
    if ~isempty(temp)
        zoomLevel = temp;
    end
end


% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (4<=nargin)
    temp = varargin{end};
    if ~isempty(temp)
        LLA_fig_num = temp;
        flag_do_plots = 1;
    end
end

%% Solve for the gps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(reference_LLA)
    gps_object = GPS(reference_LLA(1), reference_LLA(2), reference_LLA(3)); % Load the GPS class
else
    gps_object = GPS();
end

% Find the vehicle's LLA position using ENU inputs
LLA_VehiclePose = gps_object.ENU2WGSLLA(VehiclePose(:,1:3)); % convert data position from ENU to LLA coordination

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
% Do plot, and set it up so that zoom is correct, tick marks are correct,
% map is centered where we want, etc. To see options, do the geoplot, then
% zoom into the location we want, and then do:
%
% format long % This sets the format so we can see all the digits
% temp = gca  % This Gets Current Axis (GCA) and saves it as temp
%
% Next, click on "show all properties" and copy out the ZoomLevel and
% MapCenter values into the correct locations below. This way, every time
% we run this script, it automatically zooms and centers onto the correct
% location.
if flag_do_plots
    figure(LLA_fig_num)
    clf

    h_geoplot = geoplot(LLA_VehiclePose(:,1),LLA_VehiclePose(:,2),'-','Color',[0 0 1],'MarkerSize',10);
    hold on;
    h_parent =  get(h_geoplot,'Parent');

    if ~isempty(zoomLevel)
        set(h_parent,'ZoomLevel',zoomLevel);
    end

    if ~isempty(zoom_in_location)
        set(h_parent, 'MapCenter',zoom_in_location);
    end

    geobasemap satellite
    geotickformat -dd  % Sets the tick marks to decimal format, not degrees/minutes/seconds which is default


    % NOTE: transition the geoplotting to use the PlotTestTrack library
    % Plot start and end points
    geoplot(LLA_VehiclePose(1,1),LLA_VehiclePose(1,2),'.','Color',[0 1 0],'MarkerSize',10);
    geoplot(LLA_VehiclePose(end,1),LLA_VehiclePose(end,2),'o','Color',[1 0 0],'MarkerSize',10);

end % Ends check if plotting

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end
end
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
