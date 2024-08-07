function fcn_findEdge_plotLIDARXY(LIDAR_ENU,varargin)
%% fcn_findEdge_plotLIDARXY
% 
% FORMAT: 
% fcn_findEdge_plotLiDARXY(LIDAR_ENU,(Color Triplet),(Marker_Size),(fig_num))
%
% INPUTS:
%       
%       LIDAR_ENU: Dataset that will be plotted
%       
% (OPTIONAL INPUTS):
%
%       (Color Triplet): A color triplet which determines plotted color
%       default is [0 0 1] which is blue.
%
%       (Marker Size): Size of the points. Default is 5.
%
%       (fig_num): a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
% OUTPUTS: 
%   (None)

% DEPENDENCIES:
%   (None)
% EXAMPLES:
%   
%   See the script:
%   
%   script_test_fcn_findEdge_plotLIDARXY
%
%   for a full test suite
%
% This function was written on 2024_08_06 by Aleksandr Goncharov
% Questions or comments? opg5041@psu.edu -  267-304-8354
%
% Revision History
% 2024_08_06 - Aleksandr Goncharov
% -- Created the function by taking parts of the code from the
% script_demo_FindEdge and functionalizing it.


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
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

%Check correct number of inputs.
if flag_max_speed == 0
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(1,4);
    end
end 

%Does user want to specify color_triplet?
color_triplet=[0 0 1];
if (2<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        color_triplet = temp;
    end
end

%Does user want to specify marker size?
marker_size=5;
if (3<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        marker_size = temp;
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (4<=nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(NO MAIN CODE, JUST PLOTS)

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
    figure(fig_num);
    clf;

    % Plot the LIDAR in XZ ENU
    subplot(1,2,1)

    hold on;
    grid on;
    axis equal

    xlabel('East position [m]');
    ylabel('Up position [m]');

    % Plot the LIDAR data
    plot(LIDAR_ENU(:,1),LIDAR_ENU(:,3), '.','Color',color_triplet,'MarkerSize',marker_size);


    % Plot the LIDAR in YZ ENU
    subplot(1,2,2)
    hold on;
    grid on;
    axis equal

    xlabel('North position [m]');
    ylabel('Up position [m]');

    % Plot the LIDAR data
    plot(LIDAR_ENU(:,2),LIDAR_ENU(:,3), '.','Color',color_triplet,'MarkerSize',marker_size);

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
end