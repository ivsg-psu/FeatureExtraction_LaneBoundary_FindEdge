function fcn_findEdge_plotLIDARXYZ(LIDAR_intensity,LIDAR_ENU,varargin)
%% fcn_findEdge_plotLIDARXYZ
% Plots the LIDAR points on an XYZ scale
%
% FORMAT: 
% fcn_findEdge_plotLIDARXYZ(LIDAR_intensity,LIDAR_ENU,(simple_flag),(scaling),(color_map),(fig_num))
%
% INPUTS:
%       
%       LIDAR_intensity: LiDAR intensity scalar
%
%       LIDAR_ENU: LiDAR data points corresponding to the intensity
%      
%       (OPTIONAL INPUTS)
%
%       simple_flag: a toggle flag which allows for the most basic plotting
%       with no color input. Default is off. 
%
%       scaling: Scaling for the intensity. Default is a factor of 3.
%
%       color_map: Color map for the plot, default is "hot".
%
%       fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
%
% OUTPUTS: 
%   (None)
% DEPENDENCIES:
%   (none)
%
% EXAMPLES:
%   
%   See the script:
%   
%   script_test_fcn_findEdge_plotLIDARXYZ.m
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
if (nargin==6 && isequal(varargin{end},-1))
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

if flag_max_speed == 0
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(2,6);
    end
end 

%Does user want to toggle simple_flag?
simple_flag=0;
if (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        simple_flag = temp;
    end
end

%Does user want to specify scaling?
scaling=3;
if (4<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        scaling = temp;
    end
end

%Does user want to specify color_map?
color_map="jet";
if (5<=nargin)
    temp = varargin{3};
    if ~isempty(temp)
        color_map = temp;
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (6<=nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

%% calculate the intensity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the intensity ranges
intensity_min = min(LIDAR_intensity);
intensity_max = max(LIDAR_intensity);



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

    if 1==simple_flag
        % Plot only the LIDAR data simply as blue points
        plot3(LIDAR_ENU(:,1),LIDAR_ENU(:,2),LIDAR_ENU(:,3), '.','Color',[0 0 1],'MarkerSize',1);
    else
        intensity_fraction = scaling*LIDAR_intensity/(intensity_max - intensity_min);

        % Use user-defined colormap_string to map intensity to colors. For a
        % full example, see fcn_geometry_fillColorFromNumberOrName
        old_colormap = colormap;
        color_ordering = colormap(color_map);
        colormap(old_colormap);
        N_colors = length(color_ordering(:,1));

        % Make sure the plot number is a fraction between 0 and 1
        plot_number = min(max(0,intensity_fraction),1);

        % Convert the plot number to a row
        color_row = floor((N_colors-1)*plot_number) + 1;


        % Plot the LIDAR data with intensity
        for ith_color = min(color_row):max(color_row)
            % Find the color
            color_vector = color_ordering(ith_color,:);

            % Find all the points that are in this color
            index_in_this_color = find(color_row==ith_color);
            plot3(...
                LIDAR_ENU(index_in_this_color,1),...
                LIDAR_ENU(index_in_this_color,2),...
                LIDAR_ENU(index_in_this_color,3), '.','Color',color_vector,'MarkerSize',5);
        end
    end
end



end


