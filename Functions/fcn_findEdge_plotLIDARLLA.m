function fcn_findEdge_plotLIDARLLA(LIDAR_ENU,varargin)
%% fcn_findEdge_plotLIDARLLA      
% plot the LIDAR data in LLA coordinates
% 
% FORMAT: 
%
%       fcn_findEdge_plotLIDARLLA(LIDAR_ENU,(LIDAR_intensity),(scaling),(color_map),(marker_size),(reference_LLA),(fig_num))
%
% INPUTS:
%
%       LIDAR_ENU: ENU data of the LIDAR points which will be converted to
%       LLA coordinates.
%
%       (OPTIONAL INPUTS)
%
%       (LIDAR_intensity): The intensity of the LIDAR data during mapping as a
%       [Nx1] vector. Higher intensity value means the target surface is
%       more reflective.
%
%       (scaling): scales the intensity of the points
%
%       (color_map): allows for user input color maps, default is 'sky'
%
%       (marker_size): size of the markers, default is 5
%  
%       (reference_LLA):the [reference latitude reference_longitude
%       reference_altitude] of the mapping vehicle during mapping. This
%       have to be in LLA coordinates.
%
%       (format): A format string, e.g. 'b-', that dictates the plot style or
%       a color vector, e.g. [1 0 0.23], that dictates the line color.
%
%       (fig_num): a figure number to plot results. If set to -1, skips any
%       input checking or debugging, no figures will be generated, and sets
%       up code to maximize speed.
%
% OUTPUTS: 
%       (none, just the plot)
% 
% DEPENDENCIES:
%
% EXAMPLES:
%   
%   See the script:
%   
%   script_test_fcn_findEdge_plotLIDARLLA
%
%   for a full test suite
%
% This function was written on 2024_08_07 by Aleksandr Goncharov
% Questions or comments? -- Contact opg5041@psu.edu or 267-304-8354
%
% Revision History:
% 2024_08_06 -Aleksandr Goncharov
% -- Functionalized code from the FindEdge Demo 
% 2024_08_07 -Jiabao Zhao
% -- reordered and simplified the inputs, allowing variable input arguments
% -- minor clean-up of comments.
% 2024_08_09 -Jiabao Zhao
% -- Added format string as optional input. String could be marker size,
% shape or color.

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==8 && isequal(varargin{end},-1))
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
%
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
        narginchk(1,8);
    end
end 

%Does user want to enter LIDAR_intensity?
LIDAR_intensity = [];
flag_simplePlot = 1;
if (2<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        LIDAR_intensity = temp;
        flag_simplePlot = 0;
    end
end

%Does user want to specify scaling?
scaling=3;
if (3<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        scaling = temp;
    end
end

%Does user want to specify color_map?
color_map='jet';
if (4<=nargin)
    temp = varargin{3};
    if ~isempty(temp)
        color_map = temp;
    end
end

%Does user want to specify marker_size?
marker_size=5;
if (5<=nargin)
    temp = varargin{4};
    if ~isempty(temp)
        marker_size = temp;
    end
end

%Does user want to specify reference_LLA?
reference_LLA = [];
if (6<=nargin)
    temp = varargin{5};
    if ~isempty(temp)
        reference_LLA = temp;
    end
end

%Does user want to specify format?
plot_str1 = 'k.';
plot_str2 = 'mo';
plot_type = 1;  % Plot type refers to 1: a string is given or 2: a color is given - default is 1

% Check to see if user passed in a string or color style?
if 7 <= nargin
    input = varargin{6};
    if ~isempty(input)
        plot_str1 = input;
        if isnumeric(plot_str1)  % Numbers are a color style
            plot_type = 2;
        end
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (7<=nargin)
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

% Use the class to convert LLA to ENU
if ~isempty(reference_LLA)
    gps_object = GPS(reference_LLA(1), reference_LLA(2), reference_LLA(3)); % Load the GPS class
else
    gps_object = GPS();
end
concatenate_LiDAR_LLA_points = gps_object.ENU2WGSLLA(LIDAR_ENU(:,1:3));

if 0==flag_simplePlot
    intensity_max=max(LIDAR_intensity);
    intensity_min=min(LIDAR_intensity);
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
    figure(fig_num)
    clf

    if 1==flag_simplePlot
        if plot_type==1
            if length(plot_str1)>3
                eval_string = sprintf('concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),%s)',plot_str1);
                eval(eval_string);
                hold on
                eval_string = sprintf('concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),%s)',plot_str2);
                eval(eval_string);
            else
                geoplot(concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),plot_str1);
            end
        elseif plot_type==2
            geoplot(concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),'Color',plot_str1);
        end
        % % Plot the LIDAR data simply as magenta and black points
        % geoplot(concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),'mo','MarkerSize',10);
        % hold on
        % geoplot(concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),'k.','MarkerSize',10);

    else
        intensity_fraction =  scaling*(LIDAR_intensity - intensity_min)/(intensity_max - intensity_min);

        % Fix intensity fraction to be between 0 and 1
        intensity_fraction = min(max(intensity_fraction,0),1);
        % Use user-defined colormap_string to map intensity to colors. For a
        % full example, see fcn_geometry_fillColorFromNumberOrName
        old_colormap = colormap;
        % color_ordering = colormap('hot');
        color_ordering = flipud(colormap(color_map));
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
            geoplot(concatenate_LiDAR_LLA_points(index_in_this_color,1),concatenate_LiDAR_LLA_points(index_in_this_color,2), '.','Color',color_vector,'MarkerSize',marker_size);
        end
    end
    geobasemap satellite
    geotickformat -dd  % Sets the tick marks to decimal format, not degrees/minutes/seconds which is default
    if flag_do_debug
        fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
    end

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