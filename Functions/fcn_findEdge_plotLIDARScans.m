function fcn_findEdge_plotLIDARScans(LIDAR_ENU, LIDAR_intensity, varargin)
%% fcn_findEdge_plotLIDARScans   
% plot the LIDAR data in LLA coordinates
%    
% NOTE: colormap will be used when LIDAR intensity data is available,
% otherwise debug section will use format as the color of plot.
% FORMAT: 
%
%       fcn_findEdge_plotLIDARLLA_Aneesh(LIDAR_ENU, LIDAR_intensity ,(scaling),(color_map),(marker_size),(reference_LLA),(fig_num))
%
% INPUTS:
%
%       LIDAR_ENU: ENU data of the LIDAR points which will be converted to
%       LLA coordinates.
%
%       LIDAR_intensity：The intensity of LIDAR increases with the strength
%       of the reflection; a greater value indicates a stronger reflection.
%
%       (OPTIONAL INPUTS)
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
%   script_test_fcn_findEdge_plotLIDARLLA_Aneesh
%
%   for a full test suite
%
% This function was written on 2024_08_12 by Aneesh Batchu
% Questions or comments? 
%
% Revision History:

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==7 && isequal(varargin{end},-1))
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
        narginchk(2,7);
    end
end 

%Does user want to specify scaling?
scaling=3;
if (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        scaling = temp;
    end
end

%Does user want to specify color_map?
color_map='sky';
if (4<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        color_map = temp;
    end
end

%Does user want to specify marker_size?
marker_size=5;
if (5<=nargin)
    temp = varargin{3};
    if ~isempty(temp)
        marker_size = temp;
    end
end

%Does user want to specify reference_LLA?
reference_LLA = [];
if (6<=nargin)
    temp = varargin{4};
    if ~isempty(temp)
        reference_LLA = temp;
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
    if 1==0
        % Plot the LIDAR data simply as magenta and black points
        geoplot(concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),'mo','MarkerSize',10);
        geoplot(concatenate_LiDAR_LLA_points(:,1),concatenate_LiDAR_LLA_points(:,2),'k.','MarkerSize',10);
    else
        intensity_min = min(LIDAR_intensity);
        intensity_max = max(LIDAR_intensity);
        intensity_fraction = scaling*LIDAR_intensity/(intensity_max - intensity_min);

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

            % geoplot(concatenate_LiDAR_LLA_points_new(index_in_this_color,1),concatenate_LiDAR_LLA_points_new(index_in_this_color,2), '.','Color',color_vector,'MarkerSize',marker_size);
        end
        geobasemap satellite
        geotickformat -dd  % Sets the tick marks to decimal format, not degrees/minutes/seconds which is default
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