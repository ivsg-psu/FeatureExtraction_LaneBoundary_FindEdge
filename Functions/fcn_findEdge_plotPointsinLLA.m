function LLA_data = fcn_findEdge_plotPointsinLLA(ENU_data,marker_size,RGB_triplet,varargin)
% This function takes in various boundary points and then generates a plot
% that shows where the boundary,driveable/non-driveable area, and unmapped area
% 
% FORMAT:  
% [LLA_data] =
%       fcn_findEdge_plotPointsinLLA(ENU_data,marker_size,RGB_triplet,(marker_type),(legend_options),...
%       (legend_name),(legend_position),(reference_latitude),(reference_longitude),(reference_altitude),(fig_num))
%
% INPUTS:
% ENU_data: Data points 
% marker_size: size of the points that are going to be plotted
% RGB_triplet: Color of the markers
% (OPTIONAL INPUT)
% marker_type: change the marker type, default is "."
% legend_options: enable the plot to display a legend
% legend_name: name of the legend
% legend_position: position of the legend
% reference_latitude
% reference_longitude
% reference_altitude
% fig_num: figure number
%
% OUTPUTS:
% LLA coordinates of each point
%
% DEPENDENCIES:
% GPS CLASS
%
% EXAMPLES:
% See script: script_test_fcn_findEdge_plotPointsinLLA
%
% Revision History
% 2024_07_10 - Aneesh Batchu 
% -- wrote the code originally
% 2024_07_10 - Aleksandr Goncharov
% -- Funtionlized this code
% 2024_07_11 - Aleksandr Goncharov
% -- Shortened and cleaned up the function to be more universal
% 2024_07-12 - Aneesh Batchu
% -- Changed the name of the function to "fcn_geometry_plotGridCenters"
% -- Made minor changes in the legend position options
% 2024_07_17 - Aleksandr Goncharov
% -- Added max speed flag
% 2024_07_29 - Aleksandr Goncharov
% -- Changed name to "fcn_geometry_plotPointsinLLA"
% -- Updated function to include more inputs such as marker size/color
% 2024_08_07 - Jiabao Zhao
% -- pull this code from geometry class and rename it

%%
flag_max_speed = 0;
if (nargin==11 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS");
    MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG = getenv("MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG);
    end
end
if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838; %#ok<NASGU>
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

if flag_max_speed==0
    % Are there the right number of inputs?
    narginchk(3,11);
end

%does user want a different marker type
marker_type=".";
if 4<=nargin
    temp=varargin{1};
    if ~isempty(temp)
        marker_type=temp;
    end
end

%does user want a legend
flag_create_legend=0;
if 5<=nargin
    temp=varargin{2};
    if temp==1
        flag_create_legend=1;
    end
end

%legend name
if 6<=nargin
    temp=varargin{3};
    if ~isempty(temp)
        legend_name=temp;
    end
end

%legend position
legend_position = [];
if 7<=nargin
    temp=varargin{4};
    if ~isempty(temp)
        legend_position=temp;
    end
end

%reference_latitude
reference_latitude = 40.86368573;
if 8<=nargin
    temp=varargin{5};
    if ~isempty(temp)
        reference_latitude=temp;
    end
end

%reference_longitude
reference_longitude = -77.83592832;
if 9<=nargin
    temp=varargin{6};
    if ~isempty(temp)
        reference_longitude=temp;
    end
end

%reference_altitude
reference_altitude = 344.189;
if 10<=nargin
    temp=varargin{7};
    if ~isempty(temp)
        reference_altitude=temp;
    end
end


% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (11<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end


%% Write main code for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % True boundary points
% Trace_coordinates = [true_boundary_points,zeros(length(true_boundary_points),1)]; 

 % Define GPS object
gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Load the GPS class

% Use the class to convert LLA to ENU
LLA_data = gps_object.ENU2WGSLLA(ENU_data);


%% Any debugging?
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

if flag_do_plots && flag_max_speed==0

% Plot the LLA boundary points
figure(fig_num);

% Plot
geoplot(LLA_data(:,1),LLA_data(:,2),marker_type,'MarkerSize',marker_size,'Color',RGB_triplet,'DisplayName',legend_name);
hold on

%Adding Legend
if flag_create_legend == 1
    if ~isempty(legend_position)
        legend('Position',legend_position);
    else
        legend()
    end

end

title('Grid centers in LLA ')
geobasemap satellite
geotickformat -dd

end
end
