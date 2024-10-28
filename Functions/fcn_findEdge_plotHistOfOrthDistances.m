function transverse_distances = fcn_findEdge_plotHistOfOrthDistances(reference_path, computed_boundary_path, fig_num1, fig_num2, varargin)
%% fcn_findEdge_plotHistOfOrthDistances
% The histogram of the orthogonal distances from the reference points to
% computed boundary points is plotted. 
%
% FORMAT: 
%
% transverse_distances = fcn_findEdge_plotHistOfOrthDistances(reference_path, computed_boundary_path, fig_num1, fig_num2, (fig_num))
%
% INPUTS:
%
%      reference_path: The transverse distances of the computed boundary
%      path is computed from this "reference path". The reference path must
%      be arranged in counter clockwise direction
% 
%      computed_boundary_path: The path of the computed boundary points 
%
%      fig_num1: Figure number of histogram
%
%      fig_num2: Figure number of LLA plot
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: figure number, can be set to -1 for fast mode.
%
% OUTPUTS:
%
%      transverse_distances: the orthogonal distances from the reference
%      path to computed boundary path
%
% DEPENDENCIES:
%
%      fcn_Path_convertXY2St
%
% EXAMPLES:
%      
% See the script: script_test_fcn_findEdge_plotHistOfOrthDistances
% for a full test suite.
%
% This function was written on 2024_10_24 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu or sbrennan@psu.edu 

% Revision history:
% 2024_10_24 - Aneesh Batchu
% - wrote the code originally


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS");
    MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG = getenv("MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_FINDEDGE_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_FINDEDGE_FLAG_CHECK_INPUTS);
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

if 0==flag_max_speed
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(4,5);
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (5<= nargin)
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

St_points = fcn_Path_convertXY2St(reference_path,computed_boundary_path,[],[]);
transverse_distances = real(St_points(:,2));

max_transverse_error = 2; % Set by user - this determines the color scale
normalized_abs_transverse_error = min(2,abs(transverse_distances))/max_transverse_error;


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

    %% Hist Plot
    figure(fig_num1);

    histogram(transverse_distances,30);

    temp = axis;
    %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];

    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.1;
    axis([temp(1), temp(2),  temp(3), temp(4)+percent_larger*axis_range_y]);

    sigma_std = std(transverse_distances);
    mu_mean = mean(transverse_distances);

    text(mu_mean - 3.2*sigma_std, 50, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
    text(mu_mean - 3.2*sigma_std, 40, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');

    %% LLA Plot

    XYI_data = [computed_boundary_path normalized_abs_transverse_error];

    % Call the plotting function
    clear plotFormat
    plotFormat.LineStyle = 'none';
    plotFormat.LineWidth = 5;
    plotFormat.Marker = '.';
    plotFormat.MarkerSize = 10;

    colorMapMatrixOrString = colormap('turbo');
    Ncolors = 10;
    reducedColorMap = fcn_plotRoad_reduceColorMap(colorMapMatrixOrString, Ncolors, -1);

    % Convert ENU data to LLA
    % Define reference location at test track
    reference_latitude = 40.86368573;
    reference_longitude = -77.83592832;
    reference_altitude = 344.189;

    % Location for Test Track base station
    setenv('MATLABFLAG_PLOTROAD_REFERENCE_LATITUDE','40.86368573');
    setenv('MATLABFLAG_PLOTROAD_REFERENCE_LONGITUDE','-77.83592832');
    setenv('MATLABFLAG_PLOTROAD_REFERENCE_ALTITUDE','344.189');

    gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Initiate the class object for GPS
    % ENU_full_data = gps_object.WGSLLA2ENU(LLA_full_data(:,1), LLA_full_data(:,2), LLA_full_data(:,3));
    ENU_full_data = [XYI_data(:,1), XYI_data(:,2), XYI_data(:,2)*0];
    LLA_full_data  =  gps_object.ENU2WGSLLA(ENU_full_data);

    figure(fig_num2);
    clf;
    LLI_data = [LLA_full_data(:,1:2) normalized_abs_transverse_error];

    [~, ~]  = fcn_plotRoad_plotLLI(LLI_data, (plotFormat),  (reducedColorMap), (fig_num2));
    title('Boundary error, ENU');

    % Force the plot to fit
    geolimits('auto');
    %% ENU plot

    figure(fig_num);
   
    fcn_plotRoad_plotXYI(XYI_data, (plotFormat),  (reducedColorMap), (fig_num));
    title('Boundary error, ENU');
end

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
