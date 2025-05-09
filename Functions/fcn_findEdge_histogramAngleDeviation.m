function theta_threshold = fcn_findEdge_histogramAngleDeviation(angle_btw_unit_normals_and_vertical, angle_btw_unit_normals_and_vertical_driven_path,mean_angle_btw_unit_normals_and_vertical_driven_path, chosen_theta_threshold, varargin)
%% fcn_findEdge_histogramAngleDeviation   Histogram of angle deviation
% 
% FORMAT:
%
%      fcn_findEdge_histogramAngleDeviation(angle_btw_unit_normals_and_vertical, ...
%      angle_btw_unit_normals_and_vertical_driven_path,(fig_num))
%
% INPUTS:     
%       
%      angle_btw_unit_normals_and_vertical: angle between normal and
%      vertical direction.  
%
%      angle_btw_unit_normals_and_vertical_driven_path: angle between normal and
%      vertical driven path.  
%
%      mean_angle_btw_unit_normals_and_vertical_driven_path:  mean angle 
%      between normal and vertical driven path.  
%      
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%  
%       (no output, just plot)
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_findEdge_histogramAngleDeviation for a full
%       test suite.
%
% This function was written on 2024_08_13 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu

% Revision history: 
% 2024_07_20 - Aneesh Batchu
% -- Wrote the code originally
% 2024_08_13 - Jiabao Zhao
% -- Functionalized code from the FindEdge Demo 
% 2024_08_19 - Aneesh Batchu
% -- Replaced std_threshold (input) to chosen_std_threshold
% -- std_threshold is calculated here when chosen_std_threshold is empty


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
        narginchk(4,5);
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
% just plot
% theta_threshold = 0.1745;

theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 3*std(angle_btw_unit_normals_and_vertical_driven_path);

if ~isempty(chosen_theta_threshold)
    theta_threshold = chosen_theta_threshold;
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
    figure(fig_num);clf;
    hold on
    grid on
    xlabel('Angle deviation in z error after plane fit')
    ylabel('Frequency')
    title('Determining suitable angle deviation')
    % histogram(standard_deviation_in_z)
    histogram(angle_btw_unit_normals_and_vertical,100)
    histogram(angle_btw_unit_normals_and_vertical_driven_path,5)

    x_coord = mean_angle_btw_unit_normals_and_vertical_driven_path + 20*std(angle_btw_unit_normals_and_vertical_driven_path);
    % theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 3.3*std(angle_btw_unit_normals_and_vertical_driven_path);
    % theta_threshold = mean_angle_btw_unit_normals_and_vertical_driven_path + 3*std(angle_btw_unit_normals_and_vertical_driven_path);

    % plot(theta_threshold,max(angle_btw_unit_normals_and_vertical), 'k.', 'MarkerSize',20)
    % theta_threshold = 0.1504;
    % theta_threshold = 0.3;

    disp('mean of angle_deviation_driven_path')
    disp(mean_angle_btw_unit_normals_and_vertical_driven_path*(180/pi))
    disp('theta threshold')
    disp(theta_threshold*(180/pi))

    plot(theta_threshold,0, 'k.', 'MarkerSize',18)
    current_text = sprintf('theta threshold = %.1f',theta_threshold*(180/pi));
    text(x_coord, 30,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');
end
end