function [Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_findEdge_findMaxMinOfXYZ(N_points,varargin)
%% fcn_findEdge_findMaxMinOfXYZ
% Find the maximum and minimum numbers of x, y, and z of the given points 
% 
% FORMAT:
%
%      [Min_x,Max_x,Min_y,Max_y,Min_z,Max_z] = fcn_findEdge_findMaxMinOfXYZ(N_points,(fig_num))
%
% INPUTS:     
%       
%      N_points: Number of points (in x,y,z coordinates)
%      
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      Min_x: minimum number of x-coordinate of points
%
%      Max_x  maximum number of x-coordinate of points
%
%      Min_y: minimum number of y-coordinate of points
%
%      Max_y: maximum number of y-coordinate of points
%
%      Min_z: minimum number of z-coordinate of points
%
%      Max_z: maximum number of z-coordinate of points
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_findEdge_findMaxMinOfXYZ.m for a full
%       test suite.
%
% This function was written on 2024_08_05 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu

% Revision history
% 2024_08_05 - Jiabao Zhao
% -- Created function by copying out of Geometry library


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
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
        narginchk(1,2);

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

% Does user want to specify fig_num?
flag_do_plots = 0;
if 2<= nargin && 0==flag_max_speed
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
if isscalar(N_points(1,:)) %% N by 1 
    Min_x = min(N_points(:,1));
    Max_x = max(N_points(:,1));
    Min_y = [];
    Max_y = [];
    Min_z = [];
    Max_z = [];
elseif 2 == length(N_points(1,:)) % N by 2
    Min_x = min(N_points(:,1));
    Max_x = max(N_points(:,1));
    Min_y = min(N_points(:,2));
    Max_y = max(N_points(:,2));
    Min_z = [];
    Max_z = [];
elseif 3 == length(N_points(1,:)) % N by 3 
    Min_x = min(N_points(:,1));
    Max_x = max(N_points(:,1));
    Min_y = min(N_points(:,2));
    Max_y = max(N_points(:,2));
    Min_z = min(N_points(:,3));
    Max_z = max(N_points(:,3));
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
    % axis equal
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    title('ENU Trace Geometry')

    %plot of points
    if isscalar(N_points(1,:))  % N by 1
        [~, idx_min_x] = min(N_points(:,1));
        [~, idx_max_x] = max(N_points(:,1));
        x = (N_points(:,1));
        hold on
        plot(x,'o');
        plot(idx_min_x,x(idx_min_x),'ro'); % Min x
        plot(idx_max_x,x(idx_max_x),'ro'); % max x
        hold off
    elseif 2 == length(N_points(1,:)) % N by 2
        [~, idx_min_x] = min(N_points(:,1));
        [~, idx_max_x] = max(N_points(:,1));
        [~, idx_min_y] = min(N_points(:,2));
        [~, idx_max_y] = max(N_points(:,2));
        x = (N_points(:,1));
        y = (N_points(:,2));
        hold on
        scatter(x,y);
        scatter(x(idx_min_x),y(idx_min_x),100,'r'); % Min x
        scatter(x(idx_max_x),y(idx_max_x),100,'r'); % max x
        scatter(x(idx_min_y),y(idx_min_y),200,'g'); % min y
        scatter(x(idx_max_y),y(idx_max_y),200,'g'); % max y
        hold off
    elseif 3 == length(N_points(1,:))  % N by 3 
        [~, idx_min_x] = min(N_points(:,1));
        [~, idx_max_x] = max(N_points(:,1));
        [~, idx_min_y] = min(N_points(:,2));
        [~, idx_max_y] = max(N_points(:,2));
        [~, idx_min_z] = min(N_points(:,3));
        [~, idx_max_z] = max(N_points(:,3));
        x = (N_points(:,1));
        y = (N_points(:,2));
        z = (N_points(:,3));
        hold on 
        colors = z;
        scatter3(x,y,z,20,colors);
        % scatter3(x,y,z,20)
        scatter3(x(idx_min_x),y(idx_min_x),z(idx_min_x),100, 'r'); % Min x
        scatter3(x(idx_max_x),y(idx_max_x),z(idx_max_x),100, 'r'); % max x
        scatter3(x(idx_min_y),y(idx_min_y),z(idx_min_y),200, 'g'); % min y
        scatter3(x(idx_max_y),y(idx_max_y),z(idx_max_y),200, 'g'); % max y
        scatter3(x(idx_min_z),y(idx_min_z),z(idx_min_z),300, 'y'); % min z
        scatter3(x(idx_max_z),y(idx_max_z),z(idx_max_z),300, 'y'); % max z
        colormap('jet'); % You can use 'viridis' if you have it, or other colormaps
        colorbar; % Add a color bar to show the color mapping
        view(3)
        hold off
    end
    hold off
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



