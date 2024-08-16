function [point_density] = fcn_findEdge_determineGridPointDensity(total_points_in_each_grid_with_points_greater_than_zero,total_points_in_each_grid_in_the_driven_path,grid_size,varargin)
%% fcn_findEdge_determineGridPointDensity
% Determine the suitable "point density" for the analysis by comparing the point densities of driven grids with those of neighboring grids
% 
% FORMAT:
%
%[point_density] = fcn_findEdge_determineGridPointDensity(total_points_in_each_grid_with_points_greater_than_zero ...
% ,total_points_in_each_grid_in_the_driven_path,grid_size,(N_bins_grid_with_points_greater_than_zero)...
% (N_bins_grid_in_the_driven_path),(fig_num))
%
% INPUTS:
%      total_points_in_each_grid_with_points_greater_than_zero   
%
%      total_points_in_each_grid_in_the_driven_path
%
%      grid_size
%
% (OPTIONAL INPUTS):
%
%      (N_bins_grid_with_points_greater_than_zero)
% 
%      (N_bins_grid_in_the_driven_path)
%
%      (fig_num): figure number, can be set to -1 for fast mode
%
% OUTPUTS: 
%
%      point_density
%
%
% DEPENDENCIES:
%
% EXAMPLES:
%   See the script:
%   
%   script_test_fcn_findEdge_determineGridPointDensity
%
%   for a full test suite
%
% This function was written by Aneesh Batchu
%
% REVISION HISTORY:
% 2024_07_25 - Aneesh Batchu
% -- wrote the code originally
% 2024_08_13 - Aleksandr Goncharov
% -- Functionalized code from the FindEdge Demo 

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
        narginchk(3,6);
    end
end 

% Does the user want to specify N_bins_grid_with_points_greater_than_zero
N_bins_grid_with_points_greater_than_zero=20;
if (4<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
    N_bins_grid_with_points_greater_than_zero=temp;
    end
end

% Does the user want to specify N_bins_grid_in_the_driven_path
N_bins_grid_in_the_driven_path=10;
if (5<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
    N_bins_grid_in_the_driven_path=temp;
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

[counts1,binEdges] = histcounts(total_points_in_each_grid_with_points_greater_than_zero,N_bins_grid_with_points_greater_than_zero);
[counts2, ~] = histcounts(total_points_in_each_grid_in_the_driven_path,N_bins_grid_in_the_driven_path);

[~,index_max_counts2] = max(counts1);

% point_density = floor(mean(total_points_in_each_grid_in_the_driven_path) - 7.5*(std(total_points_in_each_grid_in_the_driven_path)));
x_location = floor(mean(total_points_in_each_grid_in_the_driven_path) - 1.5*(std(total_points_in_each_grid_in_the_driven_path)));

% Calculate the overlapping counts
% % overlapCounts = min(counts2, counts1);

% Find a ratio
% point_density = sum(binEdges(index_max_counts2:(index_max_counts2+1)))/2; 
% point_density = floor(mean(total_points_in_each_grid_in_the_driven_path) - 7.5*(std(total_points_in_each_grid_in_the_driven_path)));
% point_density = floor(mean(total_points_in_each_grid_in_the_driven_path));

% Minimum number of points required 
point_density = floor(20*((grid_size^2)/(0.3^2)));




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


    % Figure number of histogram
    figure(fig_num);
    % Add labels and title
    hold on
    grid on
    xlabel('Points per grid');
    ylabel('Frequency');
    title('Statistic 1: Determining suitable point density');

    histogram(total_points_in_each_grid_with_points_greater_than_zero,N_bins_grid_with_points_greater_than_zero,'Visible','on');
    histogram(total_points_in_each_grid_in_the_driven_path,N_bins_grid_in_the_driven_path,'Visible','on');

    disp('Chosen point density')
    disp(point_density)

    plot(point_density,0, 'k.', 'MarkerSize',20)
    current_text = sprintf('point density = %.2d',point_density);
    text(x_location, 200,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');

end


end