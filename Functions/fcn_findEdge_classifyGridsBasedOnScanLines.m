function [grid_indices_with_more_than_one_scan_line, gridCenters_with_one_scan_line] = fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points, total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters,varargin)
%% fcn_findEdge_classifyGridsBasedOnScanLines
%Grids with required point density and low point density
%
% FORMAT:
%
%      fcn_findEdge_classifyGridsBasedOnScanLines(grids_greater_than_zero_points...
%      ,total_scan_lines_in_each_grid_with_more_than_zero_points,gridCenters, (format_1),(format_2),
%      (fig_num))
%
% INPUTS:
%      grids_greater_than_zero_points: mask that shows every grid with more than 0 points.      
%
%      total_scan_lines_in_each_grid_with_more_than_zero_points
%
%      gridCenters
%
% (OPTIONAL INPUTS):
%
%      (format_1): string of the format for plot 1
% 
%      (format_2): string of the format for plot 2
%
%      (fig_num): figure number, can be set to -1 for fast mode
%
% OUTPUTS: 
%
%      gridCenters_with_one_scan_line
%
%
% DEPENDENCIES:
%       
%       This function calls fcn_findEdge_plotLIDARLLA within it.
%
% EXAMPLES:
%   See the script:
%   
%   script_test_fcn_findEdge_classifyGridsBasedOnScanLines
%
%   for a full test suite
%
% This function was written by Aneesh Batchu
%
% REVISION HISTORY:
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

%Does user want to enter format_1?
format_1=sprintf(' ''.'',''Color'',[0.8, 0.8, 0.8],''MarkerSize'', 25,''DisplayName'', "Grids with one scan line" ');
if (4<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
    format_1=temp;
    end
end

%Does user want to enter format_2?
format_2=sprintf(' ''.'',''Color'',[0.2, 0.2, 0.2],''MarkerSize'', 25,''DisplayName'', "Grids with more than one scan line" ');
if (5<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
    format_2=temp;
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

% These grids contain more than one scan line
original_grids_with_more_than_one_scan_line = grids_greater_than_zero_points(total_scan_lines_in_each_grid_with_more_than_zero_points>1);
grid_indices_with_more_than_one_scan_line = (total_scan_lines_in_each_grid_with_more_than_zero_points>1);

% These grids contain only one scan line
original_grids_with_one_scan_line = grids_greater_than_zero_points(~grid_indices_with_more_than_one_scan_line);

% Current grid numbers of the grids with more than one scan line
% These are the grid numbers of the grids with respect to grids_greater_than_zero_points
% The numbering starts with "1"

% Grid centers of the grids with more than zero points and more than one scan line
gridCenters_with_more_than_one_scan_line = gridCenters(original_grids_with_more_than_one_scan_line,1:2);

% Grid centers of the grids with one scan line
gridCenters_with_one_scan_line = gridCenters(original_grids_with_one_scan_line,1:2);

plot_gridCenters_with_one_scan_line = [gridCenters_with_one_scan_line, zeros(length(gridCenters_with_one_scan_line(:,1)),1)]; 

plot_gridCenters_with_more_than_one_scan_line = [gridCenters_with_more_than_one_scan_line, zeros(length(gridCenters_with_more_than_one_scan_line(:,1)),1)]; 

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

    %plot_gridCenters_with_one_scan_line

    fcn_findEdge_plotLIDARLLA(plot_gridCenters_with_one_scan_line,[],[],[],[],[],format_1,fig_num)
    hold on

    %plot_gridCenters_with_one_scan_line
    fcn_findEdge_plotLIDARLLA(plot_gridCenters_with_more_than_one_scan_line,[],[],[],[],[],format_2,fig_num)

    legend()
end


end

