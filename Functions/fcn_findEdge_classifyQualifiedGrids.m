function [current_qualified_grids,current_unqualified_grids,original_qualified_grids,gridCenters_qualified_grids,gridCenters_unqualified_grids] = ...
    fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,grid_indices_with_more_than_transverse_span_threshold,...
    grids_greater_than_zero_points,gridCenters,varargin)
%% fcn_findEdge_classifyQualifiedGrids
%
% FORMAT:
%   [original_qualified_grids,gridCenters_qualified_grids,gridCenters_unqualified_grids] = ...
%   fcn_findEdge_classifyQualifiedGrids(grid_indices_with_required_point_density,grid_indices_with_more_than_one_scan_line,...
%   grid_indices_with_more_than_transverse_span_threshold, grids_greater_than_zero_points,gridCenters,(format_unqualified),(format_qualified),(fig_num))
%
%
% INPUTS:
%   grid_indices_with_required_point_density
%   grid_indices_with_more_than_one_scan_line
%   grid_indices_with_more_than_transverse_span_threshold
%   grids_greater_than_zero_points
%   gridCenters
%
% (OPTIONAL INPUTS)
%   
%   format_unqualified
%   
%   format_qualified
%   
%   fig_num
%
% OUTPUTS:
%   
%   current_qualified_grids
%   current_unqualified_grids 
%   original_qualified_grids
%   gridCenters_qualified_grids
%   gridCenters_unqualified_grids
%
% DEPENDENCIES:
%
% EXAMPLES:
%   See the script:
%   
%   script_test_fcn_findEdge_classifyQualifiedGrids
%
%   for a full test suite
%
% This function was written by Aneesh Batchu
%
% REVISION HISTORY:
% 2024_08_14 - Aleksandr Goncharov
% -- Functionalized code from the FindEdge Demo 

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
        narginchk(5,8);
    end
end 

% Does the user want to specify format_unqualified?
format_unqualified=sprintf(' ''.'',''Color'',[0.8, 0.8, 0.8],''MarkerSize'', 25,''DisplayName'',''Unqualified Grids''); (legend ');
if (6<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
    format_unqualified=temp;
    end
end


% Does the user want to specify format_qualified?
format_qualified=sprintf(' ''.'',''Color'',[0.2, 0.2, 0.2],''MarkerSize'', 25,''DisplayName'',''Qualified grids''); (legend ');
if (7<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
    format_qualified=temp;
    end
end


% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (8<=nargin)
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


% Grid indices of original grids with required point density, more than one
% scan line, and more than transverse span threshold
grid_indices_qualified_grids = (grid_indices_with_required_point_density == 1) & (grid_indices_with_more_than_one_scan_line == 1) & (grid_indices_with_more_than_transverse_span_threshold == 1);

% --------------------------Qualified grid numbers-----------------------------------
% The grids with more than one scan line, more/equal to point density and
% greater than minimum transverse span threshold
original_qualified_grids = grids_greater_than_zero_points(grid_indices_qualified_grids); 

% --------------------------Unqualified grid numbers-----------------------------------
% Find grids with low point density, grids with only one LiDAR scan and
% lesser than minimum transverse span threshold
% The grids that are not qualified grids are clasified as unqualified grids
original_unqualified_grids = grids_greater_than_zero_points(~grid_indices_qualified_grids); 

% Current qualified grid numbers 
current_qualified_grids = find(grid_indices_qualified_grids); 

% Current unqualified grid numbers 
current_unqualified_grids = find(~grid_indices_qualified_grids); 

% Grid centers of qualified grids
gridCenters_qualified_grids = gridCenters(original_qualified_grids,1:2); 

% Grid centers of unqualified grids
gridCenters_unqualified_grids = gridCenters(original_unqualified_grids,1:2); 

%

plot_gridCenters_unqualified_grids = [gridCenters_unqualified_grids, zeros(length(gridCenters_unqualified_grids(:,1)),1)]; 

plot_gridCenters_qualified_grids = [gridCenters_qualified_grids, zeros(length(gridCenters_qualified_grids(:,1)),1)]; 

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
   

    fcn_findEdge_plotLIDARLLA(plot_gridCenters_unqualified_grids,[],[],[],[],[],format_unqualified,fig_num);

    hold on

    fcn_findEdge_plotLIDARLLA(plot_gridCenters_qualified_grids,[],[],[],[],[],format_qualified,fig_num);

end

end
