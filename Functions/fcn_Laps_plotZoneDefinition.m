function h_plot = fcn_Laps_plotZoneDefinition(zone_definition,varargin)
% fcn_Laps_plotZoneDefinition - Plots the zone definition for Laps codes
%
% Accepts as an optional second input standard plot style specifications,
% for example 'r-' for a red line.
%
% As optional third input, plots this in a user-specified figure. Returns
% the plot handle as the output. 
%
% FORMAT: 
%
%       h_plot = fcn_Laps_plotZoneDefinition(zone_definition,{fig_num})
%
% INPUTS:
%
%      zone_definition: the definition of the zone as given in a point zone
%      or segment zone style. See fcn_Laps_breakDataIntoLapIndices for
%      details.
%
%      (OPTIONAL INPUTS)
%   
%      plot_style: the standard plot pecification style allowing line and
%      color, for example 'r-'. Type "help plot" for a listing of options.
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      h: a handle to the resulting figure
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Laps_checkZoneType
%      fcn_Laps_plotSegmentZoneDefinition
%      fcn_Laps_plotPointZoneDefinition
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Laps_plotZoneDefinition.m for
%       a full test suite.
%
% This function was written on 2022_07_23 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2022_07_23 
%     -- wrote the code


flag_do_debug = 0; % Flag to plot the results for debugging
flag_this_is_a_new_figure = 1; %#ok<NASGU> % Flag to check to see if this is a new figure
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 1 || nargin > 3
        error('Incorrect number of input arguments')
    end
          
    % note: zone definition type is checked within code below
    
end

% Check for plot style input
plot_style = []; % Leave it empty for defaults
if 2 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        plot_style = temp;
    end
end

% Does user want to show the plots?
if 3 == nargin
    temp = varargin{2};
    if ~isempty(temp)
        fig_num = temp;
        figure(fig_num);
        flag_this_is_a_new_figure = 0; %#ok<NASGU>
    end
else    
    fig = figure;
    fig_num = fig.Number;
end

%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% What type is it?
[flag_is_a_point_zone_type, new_zone_definition] = fcn_Laps_checkZoneType(zone_definition, 'plot_zone_definition');

% Is it a point-zone type?
if flag_is_a_point_zone_type
    h_plot = fcn_Laps_plotPointZoneDefinition(new_zone_definition,plot_style,fig_num);
else % No, it's a segment zone
    h_plot = fcn_Laps_plotSegmentZoneDefinition(new_zone_definition,plot_style,fig_num);
end

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
if flag_do_debug
    % Nothing in here yet
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end

