function h_plot = fcn_Laps_plotSegmentZoneDefinition(zone_definition,varargin)
% fcn_Laps_plotSegmentZoneDefinition
% Plots the segment zone definition given by the input variable, zone_definition.
%
% Accepts as an optional second input standard plot style specifications,
% for example 'r-' for a red line.
%
% As optional third input, plots this in a user-specified figure. Returns
% the plot handle as the output. 
%
% FORMAT: 
%
%       h_plot = fcn_Laps_plotSegmentZoneDefinition(zone_definition,{fig_num})
%
% INPUTS:
%
%      zone_definition: the definition of the zone as given in a 2x2 matrix
%      with the following format [Xstart Ystart; Xend Yend], where X,Y
%      denote the coordinates of the starting and ending point
%      respectively.
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
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Laps_plotSegmentZoneDefinition.m for
%       a full test suite.
%
% This function was written on 2022_04_10 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2022_04_10 
%     -- wrote the code
%     2022_04_12
%     -- fixed minor bug with showing the arrowhead
%     2022_07_23
%     -- fixed minor bug with arrow being different color than line segment
%     -- fixed major bug with arrow being in the wrong direction (!)


flag_do_debug = 0; % Flag to plot the results for debugging
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
    
    % Check the zone_definition input that it has 2 columns, and exactly 2
    % rows
    fcn_DebugTools_checkInputsToFunctions(zone_definition, '2column_of_numbers',[2 2])
end

% Check for plot style input
flag_plot_style_is_specified = 0; % Set default flag value
if 2 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        plot_style = temp;
        flag_plot_style_is_specified = 1;
    end
end

% Does user want to show the plots?
if 3 == nargin
    fig_num = varargin{2};
    figure(fig_num);
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

figure(fig_num);
axis equal;
grid on;

% Check to see if hold is already on. If it is not, set a flag to turn it
% off after this function is over so it doesn't affect future plotting
flag_shut_hold_off = 0;
if ~ishold
    flag_shut_hold_off = 1;
    hold on
end

% Fill in the key data
Xstart = zone_definition(1,1);
Ystart = zone_definition(1,2);
Xend   = zone_definition(2,1);
Yend   = zone_definition(2,2);

% Calculate the unit vector
start_point = mean(zone_definition,1);
magnitude = sum((zone_definition(2,:)-zone_definition(1,:)).^2,2).^0.5;
unit_direction = (zone_definition(2,:)-zone_definition(1,:))./magnitude*[cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)];
vector_with_magnitude = unit_direction*magnitude/2; % Use smaller ratio to make it look good.

% Plot the line
if flag_plot_style_is_specified
    h_plot{1} = plot([Xstart Xend],[Ystart Yend],plot_style,'Markersize',20,'Linewidth',2);
else
    h_plot{1} = plot([Xstart Xend],[Ystart Yend],'-','Linewidth',2);
    main_color = get(h_plot{1},'Color');
end

% Plot the direction arrow
U = vector_with_magnitude(1,1);
V = vector_with_magnitude(1,2);
if flag_plot_style_is_specified
    h_plot{2} = quiver(start_point(1,1),start_point(1,2),U,V,0,plot_style,'Linewidth',2,'Markersize',20,'ShowArrowHead','on','MaxHeadSize',2);
else
    h_plot{2} = quiver(start_point(1,1),start_point(1,2),U,V,0,'color',main_color,'Linewidth',2,'Markersize',20,'ShowArrowHead','on','MaxHeadSize',2);
end



% Shut the hold off?
if flag_shut_hold_off
    hold off;
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

