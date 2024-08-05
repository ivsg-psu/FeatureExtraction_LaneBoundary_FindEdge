function h = fcn_Laps_plotLapsXY(traversals,varargin)
% fcn_Laps_plotLapsXY
% Plots the XY positions of all laps existing in a data structure. This is
% just a modified version of fcn_Path_plotTraversalsXY from the "Path"
% library of functions. The "laps" type and "tranversals" type are the same
% type, for consistency.
%
% FORMAT: 
%
%       h = fcn_Laps_plotLapsXY(traversals,{fig_num})
%
% INPUTS:
%
%      data: a structure containing subfields of station and Yaw in the
%      following form
%           data.traversal{ith_lap}.X
%           data.traversal{ith_lap}.Y
%      Note that ith_lap denotes an array of laps. Each lap will be
%      plotted separately.
%
%      (OPTIONAL INPUTS)
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
%       See the script: script_test_fcn_Laps_plotLapsXY.m for a full test
%       suite. 
%
% This function was written on 2022_04_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2022_04_02 
%     -- wrote the code


flag_do_debug = 0; % Flag to plot the results for debugging
flag_this_is_a_new_figure = 1; % Flag to check to see if this is a new figure
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
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end
    
    % Check the data input
    % fcn_Path_checkInputsToFunctions(traversals, 'traversals');

end

% Does user want to show the plots?
if 2 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_this_is_a_new_figure = 0;
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

% Check to see if hold is already on. If it is not, set a flag to turn it
% off after this function is over so it doesn't affect future plotting
flag_shut_hold_off = 0;
if ~ishold
    flag_shut_hold_off = 1;
    hold on
end

NumTraversals = length(traversals.traversal);
h = zeros(NumTraversals,1);
for ith_lap= 1:NumTraversals
    if ~isempty(traversals.traversal{ith_lap})
        h(ith_lap) = plot(traversals.traversal{ith_lap}.X,traversals.traversal{ith_lap}.Y,'-o');
    else
        h(ith_lap) = plot(NaN,NaN,'-o');
    end
end

% Plot the start and end values as green and red respectively
for ith_lap= 1:NumTraversals
    if ~isempty(traversals.traversal{ith_lap})
        plot(...
            traversals.traversal{ith_lap}.X(1,1),traversals.traversal{ith_lap}.Y(1,1),'go',...
            traversals.traversal{ith_lap}.X(end,1),traversals.traversal{ith_lap}.Y(end,1),'ro');
    else
        plot(NaN,NaN,'go',NaN,NaN,'ro');
    end
end

% Shut the hold off?
if flag_shut_hold_off
    hold off;
end

% Add labels? 
if flag_this_is_a_new_figure == 1
    title('X vs Y')
    xlabel('X [m]')
    ylabel('Y [m]')
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

