function [zone_start_indices, zone_end_indices] = ...
    fcn_Laps_findSegmentZoneStartStop(...
    query_path,...
    segment_definition,...
    varargin)

% fcn_Laps_findSegmentZoneStartStop
% Given a path and a segment zone definition defined by a start point and
% end point, finds the indices in the path for each crossing of the
% segment, noting the points before and after each crossing as the start
% and end respectively.
%
% A segment is defined by the start and end points in the format of a 2x2
% array with the following entries:
%
%  [X_start Y_start; X_end Y_end];
% 
% Each time a path crosses the segment, the start and end of the path
% crossing are returned by this function. These are defined as follows:
%    * The start of the zone is the index of the first point strictly
%    before the segment crossing.
%    * The end of the zone is the point after the crossing, immediately
%    after the previous start
% 
% Each crossing is direction sensitive, where the crossing must be in the
% positive cross-product direction from the path to the segment. Thus, in
% the positive travel direction along a path, a segment line starts to the
% right of the path and ends in the left of the path for the crossing to be
% positive and be counted.
%
% If a path crosses through a zone repeatedly, each start/end is
% recorded for each path through the zone as another row. Thus, if a path
% crosses through the zone M times (and each time meets the criteria), the
% start, end, and minimum indices will be an M x 1 column.
%
% FORMAT:
%
%      [zone_start_indices, zone_end_indices] = ...
%      fcn_Laps_findPointZoneStartStopAndMinimum(...
%      query_path,...
%      segment_definition,...
%      (fig_num))
%
% INPUTS:
%
%      query_path: the path, in format [X Y] where the matrix is [N by 2].
%      The path definition here is consistent with the Paths library of
%      functions. Note: the function does not yet support 3D paths but can
%      easily be modified for this.
%
%      segment_definition: the segment defining the crossing criteria
%      defined by the start and end points in the format of a 2x2 array
%      with the following entries:  [X_start Y_start; X_end Y_end];
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      zone_start_indices: an array of [M x 1] of the indices where each of
%      the M zones starts. If no zones are detected, the array is empty.
%
%      zone_end_indices: an array of [M x 1] of the indices where each of
%      the M zones ends. If no zones are detected, the array is empty.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_convertPathToTraversalStructure
%      fcn_DebugTools_debugPrintStringToNCharacters
%
% EXAMPLES:
%
%     See the script: script_test_fcn_Laps_findSegmentZoneStartStop
%     for a full test suite.
%
% This function was written on 2022_07_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     
%     2022_07_12: 
%     -- wrote the code originally 
%     2022_11_10: 
%     -- fixed bug in plotting

% TO DO
% 

flag_do_debug = 0; % Flag to show the results for debugging
flag_do_plots = 0; % % Flag to plot the final results
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

if flag_check_inputs
    % Are there the right number of inputs?
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end
        
    % Check the query_path input, 2 or 3 columns, 1 or more rows
    fcn_DebugTools_checkInputsToFunctions(query_path, '2or3column_of_numbers',[1 2]);
    
    % Check the segment_definition input, 2 or 3 columns, 2 rows
    fcn_DebugTools_checkInputsToFunctions(segment_definition, '2or3column_of_numbers',[2 2]);
    
end
        
% Does user want to show the plots?
fig_num = [];
if 3 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plots = 1;
    end
end

%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default values
zone_start_indices = [];
zone_end_indices = [];

% Make sure the segment is 2D
if isequal(size(segment_definition),[2 3])
    warning('The function: fcn_Laps_findSegmentZoneStartStop does not yet support 3D zone definitions. The zone, specified as a 3D point, is being flattened into 2D by ignoring the z-axis value.');
    segment_definition = segment_definition(1:2,1:2);
end

% Calculate the intersection points
if ~isempty(fig_num)
    [distance,~, segment_numbers] = ...
        fcn_Path_findProjectionHitOntoPath(...
        query_path,segment_definition(1,:),segment_definition(2,:),...
        2,fig_num);
else
    [distance,~, segment_numbers] = ...
        fcn_Path_findProjectionHitOntoPath(...
        query_path,segment_definition(1,:),segment_definition(2,:),...
        2);
end

% For debugging:
if flag_do_debug
    INTERNAL_print_results(distance,location, segment_numbers);
end


% Check each of the zones to see if they are empty, and if not, whether
% they are of correct length

% Are zones empty?
if ~isempty(distance)  % empty distances?
    
    % Find the point where segments end (its the one after the segment
    % number)
    next_path_numbers = min(length(query_path(:,1)),segment_numbers+1);
    
    % Convert these to vectors, 
    segment_vector = segment_definition(2,:) - segment_definition(1,:);
    vectors_of_path_hits = query_path(next_path_numbers,:)-query_path(segment_numbers,:);
    
    % Make sure each of the intersections is crossed the correct way by
    % checking the cross product
    cross_result_positive = INTERNAL_crossProduct(vectors_of_path_hits,ones(length(vectors_of_path_hits(:,1)),1)*segment_vector)>0;
    good_segment_numbers = segment_numbers(cross_result_positive);
    
    % Figure out which ones are possible repeats (end to end) by seeing
    % when the difference between them is equal to 1. This only happens
    % when there's end-to-end repeats. These are the zone start indices!
    if ~isempty(good_segment_numbers)
        zone_start_indices = good_segment_numbers([2; diff(good_segment_numbers)]~=1);
        zone_end_indices = zone_start_indices+1;
    else
        zone_start_indices = [];
        zone_end_indices = [];
    end

    num_good_zones = length(zone_start_indices);
    if ~isempty(zone_start_indices)
        % For debugging:
        if flag_do_debug
            fprintf(1,'\nStart, end indices for good zones: \nIstart \t Iend \n');
            for ith_index = 1:num_good_zones
                fprintf(1,'%d\t\t %d\t\t %d\n',zone_start_indices(ith_index,1), zone_end_indices(ith_index,1));
            end
        end
        
    end
   
end % Ends if check to see if zones are empty



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
    
    % plot the final XY result
    figure(fig_num);
    clf;
    hold on;
    grid on
    axis equal
    
    % Plot the query path in blue dots
    plot(query_path(:,1),query_path(:,2),'b.-','Markersize',10);

    % Plot the start and end of query in green and red, respectively
    plot(query_path(1,1),query_path(1,2),'go','Markersize',10);

    % Plot the start and end of query in green and red, respectively
    plot(query_path(end,1),query_path(end,2),'ro','Markersize',10);

    
    % Plot the zone segment in green
    h1 = fcn_Laps_plotZoneDefinition(segment_definition,'g-',fig_num);
    set(h1{1},'Markersize',10);
    
    % Plot the results
    if ~isempty(distance)  % empty
        for ith_zone = 1:num_good_zones
            % Plot the zone
            data_to_plot = query_path(...
                zone_start_indices(ith_zone,1):zone_end_indices(ith_zone,1),:);
            h_fig = plot(data_to_plot(:,1),data_to_plot(:,2),'o-','Markersize',15,'Linewidth',3);
            color_value = get(h_fig,'Color'); %#ok<NASGU>
            
            % Plot the minimum
            %             plot(...
            %                 query_path(zone_min_indices(ith_zone,1),1),...
            %                 query_path(zone_min_indices(ith_zone,1),2),...
            %                 'x','Markersize',10,...
            %                 'Color',color_value,'Linewidth',3);
            
        end
    end
        
    
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

function INTERNAL_print_results(distances,location,segment_number) %#ok<DEFNU>
N_chars = 25;

% Print the header
header_1_str = sprintf('%s','Intersection #');
fixed_header_1_str = fcn_DebugTools_debugPrintStringToNCharacters(header_1_str,N_chars);
header_2_str = sprintf('%s','s-coord along segment');
fixed_header_2_str = fcn_DebugTools_debugPrintStringToNCharacters(header_2_str,N_chars);
header_3_str = sprintf('%s','Location X');
fixed_header_3_str = fcn_DebugTools_debugPrintStringToNCharacters(header_3_str,N_chars);
header_4_str = sprintf('%s','Location Y');
fixed_header_4_str = fcn_DebugTools_debugPrintStringToNCharacters(header_4_str,N_chars);
header_5_str = sprintf('%s','Path segment hit');
fixed_header_5_str = fcn_DebugTools_debugPrintStringToNCharacters(header_5_str,N_chars);

fprintf(1,'\n\n%s %s %s %s %s\n',...
    fixed_header_1_str,...
    fixed_header_2_str,...
    fixed_header_3_str,...
    fixed_header_4_str,...
    fixed_header_5_str);

% Print the results
if ~isempty(distances)
    for ith_intersection =1:length(distances(:,1))
        results_1_str = sprintf('%.0d',ith_intersection);
        fixed_results_1_str = fcn_DebugTools_debugPrintStringToNCharacters(results_1_str,N_chars);
        results_2_str = sprintf('%.2f',distances(ith_intersection,1));
        fixed_results_2_str = fcn_DebugTools_debugPrintStringToNCharacters(results_2_str,N_chars);
        results_3_str = sprintf('%.2f',location(ith_intersection,1));
        fixed_results_3_str = fcn_DebugTools_debugPrintStringToNCharacters(results_3_str,N_chars);
        results_4_str = sprintf('%.2f',location(ith_intersection,2));
        fixed_results_4_str = fcn_DebugTools_debugPrintStringToNCharacters(results_4_str,N_chars);
        results_5_str = sprintf('%.0d',segment_number(ith_intersection));
        fixed_results_5_str = fcn_DebugTools_debugPrintStringToNCharacters(results_5_str,N_chars);
        
        fprintf(1,'%s %s %s %s %s\n',...
            fixed_results_1_str,...
            fixed_results_2_str,...
            fixed_results_3_str,...
            fixed_results_4_str,...
            fixed_results_5_str);
    
    end % Ends for loop
else
    fprintf(1,'(no intersections detected)\n');
end % Ends check to see if isempty
end % Ends function

%% Calculate cross products
function result = INTERNAL_crossProduct(v,w)
result = v(:,1).*w(:,2)-v(:,2).*w(:,1);
end