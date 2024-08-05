
%% Introduction to and Purpose of the Code
% This is the explanation of the code that can be found by running
%       script_demo_Laps.m
% This is a script to demonstrate the functions within the Laps code
% library. This code repo is typically located at:
%   https://github.com/ivsg-psu/FeatureExtraction_DataClean_BreakDataIntoLaps
%
% If you have questions or comments, please contact Sean Brennan at
% sbrennan@psu.edu
%
% The purpose of the code is to break data into "laps", namely portions of
% data defined by start and end points, and in some cases, even allowing
% excursion points that must be "hit" between start and end points. The
% reason for this code is that it is very common that data collection in
% the field passes repeatedly over a test area, even in one data set, and
% thus one must be able to quickly break the code into individual data
% groups with one grouping, or "lap", per traversal.

%% Prep workspace
clc % Clear the console
close all % Close all figures

%% Dependencies and Setup of the Code
% The code requires several other libraries to work, namely the following
%%
% 
% * DebugTools - the repo can be found at: https://github.com/ivsg-psu/Errata_Tutorials_DebugTools
% * PathClassLibrary - the repo can be found at: https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary
% 
% Each is automatically installed in a folder called "Utilities" under the
% root folder, namely ./Utilities/DebugTools/ ,
% ./Utilities/PathClassLibrary/ .
% 
% For ease of transfer, zip files of the directories used - without the
% .git repo information, to keep them small - are referenced and are NOT
% included in this repo. These dependencies are to open code repos, and
% this code accesses these and downloads, if needed, the appropriate
% releases.


% 
% The following code checks to see if the folders flag has been
% initialized, and if not, it calls the DebugTools function that loads the
% path variables. It then loads the PathClassLibrary functions as well.
% Note that the PathClass Library also has sub-utilities that are included.

% USE THE FOLLOWING CODE TO ALLOW MANUAL INSTALLS OF LIBRARIES (IF PRIVATE)
% if ~exist('flag_Laps_Folders_Initialized','var')
%     
%     % add necessary directories for function creation utility 
%     %(special case because folders not added yet)
%     debug_utility_folder = fullfile(pwd, 'Utilities', 'DebugTools');
%     debug_utility_function_folder = fullfile(pwd, 'Utilities', 'DebugTools','Functions');
%     debug_utility_folder_inclusion_script = fullfile(pwd, 'Utilities', 'DebugTools','Functions','fcn_DebugTools_addSubdirectoriesToPath.m');
%     if(exist(debug_utility_folder_inclusion_script,'file'))
%         current_location = pwd;
%         cd(debug_utility_function_folder);
%         fcn_DebugTools_addSubdirectoriesToPath(debug_utility_folder,{'Functions','Data'});
%         cd(current_location);
%     else % Throw an error?
%         error('The necessary utilities are not found. Please add them (see README.md) and run again.');
%     end
%     
%     % Now can add the Path Class Library automatically
%     utility_folder_PathClassLibrary = fullfile(pwd, 'Utilities', 'PathClassLibrary');
%     fcn_DebugTools_addSubdirectoriesToPath(utility_folder_PathClassLibrary,{'Functions','Utilities'});
%     
%     % utility_folder_GetUserInputPath = fullfile(pwd, 'Utilities', 'GetUserInputPath');
%     % fcn_DebugTools_addSubdirectoriesToPath(utility_folder_GetUserInputPath,{'Functions','Utilities'});
% 
%     % Now can add all the other utilities automatically
%     folder_LapsClassLibrary = fullfile(pwd);
%     fcn_DebugTools_addSubdirectoriesToPath(folder_LapsClassLibrary,{'Functions'});
% 
%     % set a flag so we do not have to do this again
%     flag_Laps_Folders_Initialized = 1;
% end

% Use the following code to install public libraries

% List what libraries we need, and where to find the codes for each
clear library_name library_folders library_url

ith_library = 1;
library_name{ith_library}    = 'DebugTools_v2023_01_29';
library_folders{ith_library} = {'Functions','Data'};
library_url{ith_library}     = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_29.zip?raw=true';

ith_library = ith_library+1;
library_name{ith_library}    = 'PathClass_v2023_02_01';
library_folders{ith_library} = {'Functions'};                                
library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary/blob/main/Releases/PathClass_v2023_02_01.zip?raw=true';

ith_library = ith_library+1;
library_name{ith_library}    = 'GetUserInputPath_v2023_02_01';
library_folders{ith_library} = {''};
library_url{ith_library}     = 'https://github.com/ivsg-psu/PathPlanning_PathTools_GetUserInputPath/blob/main/Releases/GetUserInputPath_v2023_02_01.zip?raw=true';

% Do we need to set up the work space?
if ~exist('flag_Laps_Folders_Initialized','var')

    % Reset all flags for installs to empty
    clear global FLAG*

    fprintf(1,'Installing utilities necessary for code ...\n');

    % Dependencies and Setup of the Code
    % This code depends on several other libraries of codes that contain
    % commonly used functions. We check to see if these libraries are installed
    % into our "Utilities" folder, and if not, we install them and then set a
    % flag to not install them again.
    
    % Set up libraries
    for ith_library = 1:length(library_name)
        dependency_name = library_name{ith_library};
        dependency_subfolders = library_folders{ith_library};
        dependency_url = library_url{ith_library};

        fprintf(1,'\tAdding library: %s ...',dependency_name);
        fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url);
        clear dependency_name dependency_subfolders dependency_url
        fprintf(1,'Done.\n');
    end

    % Set dependencies for this project specifically
    fcn_DebugTools_addSubdirectoriesToPath(pwd,{'Functions','Data'});
    
    disp('Done setting up libraries, adding each to MATLAB path, and adding current repo folders to path.');
    
end

%% Using Zone Definitions to Define Start, End, and Excursion Locations
% To define the start, end, and excursion locations for data, the data must
% pass through or nearby a geolocation which is hereafter called a "zone
% definition". There are two types of zone definitions used in this code:
%%
% 
% * Point methods of zone definitions - this is when a start, stop, or
% excursion is defined by "passing by" a point. For example, if a journey
% is said to start at someone's house and go to someone's office, then the
% location of the house and office define the start and end of the journey.
% The specification is given by an X,Y location and a radius in the form of
% [X Y radius], as a 3x1 matrix. Whenever the path passes within the radius
% with a specified number of points within that radius, the minimum
% distance point then "triggers" the zone.
% * Line segment methods of zone definitions - this when a start, stop, or
% excursion condition is defined by a path passing through a line segment.
% The line segment is given by the X,Y coordinates of the start and stop of
% the line segment, in the form [Xstart Ystart; Xend Yend], thus producing
% a 2x2 matrix. An example of a line segment definition is the start line
% and finish line of a race.
% 
% To illustrate both definitions, we first create some data to plot:

full_steps = (-1:0.1:1)';
zero_full_steps = 0*full_steps; %#ok<NASGU>
ones_full_steps = ones(length(full_steps(:,1)),1);
half_steps = (-1:0.1:0)';
zero_half_steps = 0*half_steps;
ones_half_steps = ones(length(half_steps(:,1)),1); %#ok<PREALL>
path_examples{1} = [-1*ones_full_steps full_steps];
path_examples{2} = [1*ones_full_steps full_steps];

%% 
% Each of the path_example matrices above can be plotted easily using the
% "plotLapsXY" subfunction, but this function expects the paths to be in a
% traversal type so that it is compatible with the Path library of
% functions. To convert them, we use the conversion utility from the Path
% library, convert each to "traversal" types stored in a variable called
% path_data. We then plot the paths.

clear path_data
for i_Path = 1:length(path_examples)
    traversal = fcn_Path_convertPathToTraversalStructure(path_examples{i_Path});
    path_data.traversal{i_Path} = traversal;
end

% Plot the results via fcn_Laps_plotLapsXY
fig_num = 222;
fcn_Laps_plotLapsXY(path_data,fig_num);    

%%
% Now, use a zone plotting tool to show the point and line-segment types of
% zone definitions. The point definition is shown in green, and the segment
% definition is shown in blue. The segment definition includes an arrow
% that points in the direction of an allowable crossing.

fig_num = 444;

zone_center = [0.8 0];
zone_radius = 2;
num_points = 3;
point_zone_definition = [zone_radius num_points zone_center];
fcn_Laps_plotPointZoneDefinition(point_zone_definition,'g',fig_num);

segment_zone_definition = [0.8 0; 1.2 0];
fcn_Laps_plotSegmentZoneDefinition(segment_zone_definition,'b',fig_num);


%%
% Show we can get the same plot now via a combined function

fig_num = 4443;
fcn_Laps_plotZoneDefinition(point_zone_definition,'g',fig_num);
fcn_Laps_plotZoneDefinition(segment_zone_definition,'b',fig_num);

%% Point zone evaluations
% The function, fcn_Laps_findPointZoneStartStopAndMinimum, uses a point
% zone evaluation to determine portions of a segment that are within a
% point zone definition. For example, if the path does not cross into the
% zone, nothing is returned:
fig_num = 1;

query_path = ...
    [full_steps 0.4*ones_full_steps];

zone_center = [0 0 0.2]; % Located at [0,0] 
zone_radius = 0.2; % with radius 0.2
[zone_start_indices, zone_end_indices, zone_min_indices] = ...
    fcn_Laps_findPointZoneStartStopAndMinimum(...
    query_path,...
    zone_center,...
    zone_radius,...
    [],...
    fig_num);

assert(isempty(zone_start_indices));
assert(isempty(zone_end_indices));
assert(isempty(zone_min_indices));

%%
% And, the default is that three points must be within the zone. So, if a
% path only crosses one or two points, then nothing is returned.

fig_num = 2;

query_path = ...
    [full_steps 0.2*ones_full_steps];

zone_center = [0 0 0.2]; % Located at [0,0] 
zone_radius = 0.2; % with radius 0.2
[zone_start_indices, zone_end_indices, zone_min_indices] = ...
    fcn_Laps_findPointZoneStartStopAndMinimum(...
    query_path,...
    zone_center,...
    zone_radius,...
    [],...
    fig_num);

assert(isempty(zone_start_indices));
assert(isempty(zone_end_indices));
assert(isempty(zone_min_indices));


% Show that 2 points still doesn't work
query_path = ...
    [full_steps 0.2*ones_full_steps];

zone_center = [0.05 0 0.2]; % Located at [0.05,0] 
zone_radius = 0.23; % with radius 0.23
[zone_start_indices, zone_end_indices, zone_min_indices] = ...
    fcn_Laps_findPointZoneStartStopAndMinimum(...
    query_path,...
    zone_center,...
    zone_radius,...
    [],...
    fig_num);

assert(isempty(zone_start_indices));
assert(isempty(zone_end_indices));
assert(isempty(zone_min_indices));


%%
% But, if a path crosses the zone with at least three points, then the
% indices of the start, end, and minimum of the path are returned.
fig_num = 3;

query_path = ...
    [half_steps zero_half_steps];

zone_center = [-0.02 0 0.2]; % Located at [-0.02,0] 
zone_radius = 0.2; % with radius 0.2
[zone_start_indices, zone_end_indices, zone_min_indices] = ...
    fcn_Laps_findPointZoneStartStopAndMinimum(...
    query_path,...
    zone_center,...
    zone_radius,...
    [],...
    fig_num);

assert(isequal(zone_start_indices,9));
assert(isequal(zone_end_indices,11));
assert(isequal(zone_min_indices,11));

%%
% If there are multiple crossings of the zone, then indices of the
% start/stop/minimum are returned for each crossing:
full_steps = (-1:0.1:1)';
zero_full_steps = 0*full_steps;
ones_full_steps = ones(length(full_steps(:,1)),1);
half_steps = (-1:0.1:0)';
zero_half_steps = 0*half_steps;
ones_half_steps = ones(length(half_steps(:,1)),1);

minimum_number_of_indices_in_zone = 3;
fig_num = 5;


query_path = ...
    [full_steps 0*ones_full_steps; -full_steps 0.1*ones_full_steps; full_steps 0.2*ones_full_steps ];

zone_center = [0.05 0]; % Located at [0.05,0]
zone_radius = 0.23;
[zone_start_indices, zone_end_indices, zone_min_indices] = ...
    fcn_Laps_findPointZoneStartStopAndMinimum(...
    query_path,...
    zone_center,...
    zone_radius,...
    minimum_number_of_indices_in_zone,...
    fig_num);

assert(isequal(zone_start_indices,[10; 30]));
assert(isequal(zone_end_indices,  [13; 33]));
assert(isequal(zone_min_indices,  [12; 31]));


%% Create sample paths
% To illustrate the functionality of this library, we call the library
% function fillPathViaUserInputs which fills in an array of "path" types.
% Load some test data and plot it in figure 1 

% Call the function to fill in an array of "path" type
laps_array = fcn_Laps_fillSampleLaps;


% Use Path library functions to onvert paths to traversals structures. Each
% traversal instance is a "traversal" type, and the array called "data"
% below is a "traversals" type.
for i_Path = 1:length(laps_array)
    traversal = fcn_Path_convertPathToTraversalStructure(laps_array{i_Path});
    data.traversal{i_Path} = traversal;
end


% Plot all the laps
fig_num = 22323;
for ith_example = 1:length(data.traversal)
    single_lap.traversal{1} = data.traversal{ith_example};
    fcn_Laps_plotLapsXY(single_lap,fig_num);
end

% Plot the last one
fig_num = 1;
single_lap.traversal{1} = data.traversal{end};
fcn_Laps_plotLapsXY(single_lap,fig_num);

%% Show fcn_Laps_plotZoneDefinition.m 
% Plots the zone, allowing user-defined colors. For example, the figure
% below shows a radial zone for the start, and a line segment for the end.
start_definition = [10 3 0 0]; % Radius 10, 3 points must pass near [0 0]
fcn_Laps_plotZoneDefinition(start_definition,'g',fig_num);

end_definition = [40 -40; 80 -40]; % must cross a line segment starting at [40 -40], ending at [80 -40]
fcn_Laps_plotZoneDefinition(end_definition,'r',fig_num);

%% Call the fcn_Laps_breakDataIntoLaps function, plot in figure 2
% Test of fcn_Laps_breakDataIntoLaps.m : This is the core function for this
% repo that breaks data into laps. Note: for radial zone definitions, the
% image illustrates how a lap starts at the first point within a start
% zone, and ends at the last point before exiting the end zone.
start_definition = [10 3 0 0]; % Radius 10, 3 points must pass near [0 0]
end_definition = [30 3 0 -60]; % Radius 30, 3 points must pass near [0,-60]

excursion_definition = []; % empty
fig_num = 2;
lap_traversals = fcn_Laps_breakDataIntoLaps(...
    single_lap.traversal{1},...
    start_definition,...
    end_definition,...
    excursion_definition,...
    fig_num);

% Do we get 3 laps?
assert(isequal(3,length(lap_traversals.traversal)));

% Are the laps different lengths?
assert(isequal(87,length(lap_traversals.traversal{1}.X)));
assert(isequal(98,length(lap_traversals.traversal{2}.X)));
assert(isequal(79,length(lap_traversals.traversal{3}.X)));

%% Call the fcn_Laps_breakDataIntoLapIndices function, plot in figure 3
start_definition = [10 3 0 0]; % Radius 10, 3 points must pass near [0 0]
end_definition = [30 3 0 -60]; % Radius 30, 3 points must pass near [0,-60]
excursion_definition = []; % empty
fig_num = 2;
lap_traversals = fcn_Laps_breakDataIntoLaps(...
    single_lap.traversal{1},...
    start_definition,...
    end_definition,...
    excursion_definition,...
    fig_num);

% Do we get 3 laps?
assert(isequal(3,length(lap_traversals.traversal)));

% Are the laps different lengths?
assert(isequal(87,length(lap_traversals.traversal{1}.X)));
assert(isequal(98,length(lap_traversals.traversal{2}.X)));
assert(isequal(79,length(lap_traversals.traversal{3}.X)));


%% Revision History:
%      2022_03_27:
%      -- created a demo script of core debug utilities
%      2022_04_02
%      -- Added sample path creation
%      2022_04_04
%      -- Added minor edits
%      2022_04_10
%      -- Added comments, plotting utilities for zone definitions
%      2022_05_21
%      -- More cleanup
%      2022_07_23 - sbrennan@psu.edu
%      -- Enable index-based look-up
%      2023_02_01 - sbrennan@psu.edu
%      -- Enable web-based installs


function fcn_INTERNAL_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url, varargin)
%% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES - MATLAB package installer from URL
%
% FCN_DEBUGTOOLS_INSTALLDEPENDENCIES installs code packages that are
% specified by a URL pointing to a zip file into a default local subfolder,
% "Utilities", under the root folder. It also adds either the package
% subfoder or any specified sub-subfolders to the MATLAB path.
%
% If the Utilities folder does not exist, it is created.
% 
% If the specified code package folder and all subfolders already exist,
% the package is not installed. Otherwise, the folders are created as
% needed, and the package is installed.
% 
% If one does not wish to put these codes in different directories, the
% function can be easily modified with strings specifying the
% desired install location.
% 
% For path creation, if the "DebugTools" package is being installed, the
% code installs the package, then shifts temporarily into the package to
% complete the path definitions for MATLAB. If the DebugTools is not
% already installed, an error is thrown as these tools are needed for the
% path creation.
% 
% Finally, the code sets a global flag to indicate that the folders are
% initialized so that, in this session, if the code is called again the
% folders will not be installed. This global flag can be overwritten by an
% optional flag input.
%
% FORMAT:
%
%      fcn_DebugTools_installDependencies(...
%           dependency_name, ...
%           dependency_subfolders, ...
%           dependency_url)
%
% INPUTS:
%
%      dependency_name: the name given to the subfolder in the Utilities
%      directory for the package install
%
%      dependency_subfolders: in addition to the package subfoder, a list
%      of any specified sub-subfolders to the MATLAB path. Leave blank to
%      add only the package subfolder to the path. See the example below.
%
%      dependency_url: the URL pointing to the code package.
%
%      (OPTIONAL INPUTS)
%      flag_force_creation: if any value other than zero, forces the
%      install to occur even if the global flag is set.
%
% OUTPUTS:
%
%      (none)
%
% DEPENDENCIES:
%
%      This code will automatically get dependent files from the internet,
%      but of course this requires an internet connection. If the
%      DebugTools are being installed, it does not require any other
%      functions. But for other packages, it uses the following from the
%      DebugTools library: fcn_DebugTools_addSubdirectoriesToPath
%
% EXAMPLES:
%
% % Define the name of subfolder to be created in "Utilities" subfolder
% dependency_name = 'DebugTools_v2023_01_18';
%
% % Define sub-subfolders that are in the code package that also need to be
% % added to the MATLAB path after install; the package install subfolder
% % is NOT added to path. OR: Leave empty ({}) to only add 
% % the subfolder path without any sub-subfolder path additions. 
% dependency_subfolders = {'Functions','Data'};
%
% % Define a universal resource locator (URL) pointing to the zip file to
% % install. For example, here is the zip file location to the Debugtools
% % package on GitHub:
% dependency_url = 'https://github.com/ivsg-psu/Errata_Tutorials_DebugTools/blob/main/Releases/DebugTools_v2023_01_18.zip?raw=true';
%
% % Call the function to do the install
% fcn_DebugTools_installDependencies(dependency_name, dependency_subfolders, dependency_url)
%
% This function was written on 2023_01_23 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2023_01_23:
% -- wrote the code originally

% TO DO
% -- Add input argument checking

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
    narginchk(3,4);
end

%% Set the global variable - need this for input checking
% Create a variable name for our flag. Stylistically, global variables are
% usually all caps.
flag_varname = upper(cat(2,'flag_',dependency_name,'_Folders_Initialized'));

% Make the variable global
eval(sprintf('global %s',flag_varname));

if nargin==4
    if varargin{1}
        eval(sprintf('clear global %s',flag_varname));
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



if ~exist(flag_varname,'var') || isempty(eval(flag_varname))
    % Save the root directory, so we can get back to it after some of the
    % operations below. We use the Print Working Directory command (pwd) to
    % do this. Note: this command is from Unix/Linux world, but is so
    % useful that MATLAB made their own!
    root_directory_name = pwd;

    % Does the directory "Utilities" exist?
    utilities_folder_name = fullfile(root_directory_name,'Utilities');
    if ~exist(utilities_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(root_directory_name,'Utilities');

        % Did it work?
        if ~success_flag
            error('Unable to make the Utilities directory. Reason: %s with message ID: %s\n',error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The Utilities directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',error_message, message_ID);
        end

    end

    % Does the directory for the dependency folder exist?
    dependency_folder_name = fullfile(root_directory_name,'Utilities',dependency_name);
    if ~exist(dependency_folder_name,'dir')
        % If we are in here, the directory does not exist. So create it
        % using mkdir
        [success_flag,error_message,message_ID] = mkdir(utilities_folder_name,dependency_name);

        % Did it work?
        if ~success_flag
            error('Unable to make the dependency directory: %s. Reason: %s with message ID: %s\n',dependency_name, error_message,message_ID);
        elseif ~isempty(error_message)
            warning('The %s directory was created, but with a warning: %s\n and message ID: %s\n(continuing)\n',dependency_name, error_message, message_ID);
        end

    end

    % Do the subfolders exist?
    flag_allFoldersThere = 1;
    if isempty(dependency_subfolders)
        flag_allFoldersThere = 0;
    else
        for ith_folder = 1:length(dependency_subfolders)
            subfolder_name = dependency_subfolders{ith_folder};
            
            % Create the entire path
            subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
            
            % Check if the folder and file exists that is typically created when
            % unzipping.
            if ~exist(subfunction_folder,'dir')
                flag_allFoldersThere = 0;
            end
        end
    end

    % Do we need to unzip the files?
    if flag_allFoldersThere==0
        % Files do not exist yet - try unzipping them.
        save_file_name = tempname(root_directory_name);
        zip_file_name = websave(save_file_name,dependency_url);
        % CANT GET THIS TO WORK --> unzip(zip_file_url, debugTools_folder_name);

        % Is the file there?
        if ~exist(zip_file_name,'file')
            error('The zip file: %s for dependency: %s did not download correctly. This is usually because permissions are restricted on the current directory. Check the code install (see README.md) and try again.\n',zip_file_name, dependency_name);
        end

        % Try unzipping
        unzip(zip_file_name, dependency_folder_name);

        % Did this work?
        flag_allFoldersThere = 1;
        if ~isempty(dependency_subfolders)
            for ith_folder = 1:length(dependency_subfolders)
                subfolder_name = dependency_subfolders{ith_folder};
                
                % Create the entire path
                subfunction_folder = fullfile(root_directory_name, 'Utilities', dependency_name,subfolder_name);
                
                % Check if the folder and file exists that is typically created when
                % unzipping.
                if ~exist(subfunction_folder,'dir')
                    flag_allFoldersThere = 0;
                end
            end
        end
        
        if flag_allFoldersThere==0
            error('The necessary dependency: %s has an error in install, or error performing an unzip operation. Check the code install (see README.md) and try again.\n',dependency_name);
        else
            % Clean up the zip file
            delete(zip_file_name);
        end

    end


    % For path creation, if the "DebugTools" package is being installed, the
    % code installs the package, then shifts temporarily into the package to
    % complete the path definitions for MATLAB. If the DebugTools is not
    % already installed, an error is thrown as these tools are needed for the
    % path creation.
    %
    % In other words: DebugTools is a special case because folders not
    % added yet, and we use DebugTools for adding the other directories
    if strcmp(dependency_name(1:10),'DebugTools')
        debugTools_function_folder = fullfile(root_directory_name, 'Utilities', dependency_name,'Functions');

        % Move into the folder, run the function, and move back
        cd(debugTools_function_folder);
        fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        cd(root_directory_name);
    else
        try
            fcn_DebugTools_addSubdirectoriesToPath(dependency_folder_name,dependency_subfolders);
        catch
            error('Package installer requires DebugTools package to be installed first. Please install that before installing this package');
        end
    end


    % Finally, the code sets a global flag to indicate that the folders are
    % initialized.  Check this using a command "exist", which takes a
    % character string (the name inside the '' marks, and a type string -
    % in this case 'var') and checks if a variable ('var') exists in matlab
    % that has the same name as the string. The ~ in front of exist says to
    % do the opposite. So the following command basically means: if the
    % variable named 'flag_CodeX_Folders_Initialized' does NOT exist in the
    % workspace, run the code in the if statement. If we look at the bottom
    % of the if statement, we fill in that variable. That way, the next
    % time the code is run - assuming the if statement ran to the end -
    % this section of code will NOT be run twice.

    eval(sprintf('%s = 1;',flag_varname));
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

    % Nothing to do!



end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function fcn_DebugTools_installDependencies