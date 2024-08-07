function [VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
    LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
    fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, varargin)
% fcn_findEdge_extractScanLines   extract LIDAR scan lines associated with a vehicle pose location.
% 
% FORMAT:
%
%      [VehiclePose_ENU, VehiclePose_UnitOrthoVectors, ...
%      LIDAR_ENU, LIDAR_intensity, LIDAR_scanLineAndRingID] = ...
%      fcn_findEdge_extractScanLines(VehiclePose, LiDAR_Scan_ENU_Entire_Loop, (scanLineRange), (ringsRange), (fig_num))
%
% INPUTS:  
%
%      VehiclePose: the position of the mapping vehicle during mapping
%
%      LiDAR_Scan_ENU_Entire_Loop: the points from the LIDAR scan 
%      
%      (OPTIONAL INPUTS)
%       
%      scanLineRange: the [min_number max_number] of the scanlines to
%      extract. The default is to use all scan lines.
%
%      ringsRange: a vector of form [L1 L2 L3 ... LN] listing the
%      rings of the laser to analyze. The default is to use [1 ... N] where
%      N is the total number of rings.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.     
%
% OUTPUTS:
%
%      VehiclePose_ENU: the position of the mapping vehicle during mapping
%      in ENU coordinates as an [Nx3] vector. The first column is
%      EAST: the direction pointing towards the geographic east. The second
%      column is NORTH: the direction pointing towards the geographic
%      north. The third column is UP (U): the direction pointing vertically
%      upwards, perpendicular to the Earth's surface.
%
%      VehiclePose_UnitOrthoVectors: Orthogonal unit vectors of the
%      position of the mapping vehicle during mapping as a [Nx3] vector,
%      with one orthogonal vector per position. The first two columns are
%      the X and Y components of the unit orthogonal vectors, e.g. the
%      direction defined as positive - at that location - in the transverse
%      or "left" direction. The third column contains all zero unit vectors
%      because the elevation of the mapping vehicle is not considered when
%      calculating whether objects are to the right/left of the vehicle.
%      Thus, this output is used typically to determine whether coordinates
%      are to the right or left of the vehicle.
%
%      LIDAR_ENU: The converted data collected by the LIDAR on the mapping 
%      vehicle in ENU coordinates as an [Nx3] vector. The format is the
%      same as VehiclePose_ENU.
%
%      LIDAR_intensity: The intensity of the LIDAR data during mapping as a
%      [Nx1] vector. Higher intensity value means the target surface is
%      more reflective.
%
%      LIDAR_scanLineAndRingID: The scanline and ringID during mapping as a
%      [Nx2] vector. The first column is the scanline, which indicates
%      which scan instance (rotation of the LIDAR sensor) created that
%      specific LIDAR data. The second column is the ringID which denotes
%      which specific laser beam, in the array of lasers, created that data
%      point.
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
% 
%       script_test_fcn_findEdge_extractScanLines
%  
%       for a full test suite.
%
% This function was written on 2024_08_05 by Sean Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history
% 2024_08_05 - Sean Brennan
% -- Created function by copying out of load script in Geometry library
% 2024_08_06 - Jiabao Zhao
% -- Documented the output variables.


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
        narginchk(2,5);

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


% Does user want to specify scanLineRange?
scanLineRange = []; 
if (1<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        scanLineRange = temp;
    end
end

% Does user want to specify rings_to_analyze?
ringsRange = []; 
if (2<=nargin)
    temp = varargin{2};
    if ~isempty(temp)
        ringsRange = temp;
    end
end


% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (2<=nargin)
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


%% Make sure the vehicle pose is EXACTLY the same length as the LIDAR data
assert(length(LiDAR_Scan_ENU_Entire_Loop)==length(VehiclePose(:,1)));
N_scanLines = length(VehiclePose(:,1));

%% Calculate the vehicle orientation
vehicle_change_in_pose_XY = diff(VehiclePose(:,1:2));

% Repeat the last value again, since diff removes one row. We want the same
% number of vectors as the number of points, and diff removed one point.
vehicle_change_in_pose_XY = [vehicle_change_in_pose_XY; vehicle_change_in_pose_XY(end,:)];

% Convert these to unit vectors
unit_vehicle_change_in_pose_XY = fcn_INTERNAL_calcUnitVector(vehicle_change_in_pose_XY);

% Find orthogonal vetors by rotating by 90 degrees in the CCW direction by
% multiplying by the rotatio matrix
unit_ortho_vehicle_vectors_XY = unit_vehicle_change_in_pose_XY*[0 1; -1 0];

%% Get the scanLineRange
if ~isempty(scanLineRange)
    scanLineNumber_start = scanLineRange(1);
    scanLineNumber_end   = scanLineRange(2);
    if scanLineNumber_start<1
        warning('on','backtrace');
        warning('Incorrect scan line detected');
        error('Starting scan line must be greater than or equal to 1');
    end
    if scanLineNumber_start>N_scanLines
        warning('on','backtrace');
        warning('Incorrect scan line detected');
        error('Starting scan line must be less than or equal to number of scans: %.0f',N_scanLines);
    end

else
    scanLineNumber_start = 1;
    scanLineNumber_end   = N_scanLines;
end

%% Get the ringsRange
if ~isempty(ringsRange)
    if min(ringsRange)<1
        warning('on','backtrace');
        warning('Incorrect ring range detected');
        error('Lowest ring range index must be greater than or equal to 1');
    end
    if max(ringsRange)>15
        warning('on','backtrace');
        warning('Incorrect ring range detected');
        error('Highest ring range index must be less than or equal to 15');
    end
else
    ringsRange = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
end

%% Find largest scan data
max_lines = -inf;
for ith_Scan = 1:N_scanLines
    max_lines = max(max_lines,length(LiDAR_Scan_ENU_Entire_Loop{ith_Scan}(:,1)));    
end
max_lines = max_lines+16; % There are 16 rings that may have NaN values at end of each

%% Load all the data (only do this ONCE)
% persistent allScanLineNumbers allVehiclePose_ENU allVehiclePose_UnitOrthoVectors allLIDAR_ENU allLIDAR_intensity allLIDAR_scanLineAndRingID
% 
% if isempty(allScanLineNumbers) || isempty(allVehiclePose_ENU) || isempty(allVehiclePose_UnitOrthoVectors) || isempty(allLIDAR_ENU) || isempty(allLIDAR_intensity) || isempty(allLIDAR_scanLineAndRingID)

totalDataLength = N_scanLines*max_lines;

persistent allData alltotalDataLength

if isempty(allData) || isempty(alltotalDataLength) || ~isequal(totalDataLength,alltotalDataLength)
    % ScanLineNumbers              = nan(N_scanLines*max_lines,1);
    % VehiclePose_ENU              = nan(N_scanLines*max_lines,3);
    % VehiclePose_UnitOrthoVectors = nan(N_scanLines*max_lines,3);
    % LIDAR_ENU                    = nan(N_scanLines*max_lines,3);
    % LIDAR_intensity              = nan(N_scanLines*max_lines,1);
    % LIDAR_scanLineAndRingID      = nan(N_scanLines*max_lines,2);

    fillData = nan(totalDataLength,(1+3+3+3+1+2));

    empty_XYZ_nanLine = [nan nan nan];

    % Loop through the scans we want to keep
    for ith_Scan = 1:N_scanLines
        if 0==mod(ith_Scan,100) && 1==flag_do_plots
            fprintf('Loading scan: %.0d of %.0d\n',ith_Scan,N_scanLines);
        end
        LIDAR_scan = LiDAR_Scan_ENU_Entire_Loop{ith_Scan};
        LIDAR_XYZ  = LIDAR_scan(:,1:3);
        LIDAR_intensity = LIDAR_scan(:,4);
        LIDAR_ringID = round(LIDAR_scan(:,5));

        indicies_to_keep = [];
        for ith_ring = 1:length(ringsRange)
            this_ring_indicies = find(LIDAR_ringID==ringsRange(ith_ring));
            indicies_to_keep = [indicies_to_keep; this_ring_indicies]; %#ok<AGROW>
        end
        Nindicies = length(indicies_to_keep);

        this_LIDAR_XYZ = [LIDAR_XYZ(indicies_to_keep,:); empty_XYZ_nanLine];
        this_LIDAR_intensity = [LIDAR_intensity(indicies_to_keep,:); nan];
        this_scanLine_rings = [ith_Scan*ones(Nindicies,1) LIDAR_ringID(indicies_to_keep,1); nan nan];

        % Save data        
        Npoints_added = length(this_LIDAR_XYZ(:,1));
        indicies_to_fill = (1:Npoints_added)' + max_lines*(ith_Scan-1);
        this_scanInfo = ith_Scan*ones(Npoints_added,1);

        % make sure we are not overwriting any data
        temp = fillData(indicies_to_fill,1);
        if any(~isnan(temp))
            error('Stop here');
        end

        this_vehiclePose = ones(Npoints_added,1)*VehiclePose(ith_Scan,1:3);
        this_orthoVectors =  ones(Npoints_added,1)*[unit_ortho_vehicle_vectors_XY(ith_Scan,:) 0];

        % LIDAR_ENU(indicies_to_fill,:) = thisLIDAR_XYZ_to_add;
        % LIDAR_intensity(indicies_to_fill,:)  = thisLIDAR_intensity_to_add;
        % LIDAR_scanLineAndRingID(indicies_to_fill,:) = thisScanLine_rings_to_add;
        % 
        % VehiclePose_ENU(indicies_to_fill,:) = ones(Npoints_added,1)*VehiclePose(ith_Scan,1:3);
        % VehiclePose_UnitOrthoVectors(indicies_to_fill,:) = ones(Npoints_added,1)*[unit_ortho_vehicle_vectors_XY(ith_Scan,:) 0];
        % 
        % ScanLineNumbers(indicies_to_fill,:) = ith_Scan;
        fillData(indicies_to_fill,:) = [this_scanInfo this_vehiclePose this_orthoVectors this_LIDAR_XYZ this_LIDAR_intensity this_scanLine_rings];
    end
    allData = fillData;
    alltotalDataLength = totalDataLength;

end


%% Get the LiDAR data based on ring and scanLine
startIndex = find(allData(:,1)==scanLineNumber_start,1,'first');
endIndex   = find(allData(:,1)==scanLineNumber_end,1,'last');
indicies_to_pull = (startIndex:endIndex)';
indices_to_pull = indicies_to_pull(~isnan(allData(indicies_to_pull,1)));

VehiclePose_ENU              = allData(indices_to_pull,2:4);
VehiclePose_UnitOrthoVectors = allData(indices_to_pull,5:7);
LIDAR_ENU                    = allData(indices_to_pull,8:10);
LIDAR_intensity              = allData(indices_to_pull,11);
LIDAR_scanLineAndRingID      = allData(indices_to_pull,12:13);



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

   fcn_findEdge_plotVehicleXY(VehiclePose_ENU,fig_num);   

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


%% fcn_INTERNAL_calcUnitVector
function unit_vectors = fcn_INTERNAL_calcUnitVector(input_vectors)
vector_length = sum(input_vectors.^2,2).^0.5;
unit_vectors = input_vectors./vector_length;
end
