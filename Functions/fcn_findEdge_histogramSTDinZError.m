function fcn_findEdge_histogramSTDinZError(standard_deviation_in_z,N_bins_stdz,std_in_z_driven_path,N_bins_std_drivenpath,mean_std_in_z_driven_path,std_threshold,varargin)
%% fcn_findEdge_histogramSTDinZError
%
% FORMAT:
%   fcn_findEdge_histogramSTDinZError(standard_deviation_in_z,N_bins_stdz,std_in_z_driven_path,N_bins_std_drivenpath,mean_std_in_z_driven_path,std_threshold,(fig_num))
%
% INPUTS:
%
%   standard_deviation_in_z
%   N_bins_stdz
%   std_in_z_driven_path
%   N_bins_std_drivenpath
%   mean_std_in_z_driven_path
%   std_threshold
%
%
% (OPTIONAL INPUTS):
%   
%   fig_num
%
% OUTPUTS: 
%
%
% DEPENDENCIES:
%
% EXAMPLES:
%   See the script:
%   
%   script_test_fcn_findEdge_histogramSTDinZError
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
if (nargin==7 && isequal(varargin{end},-1))
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
        narginchk(6,7);
    end
end 

% Does user want to specify fig_num?
flag_do_plots = 0;
if (0==flag_max_speed) &&  (7<=nargin)
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

x_coord = mean_std_in_z_driven_path + 100*std(std_in_z_driven_path); 

current_text = sprintf('std threshold = %.4f',std_threshold);


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
    figure(fig_num);
    hold on
    grid on
    xlabel('Standard deviation in z error after plane fit')
    ylabel('Frequency')
    title('Statistic 3: Determining suitable standard deviation in z')
    % histogram(standard_deviation_in_z)
    histogram(standard_deviation_in_z,N_bins_stdz)
    histogram(std_in_z_driven_path,N_bins_std_drivenpath)


    disp('mean of std_threshold of driven path')
    disp(mean_std_in_z_driven_path)
    disp('chosen std_threshold')
    disp(std_threshold)

    % std_threshold = 0.05;
    plot(std_threshold,0, 'k.', 'MarkerSize',18)
    text(x_coord, 80,current_text,'Color',[0 0 0],'HorizontalAlignment','center','FontSize', 12, 'FontWeight','bold');
end


end