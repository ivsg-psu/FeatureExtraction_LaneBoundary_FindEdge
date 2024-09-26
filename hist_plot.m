% fcn_Path_fillRandomTraversalsAboutTraversal

% Load the data

savefile = fullfile(pwd, 'Data', 'ExampleData_findOrthoError.mat');
load(savefile, 'computed_boundary_path', 'hand_labeled_path')

% NOTE: the hand_labeled_path is backwards (CW). Need to fix it
hand_labeled_path = flipud(hand_labeled_path);

% Plot the inputs?
if 1==0
    figure(346464);
    clf;
    hold on;
    grid on;
    grid minor;
    axis equal

    plot(hand_labeled_path(:,1),hand_labeled_path(:,2),'b-');
    Npoints = length(hand_labeled_path(:,1));
    for ith_point = 2:Npoints
        connecting_vectors = hand_labeled_path(ith_point,:) - hand_labeled_path(ith_point-1,:);
        quiver(hand_labeled_path(ith_point-1,1),hand_labeled_path(ith_point-1,2),...
            connecting_vectors(1,1),connecting_vectors(1,2),...
            0,'Color',[0.5 0.5 0.5]);
    end
    plot(computed_boundary_path(:,1),computed_boundary_path(:,2),'r.');
end

fig_num_StPointConverstion = []; %1234;
% figure(fig_num_StPointConverstion);
% clf;

St_points = fcn_Path_convertXY2St(hand_labeled_path,computed_boundary_path,[],fig_num_StPointConverstion);
transverse_error = real(St_points(:,2));

max_transverse_error = 2; % Set by user - this determines the color scale
normalized_abs_transverse_error = min(2,abs(transverse_error))/max_transverse_error;


%%
fig_num_colorized_error_XYI = 345; 
figure(fig_num_colorized_error_XYI); 
clf; 
XYI_data = [computed_boundary_path normalized_abs_transverse_error];


% Call the plotting function
clear plotFormat
plotFormat.LineStyle = 'none';
plotFormat.LineWidth = 5;
plotFormat.Marker = '.';
plotFormat.MarkerSize = 10;

colorMapMatrixOrString = colormap('turbo');
Ncolors = 10;
reducedColorMap = fcn_plotRoad_reduceColorMap(colorMapMatrixOrString, Ncolors, -1);

fcn_plotRoad_plotXYI(XYI_data, (plotFormat),  (reducedColorMap), (fig_num_colorized_error_XYI));
title('Boundary error, ENU');

%% DO LLA
% Convert ENU data to LLA
% Define reference location at test track
reference_latitude = 40.86368573;
reference_longitude = -77.83592832;
reference_altitude = 344.189;

% Location for Test Track base station
setenv('MATLABFLAG_PLOTROAD_REFERENCE_LATITUDE','40.86368573');
setenv('MATLABFLAG_PLOTROAD_REFERENCE_LONGITUDE','-77.83592832');
setenv('MATLABFLAG_PLOTROAD_REFERENCE_ALTITUDE','344.189');

gps_object = GPS(reference_latitude,reference_longitude,reference_altitude); % Initiate the class object for GPS
% ENU_full_data = gps_object.WGSLLA2ENU(LLA_full_data(:,1), LLA_full_data(:,2), LLA_full_data(:,3));
ENU_full_data = [XYI_data(:,1), XYI_data(:,2), XYI_data(:,2)*0];
LLA_full_data  =  gps_object.ENU2WGSLLA(ENU_full_data);

fig_num_colorized_error_LLI = 678; 
figure(fig_num_colorized_error_LLI); 
clf; 
LLI_data = [LLA_full_data(:,1:2) normalized_abs_transverse_error];

[h_plot, indiciesInEachPlot]  = fcn_plotRoad_plotLLI(LLI_data, (plotFormat),  (reducedColorMap), (fig_num_colorized_error_LLI));
title('Boundary error, ENU');

% Force the plot to fit
geolimits('auto');

%%

figure(1331);
histogram(closestDistances(:,2),30);
title('Initial Run: Histogram of all orthogonal distance projections');

OrthoDistances = closestDistances(~isnan(closestDistances(:,2)),2); 

sigma_std = std(closestDistances(~isnan(closestDistances(:,2)),2)); 
mu_mean = mean(closestDistances(~isnan(closestDistances(:,2)),2)); 


text(mu_mean - 3.8*sigma_std, 50, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 3.8*sigma_std, 45, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
