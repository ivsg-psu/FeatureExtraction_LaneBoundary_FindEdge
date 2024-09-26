% fcn_Path_fillRandomTraversalsAboutTraversal

% Load the data

savefile = fullfile(pwd, 'Data', 'ExampleData_findOrthoError.mat');
load(savefile, 'computed_boundary_path', 'hand_labeled_path')

% Plot the inputs
figure(346464);
clf;
hold on;
grid on;
grid minor;
axis equal

plot(hand_labeled_path(:,1),hand_labeled_path(:,2),'b-');
plot(computed_boundary_path(:,1),computed_boundary_path(:,2),'r.');

fig_num_StPointConverstion = [];
% figure(fig_num_StPointConverstion);
% clf;

St_points = fcn_Path_convertXY2St(hand_labeled_path,computed_boundary_path,[],fig_num_StPointConverstion);

% clear paths_array
% % Fill in sample paths (as a starter)
% paths_array{1} = hand_labled_path;
% paths_array{2} = computed_boundary_path;

% Convert all to traversals
hand_labeled_traversal = fcn_Path_convertPathToTraversalStructure(hand_labeled_path);
computed_boundary_traversal = fcn_Path_convertPathToTraversalStructure(computed_boundary_path);

% % Convert paths to traversal structures
% for i_Path = 1:length(paths_array)
%     traversal = ...
%         fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
%     all_traversals.traversal{i_Path} = traversal;
% end


% % Plot the results? (Note: they are plotted below as well)
% if 1==1
%     fig_num = 13;
%     fcn_Path_plotTraversalsXY(all_traversals,fig_num);
% end


reference_station_points = computed_boundary_traversal.Station;
reference_traversal = hand_labeled_traversal;
flag_rounding_type = 3; % Use average of projections at end points
search_radius = 10;
fig_num = 1;
clear allTraversals
allTraversals = struct;
allTraversals.traversal = cell(1,1);
allTraversals.traversal{1} = computed_boundary_traversal;
[closestXs, closestYs, closestDistances] = fcn_Path_findOrthoScatterFromTraversalToTraversals(reference_station_points, reference_traversal, allTraversals, flag_rounding_type,search_radius,fig_num); 
   

% % Peform histogram
% Nrows = length(closestDistances(:,1));
% Ncols = length(closestDistances(1,:));
% allDistances = reshape(closestDistances,[Nrows*Ncols 1]);
fig_num_colorized_error_LLA = 345; 
figure(fig_num_colorized_error_LLA); 
clf; 
XYI_data = []
fcn_plotRoad_plotXYI()

figure(1331);
histogram(closestDistances(:,2),30);
title('Initial Run: Histogram of all orthogonal distance projections');

OrthoDistances = closestDistances(~isnan(closestDistances(:,2)),2); 

sigma_std = std(closestDistances(~isnan(closestDistances(:,2)),2)); 
mu_mean = mean(closestDistances(~isnan(closestDistances(:,2)),2)); 


text(mu_mean - 3.8*sigma_std, 50, ['std = ' num2str(sigma_std)], 'FontSize', 12, 'Color', 'k');
text(mu_mean - 3.8*sigma_std, 45, ['mean = ' num2str(mu_mean)], 'FontSize', 12, 'Color', 'k');
