
load("boundary_points_entire_test_track.mat"); 

fig_num = 408; 
figure(fig_num); 

marker_size = 30;
RGB_triplet = [0 1 1]; 
legend_option = 0;
legend_name = 'Right Boundary Points';
legend_position = [];
marker_type = []; 

plot_boundary_points_right = [boundary_points_test_track, zeros(length(boundary_points_test_track),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);


marker_size = 12;
RGB_triplet = [1 0 0]; 
legend_option = 0;
legend_name = 'Right Boundary Points';
legend_position = [];
marker_type = []; 

plot_boundary_points_right = [boundary_points_test_track, zeros(length(boundary_points_test_track),1)];
[~] = fcn_findEdge_plotPointsinLLA(plot_boundary_points_right,marker_size,RGB_triplet,marker_type,legend_option,legend_name,legend_position,[],[],[],fig_num);

