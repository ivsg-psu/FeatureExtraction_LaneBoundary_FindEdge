%% SCRIPT OF JUST RUNNING THE FUNCTION (NEEDS WORK)
fig_1=1;
fig_2=2;
fig_3=3;

[input_points,original_mapped_gridIndices_cell,total_mapped_grids,total_points_in_mapped_grids,standard_deviation_in_z,gridlines_mapped_grids,...
    driven_path_grid_indices_in_current_mapped_grids,std_in_z_driven_path,std_in_z_other_mapped_grids,mean_std_in_z_driven_path,mean_std_in_z_not_driven_path,max_std_in_z_not_driven_path] = fcn_findEdge_determineSTDInZError(LiDAR_allPoints,gridIndices_cell_array,original_qualified_grids,gridCenters_qualified_grids,gridCenters_driven_path,...
    current_qualified_grids,grid_AABBs,grid_size,gridIndices,current_grid_numbers_of_driven_path,fig_1,fig_2,fig_3);