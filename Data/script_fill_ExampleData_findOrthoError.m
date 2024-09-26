


%% Hand-labled Points
load /Users/aneeshbatchu/Documents/IVSG/Aneesh/Repeatability_VeloLiDARdata/2024_08_05/CCW/Figures/HandLabeled_BDpts/hand_labeled_path.mat

%% Computed Boundary Points
load /Users/aneeshbatchu/Documents/IVSG/Aneesh/Repeatability_VeloLiDARdata/2024_08_05/CCW/Figures/ClosestDistances/Run6/test_track_organized_points_ENU.mat

%% Ortho Distances
load /Users/aneeshbatchu/Documents/IVSG/Aneesh/Repeatability_VeloLiDARdata/2024_08_05/CCW/Figures/ClosestDistances/Run6/ortho_distances.mat

computed_boundary_path = path; 
clear path

savefile = fullfile(pwd, 'Data', 'ExampleData_findOrthoError_6.mat');
save(savefile, 'computed_boundary_path', 'hand_labeled_path', 'OrthoDistances')
