%script_test_fcn_findEdge_prepGridCentersForBoundaryDetection
% Exercises the function: fcn_findEdge_prepGridCentersForBoundaryDetection
%
% Revision History:
% 2024_08_015 -Jiabao Zhao
% -- wrote the script orginally 

%% Test 1 
gridCenters_qualified_grids = [1 2; 3 4]; % Two qualified grid centers
gridCenters_unqualified_grids = [1 2; 5 6]; % Two unqualified grid centers (one overlaps with qualified)

[X, Y, Z] = fcn_findEdge_prepGridCentersForBoundaryDetection...
    (gridCenters_qualified_grids, gridCenters_unqualified_grids);
assert(isequal(X(1,:), [1 1 1]))
assert(isequal(Y(1,:), [2 4 6]))
assert(isequal(size(Z),size(X)))


%% Test 2  
% Define very simple example grid centers for qualified and unqualified grids
gridCenters_qualified_grids = [1 2]; % One qualified grid center
gridCenters_unqualified_grids = [1 2; 3 4]; % Two unqualified grid centers (one overlaps with the qualified one)

[X, Y, Z] = fcn_findEdge_prepGridCentersForBoundaryDetection...
    (gridCenters_qualified_grids, gridCenters_unqualified_grids);
assert(isequal(X(1,:), [1 3]))
assert(isequal(Y(1,:), [2 2]))
assert(isequal(size(Z),size(X)))

%% Test 3  more complex example 
gridCenters_qualified_grids = [1 2; 3 4; 5 6; 7 8]; % Four qualified grid centers
gridCenters_unqualified_grids = [1 2; 3 5; 5 6; 9 10]; % Four unqualified grid centers, with some overlaps with qualified

[X, Y, Z] = fcn_findEdge_prepGridCentersForBoundaryDetection...
    (gridCenters_qualified_grids, gridCenters_unqualified_grids);
assert(isequal(size(X), size(Y)))
assert(isequal(size(X), size(Z)))