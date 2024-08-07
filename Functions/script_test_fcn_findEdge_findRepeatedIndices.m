% script_test_fcn_findEdge_findRepeatedIndices
% Exercises the function: fcn_findEdge_findRepeatedIndices
% Revision history:
% 2024_6_17
% Jiabao Zhao wrote the code
% 2024_06_19 - Aneesh Batchu
% -- Added a case when "length_indices_array" is more than max(indices) of
% array
% -- Changed the name of the script

%% Test 1
fig_num = 1;
indices = [1, 1, 1, 2, 4, 4]';

length_indices_array = 4; 
repeated_indices = fcn_findEdge_findRepeatedIndices(indices,length_indices_array,fig_num);

assert(isequal(repeated_indices, [3 1 0 2]'))
assert(isequal(length(repeated_indices),4))

%% Test 2
fig_num = 1;
indices = [1, 1, 1, 2, 4, 4]';

length_indices_array = 5; 
repeated_indices = fcn_findEdge_findRepeatedIndices(indices,length_indices_array,fig_num);

assert(isequal(repeated_indices, [3 1 0 2 0]'))
assert(isequal(length(repeated_indices),5))

%%  Test 3 
fig_num = [];
num_points = 50;
range = [1,10];
indices = randi(range,num_points,1);
repeated_indices = fcn_findEdge_findRepeatedIndices(indices,length_indices_array,fig_num);
assert(isequal(length(repeated_indices),10));

