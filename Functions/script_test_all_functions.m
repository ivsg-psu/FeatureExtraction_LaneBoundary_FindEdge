%% script_test_all_functions.m
% This is a wrapper script to run all the test scripts in the 
% library for the purpose of evaluating every assertion test in these
% files.
%
% NOTE: to view the output file with formatting, use the "type" command.
% For example:
% type('script_test_fcn_geometry_all_stdout.txt')

clearvars; 
close all; 
clc;



directoryQuery = fullfile(pwd,'Functions','**','*.*');

all_files = dir(directoryQuery); %cat(2,'.',filesep,filesep,'script_test_fcn_*.m'));
N_files = length(all_files);
testing_times = nan(N_files,1);

diary(fullfile(pwd,'script_test_fcn_FindEdge_all_stdout.txt'));

flags_isFile  = zeros(N_files,1);
flags_isMfile = zeros(N_files,1);
flags_isMfileRepeated = zeros(N_files,1);
flags_isMfileTestingScript = zeros(N_files,1);
flags_isMfileTestedFunction = zeros(N_files,1);
flags_isMfileTestingScriptWithMatchingFunction = zeros(N_files,1);

% Check all the files to see which ones should be tested
for i_script = 1:N_files
    % Is this a file?
    if ~all_files(i_script).isdir
        flags_isFile(i_script,1) = 1;
    end
    
    file_name_extended = all_files(i_script).name;
    
    % Is this an m-file?
    if (1==flags_isFile(i_script,1) ) && length(file_name_extended)>2 && strcmp(file_name_extended(end-1:end),'.m')
        flags_isMfile(i_script,1) = 1;

        % Is this a repeated m-file?
        for jth_file = 1:N_files
            nameToTest = all_files(jth_file).name;
            if (i_script~=jth_file) && strcmp(file_name_extended,nameToTest)
                flags_isMfileRepeated(jth_file,1) = 1;
            end
        end
    end


    % Is this an m-file testing script?
    flag_file_was_found = 0;
    if (1==flags_isMfile(i_script,1) ) && length(file_name_extended)>19 && strcmp(file_name_extended(1:16),'script_test_fcn_')
        flags_isMfileTestingScript(i_script,1) = 1;

        % Mark the corresponding m-file as tested
        testMfileName = file_name_extended(13:end);
        fullPath = which(testMfileName);
        for jth_file = 1:N_files
            listedFullName = fullfile(all_files(jth_file).folder,all_files(jth_file).name);
            if strcmp(fullPath,listedFullName)
                flags_isMfileTestedFunction(jth_file,1) = 1;
                flag_file_was_found = 1;
            end
        end
    end

    % Is this an m-file testing script with matching function?
    if (1==flags_isMfileTestingScript(i_script,1) ) && (1==flag_file_was_found)
        flags_isMfileTestingScriptWithMatchingFunction(i_script,1) = 1;
    end

end

flags_isEitherTestScriptOrTestedFunction = flags_isMfileTestedFunction+flags_isMfileTestingScript;

%% Summarize results
fprintf(1,'\nSUMMARY OF FOUND FILES: \n');
indicies_filesToTest = find(1==flags_isMfileTestingScriptWithMatchingFunction);
if ~isempty(indicies_filesToTest)
    fcn_DebugTools_cprintf('*blue','The following scripts were found that will be tested:\n');
    for ith_file = 1:length(indicies_filesToTest)
        currentFileIndex = indicies_filesToTest(ith_file);
        fprintf(1,'\t%s\n',all_files(currentFileIndex).name)
    end    
end

% List missed files
indicies_missedFiles = find(1==(flags_isFile.*(0==flags_isMfile)));
if ~isempty(indicies_missedFiles)
    fcn_DebugTools_cprintf('*red','The following files were found, but do not seem to be matlab functions or scripts:\n');
    for ith_file = 1:length(indicies_missedFiles)
        currentFileIndex = indicies_missedFiles(ith_file);
        fcn_DebugTools_cprintf('*red','\t%s\n',all_files(currentFileIndex).name)
    end    
end

% List missed mfiles
indicies_missedMfiles = find(1==(flags_isMfile.*(0==flags_isEitherTestScriptOrTestedFunction)));
if ~isempty(indicies_missedMfiles)
    fcn_DebugTools_cprintf('*red','The following m-files were found, but do not have a test script:\n');
    for ith_file = 1:length(indicies_missedMfiles)
        currentFileIndex = indicies_missedMfiles(ith_file);
        fcn_DebugTools_cprintf('*red','\t%s\n',all_files(currentFileIndex).name)
    end    
end

% List missed scripts
indicies_missedScripts = find(1==(flags_isMfileTestingScript.*(0==flags_isMfileTestingScriptWithMatchingFunction)));
if ~isempty(indicies_missedScripts)
    fcn_DebugTools_cprintf('*red','The following test scripts were found, but do not have a matching function:\n');
    for ith_file = 1:length(indicies_missedScripts)
        currentFileIndex = indicies_missedScripts(ith_file);
        fcn_DebugTools_cprintf('*red','\t%s\n',all_files(currentFileIndex).name)
    end    
end

% List repeats
indicies_repeatedFiles = find(1==flags_isMfileRepeated);
if ~isempty(indicies_repeatedFiles)
    fcn_DebugTools_cprintf('*red','The following files seem to be repeated:\n');
    for ith_file = 1:length(indicies_repeatedFiles)
        currentFileIndex = indicies_repeatedFiles(ith_file);
        fcn_DebugTools_cprintf('*red','\t%s\n',all_files(currentFileIndex).name)
    end    
end

%% Test the good scripts
NtestScripts = length(indicies_filesToTest);
fprintf(1,'\nSTARTING TESTS:\n');
allResults = cell(NtestScripts,1);
for ith_testScript = 1:NtestScripts
    currentFileIndex = indicies_filesToTest(ith_testScript);
    file_name_extended = all_files(currentFileIndex).name;
    file_name = erase(file_name_extended,'.m');
    file_name_trunc = erase(file_name,'script_');
    fcn_DebugTools_cprintf('*blue','%s\n','   ');
    fcn_DebugTools_cprintf('*blue','Testing script: %.0d of %.0d, %s\n\n',ith_testScript,NtestScripts,file_name_trunc);

    % Start the test
    tstart = tic;
    suite = testsuite(file_name);
    allResults{ith_testScript,1} = run(suite);
    telapsed = toc(tstart);
    testing_times(ith_testScript) = telapsed;
end
diary off

close all;
figure(1)
plot(testing_times);
grid on;
xlabel('Script test number');
ylabel('Elapsed time to test (sec)');



