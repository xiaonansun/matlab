function compressTifStacks(rootDir)
% This function (1) looks for .TIF files recursively within

% rootDir = 'D:\raw';

list=dir([rootDir filesep '**/*.TIF']);

num_of_files = length(list);
% num_of_files = 10;

tic

parfor i = 1:num_of_files
    [status,result]=system(['7z a ' list(i).folder filesep list(i).name '.bz2 ' list(i).folder filesep list(i).name]);
    disp(['Compressing ' list(i).folder filesep list(i).name]);
    delete([list(i).folder filesep list(i).name])
end

toc
