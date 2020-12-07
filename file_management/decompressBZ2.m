rootDir = 'D:\raw\CSP30';

list=dir([rootDir filesep '**/*.bz2']);

num_of_files = length(list);
% num_of_files = 5;

tic

parfor i = 1:num_of_files
    [status,result]=system(['7z e ' list(i).folder filesep list(i).name ' -o' list(i).folder]);
    disp(['Decompressing ' list(i).folder filesep list(i).name]);
    delete([list(i).folder filesep list(i).name]);
end

disp(['Processing completed in ' datestr(seconds(toc),'HH:MM:SS') '.']);