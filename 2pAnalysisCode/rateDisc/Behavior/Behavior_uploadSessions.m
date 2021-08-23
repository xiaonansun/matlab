function Behavior_uploadSessions(sPath, tPath)
% code to upload behavior data files from local setup computers to the
% server. This is to ensure that files are not missing when the network
% connection is disrupted. 
% sPath should point to Bpods data folder. Code checks SessionData in all
% animals and copys files that are not found in the according tPath.

if sPath(end) ~= filesep
    sPath(end + 1) = filesep;
end

if tPath(end) ~= filesep
    tPath(end + 1) = filesep;
end

%% loop through subjects, paradigms and files
subjects = dir(sPath);

for iSubs = 1 : length(subjects)
    if ~(strcmp(subjects(iSubs).name, '.') || strcmp(subjects(iSubs).name, '..'))

        paradigms = dir([sPath subjects(iSubs).name]);
        for iParams = 1 : length(paradigms)
            if ~(strcmp(paradigms(iParams).name, '.') || strcmp(paradigms(iParams).name, '..'))
                
                files = dir([sPath subjects(iSubs).name filesep paradigms(iParams).name filesep 'Session Data' filesep]);
                for iFiles = 1 : length(files)
                    if ~(strcmp(files(iFiles).name, '.') || strcmp(files(iFiles).name, '..'))
                        
                        % check if target file exists and copy if not
                        sFile = [sPath subjects(iSubs).name filesep paradigms(iParams).name ... 
                            filesep 'Session Data' filesep files(iFiles).name];
                        
                        tFile = [tPath subjects(iSubs).name filesep paradigms(iParams).name ... 
                            filesep 'Session Data' filesep files(iFiles).name];
                        
                        if ~exist(tFile, 'file')
                            if ~exist([tPath subjects(iSubs).name filesep paradigms(iParams).name filesep 'Session Data' filesep] , 'dir')
                                mkdir([tPath subjects(iSubs).name filesep paradigms(iParams).name filesep 'Session Data' filesep]);
                            end
                            
                            copyfile(sFile, tFile); %copy missing file to server
                            
                        end
                    end
                end
            end
        end
    end
end
