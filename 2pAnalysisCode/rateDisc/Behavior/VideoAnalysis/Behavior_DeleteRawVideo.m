function Behavior_DeleteRawVideo(fPath)
% code to convert behavioral videos from timestamp (.mat) and compressed movie (.mj2) files.

if fPath(end) ~= filesep
    fPath = [fPath filesep];
end

%% folders that contain animal video data
animals = dir(fPath);
animals = animals([animals.isdir] & ~strncmpi('.', {animals.name}, 1));

%% loop through animal data
for iAnimals = 1:length({animals.name})
    
    aPath = [animals(iAnimals).name filesep];
    paradigms = dir([fPath aPath]);
    paradigms = paradigms([paradigms.isdir] & ~strncmpi('.', {paradigms.name}, 1));
    
    %% loop through paradigms
    for iParadigms = 1:length({paradigms.name})
        
        pPath = [paradigms(iParadigms).name filesep];
        sessions = dir([fPath aPath  pPath]);
        sessions = sessions([sessions.isdir] & ~strncmpi('.', {sessions.name}, 1));
        
        %% loop through sessions
        for iSessions = 1:length({sessions.name})
            
            sPath = [sessions(iSessions).name filesep];
            experiments = dir([fPath aPath  pPath sPath]);
            experiments = experiments([experiments.isdir] & ~strncmpi('.', {experiments.name}, 1));
            
            %% loop through experiments
            for iExperiments = 1:length({experiments.name})
                ePath = [experiments(iExperiments).name filesep];
                files = dir([fPath aPath pPath sPath ePath '*' ePath(1:end-1) '*frameTimes*.mat']);
                if ~isempty(files)
                    temp = cat(1,files.name);
                    ind = strfind(temp(1,:),'_');ind = ind(end);
                    nrCams = length(unique(str2num(temp(:,ind+1))));
                    
                    %% loop through available webcams
                    for iCams = 1:nrCams
                        
                        cPath = [fPath aPath  pPath sPath ePath]; %path to current experiment
                        movieFiles = dir([fPath aPath pPath sPath ePath '*' ePath(1:end-1) '*Video*_' int2str(iCams) '.mj2']); %all files for current cam based on integer before .mat
                        
                        svdFile = dir([cPath 'SVD_Cam ' int2str(iCams) '*']);
                        svdCheck = dir([cPath 'SVD_Complete.mat']);
                        
                        if size(svdFile,1) > 0 && size(svdCheck,1) > 0
                            for iFiles = 1:length({movieFiles.name})
                                delete([cPath movieFiles(iFiles).name]);
                            end
                            disp(['Deleted data for cam' int2str(iCams) ' - ' experiments(iExperiments).name]);
                        end
                        
                    end
                end
            end
        end
    end
end
