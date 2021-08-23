function Behavior_ConvertVideo(cPath)
% code to convert behavioral videos from uint8 binary files to timestamps
% of indivudal frames (.mat) and compressed movie (.mj2) files.

if cPath(end) ~= filesep
    cPath = [cPath filesep];
end

%% folders that contain animal video data
animals = dir(cPath);
animals = animals([animals.isdir] & ~strncmpi('.', {animals.name}, 1));

%% loop through animal data
for iAnimals = 1:length({animals.name})
    
    aPath = [animals(iAnimals).name filesep];
    paradigms = dir([cPath aPath]);
    paradigms = paradigms([paradigms.isdir] & ~strncmpi('.', {paradigms.name}, 1));
    
    %% loop through paradigms
    for iParadigms = 1:length({paradigms.name})
            
        pPath = [paradigms(iParadigms).name filesep];
        sessions = dir([cPath aPath  pPath]);
        sessions = sessions([sessions.isdir] & ~strncmpi('.', {sessions.name}, 1));
        
        %% loop through sessions
        for iSessions = 1:length({sessions.name})
            
            sPath = [sessions(iSessions).name filesep];
            experiments = dir([cPath aPath  pPath sPath]);
            experiments = experiments([experiments.isdir] & ~strncmpi('.', {experiments.name}, 1));
                    
            %% loop through experiments
            for iExperiments = 1:length({experiments.name})
                
                ePath = [experiments(iExperiments).name filesep];
                files = dir([cPath aPath pPath sPath ePath '*' ePath(1:end-1) '*.dat']);
                files = files(~strncmpi('.', {files.name}, 1));
                
                %% loop through files
                for iFiles = 1:length({files.name})
                
                    cFile = [cPath aPath  pPath sPath ePath files(iFiles).name];
                    [a,b] = fileparts(cFile);
                    
                    %% load video data and timestamps and save as compressed file
                    disp(['Converting file ' b]);tic
                    fID = fopen(cFile);
                    hSize = fread(fID,1,'double'); %header size
                    frameTimes = fread(fID,hSize,'double'); %Metadata. Default is 1:hSize = Absolute timestamps for each frame
                    dSize = fread(fID,1,'double'); %numer of data array dimensions
                    dSize = fread(fID,dSize,'double'); %data size
                    data = fread(fID,[prod(dSize),1],'uint8=>uint8'); %get data. Last 4 header values should contain the size of the data array.
                    data = reshape(data,dSize'); %reshape data into matrix
                    fclose(fID); clear fID
                    
                    if isunix %check for behavioral data
                        exist([a1 '\' b1 '.mat'],'file');
                        copyfile([a1 '\' b1 '.mat']
                    end
                    
                    if isunix && isdir('/mnt/managed/') %setup PC with mounted network path
                        [a1, b1] = fileparts(a);
                        
                        a = strrep(a,cPath,'/mnt/managed/'); %write to network instead of local
                        if ~isdir(a); mkdir(a);end
                        
                        if exist([a1 '\' b1 '.mat'],'file') == 2 %check for behavioral data
                            copyfile([a1 '\' b1 '.mat'],[a '.mat']); %copy to server
                        end
                        
                    else
                        disp('No network drive detected. Writing to source folder instead.');
                    end
                    
                    save([a filesep strrep(b,'Video','frameTimes')],'frameTimes');
                    v = VideoWriter([a filesep b],'Archival');
                    open(v);writeVideo(v,data);close(v);
                    delete(cFile);
                    disp('Conversion complete');toc
                    disp('-–-–-–-–-–-–-–-–-–-–-–-–-–-–-–-–-–-–-–-–-–-–-–-–-–-–');
                    
                end
            end
        end
    end
end