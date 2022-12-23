function [summaryTable,error_log]= twoP_galvoFrameCount
%%

exps = twoP_getAcquisitionRecord;
colAnimal = exps(:,1);
colSession = exps(:,6);
imagingRootDir = S.dir.imagingRootDir;
imagingSubDir = S.dir.imagingSubDir;
error_log = cell(size(exps,1),1);

parfor i = 1:size(exps,1)
    %% i=143
    animal = colAnimal{i};
    session = colSession{i};
    data_path = fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'data.mat');
    sdu_Vc_savepath = fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'sduVc.mat');
    bhv_path = fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'cBhv.mat');
    try
    data = load(data_path); data = data.data;
    oAnimal{i}= animal;
    oSession{i} = session;
    oStimSampOrig{i} = data.stimSamplesOrig;

    catch ME
        disp(ME.message)
        error_log{i}= ME.message;
    end
    
end

idxMultBins = find(strcmp(cellfun(@class,oStimSampOrig,'UniformOutput',false),'cell'));
idxSingleBin = find(~strcmp(cellfun(@class,oStimSampOrig,'UniformOutput',false),'cell') & ~cell2mat(cellfun(@isempty,oStimSampOrig,'UniformOutput',false)));

idxSessions = idxSingleBin; % change this 
parfor i = 1:length(idxSessions)
%     i=5
    animal = colAnimal{idxSessions(i)};
    session = colSession{idxSessions(i)};
    data_path = fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'data.mat');
    data = load(data_path); data = data.data;
    spks_npy_path = fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir,'spks.npy');
    spks = readNPY(spks_npy_path);
    numOfBins(i) = length(data.stimFramesOrig);
    
    slowGalvoCh = 2;
%     binFiles = dir(fullfile(imagingRootDir,animal,'imaging',session,'*.bin'));
%     binFilename = {binFiles.name};
%     binFileDir = {binFiles.folder};
    npy = twoP_importSuite2pNPY(D.animal,D.session);
    if length(npy.bin_MScan_filepath) == 1
%         binFilePath = fullfile(binFileDir{1},binFilename{1});
        [volt, ~] = readMOMAnalog(npy.bin_MScan_filepath{1});
        slowGalvo = volt(slowGalvoCh, :);
        [frameStarts, incompleteFrames] = parseSlowGalvo(slowGalvo);
        numOfFrames = length(frameStarts);
        numOfIncompleteFrames = length(incompleteFrames);
    elseif length(npy.bin_MScan_filepath) > 1
        for j = 1:length(npy.bin_MScan_filepath)
            %% j = 1
%             binFilePath = fullfile(binFileDir{j},binFilename{j});
            [volt, ~] = readMOMAnalog(npy.bin_MScan_filepath{j});
            slowGalvo{j} = volt(slowGalvoCh, :);
            [frameStarts{j}, incompleteFrames{j}] = parseSlowGalvo(slowGalvo{j});
            numOfFrames(j) = length(frameStarts{j});
            numOfIncompleteFrames(j) = length(incompleteFrames{j});
            
        end
    end
    galvoFrames(i) = sum(numOfFrames);
    galvoFramesIncomplete(i) = sum(numOfIncompleteFrames);
    spksFrames(i) = length(spks);
    disp(['Total number of frames counted from Galvo: ' num2str(galvoFrames(i))]);
    disp(['Total number of incompleted frames counted from Galvo: ' num2str(galvoFramesIncomplete(i))]);
    disp(['Total number of frames counted from spks.npy: ' num2str(spksFrames(i))]);
    
end

diffFrames = galvoFrames - spksFrames;
format long;
summaryTable = table({colAnimal{idxSessions}}',{colSession{idxSessions}}',galvoFrames', galvoFramesIncomplete',spksFrames',diffFrames');
summaryTable.Properties.VariableNames={'Animal','Session','Galvo_Frames','Galvo_Frames_Incomplete','spks_Frames','diff_Frames (galvo_Frames - spks_Frames)'};

save_csv_filename = 'galvo_spks_frame_counts.csv';
save_csv_subdir = 'frame_counts';
save_csv_dir = fullfile(S.dir.imagingRootDir,save_csv_subdir);
if ~exist('save_csv_dir','dir'); mkdir(save_csv_dir); end
writetable(summaryTable,fullfile(save_csv_dir,save_csv_filename));


