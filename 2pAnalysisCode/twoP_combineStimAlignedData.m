function D= twoP_combineStimAlignedData(data)
% 2021-10-04 This script combines stimulus-aligned 2P data  to
% correct for MScan crashes. When the MScan PC crashes, image acquisition is
% terminated while the Bpod continues to run. As a result, the 2P data will
% skip one or more trials. This results in missing trial codes from the
% analog input. The script twoP_alignDetectionTask.m has been modified to
% combine multiple 2P sub-sessions into a single data struct. This data
% struct separates individual sub-sessions into cells. This script
% concatenates data.trialNumbers, data.stimFrameOrig, data.neural, data.DS,
% and data.analog.
S = twoP_settings;

if iscell(data.neural)
    D.msPerFrame = data.msPerFrame;
    D.analogIdentities = data.analogIdentities;
    D.animal = data.animal;
    D.session = data.session;
    D.analogFreq = data.analogFreq{1};
    D.trialNumbers = horzcat(data.trialNumbers{1},data.trialNumbers{2});
    D.trialStimFrame = data.trialStimFrame{1};
    D.stimFramesOrig = horzcat(data.stimFramesOrig{1},data.stimFramesOrig{2});
    D.neural = cat(3,data.neural{1},data.neural{2});
    D.sdu = cat(3,data.sdu{1},data.sdu{2});
    D.neuralTimes = data.neuralTimes{1};
    D.neuralTimesMs = data.neuralTimesMs{1};
    D.DS = cat(3,data.DS{1},data.DS{2});
    D.trialStimSample = data.trialStimSample{1};
    D.stimSamplesOrig = horzcat(data.stimSamplesOrig{1},data.stimSamplesOrig{2});
    D.analog = cat(3,data.analog{1},data.analog{2});
    D.analogTimes = data.analogTimes{1};
    D.idx_redcell = data.idx_redcell;
    D.idx_notredcell = data.idx_notredcell;
    
    slowGalvoCh = 2;
    binFiles = dir(fullfile(S.dir.imagingRootDir,D.animal,'imaging',D.session,'*.bin'));
    binFilename = {binFiles.name}; binFileDir = {binFiles.folder};
    npyFilepath = fullfile(S.dir.imagingRootDir,D.animal,'imaging',D.session,S.dir.imagingSubDir,'spks.npy');
    spks = readNPY(npyFilepath);
    
%     clear slowGalvo frameStarts incompleteFrames numOfGalvoFrames numOfIncompleteFrames
%     for j = 1:length(binFilename)
%         %% j = 1
%         [volt, ~] = readMOMAnalog(fullfile(binFileDir{j},binFilename{j}));
%         slowGalvo{j} = volt(slowGalvoCh, :);
%         [frameStarts{j}, incompleteFrames{j}] = parseSlowGalvo(slowGalvo{j});
%         numOfGalvoFrames(j) = length(frameStarts{j});
%         numOfIncompleteFrames(j) = length(incompleteFrames{j});
%     end
    
    if isfield(D,'sdu')
        D.sdu = cat(3,data.sdu{1},data.sdu{2});
    else
        disp('No z-scored trials were found, skipping.')
    end
    
    if exist('data.suite2pDir','var')
        D.suite2pDir = data.suite2pDir;
    else
        D.suite2pDir = [];
    end
    
else
    D = data;
end