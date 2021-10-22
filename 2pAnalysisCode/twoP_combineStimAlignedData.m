function D= twoP_combineStimAlignedData(data)
% 2021-10-04 This script combines stimulus-aligned 2P data primary to deal
% with MScan computer crashes. When the MScan PC crashes, imaging is
% terminated while the behavior continues, hence one or more trials will be
% skipped in the 2P data. This results in missing trial codes from the
% analog input. The script twoP_alignDetectionTask.m has been modified to
% combine multiple 2P sub-sessions into a single data struct. This data
% struct separates individual sub-sessions into cells. This script
% concatenates data.trialNumbers, data.stimFrameOrig, data.neural, data.DS,
% and data.analog.

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
    D.neuralTimes = data.neuralTimes{1};
    D.neuralTimesMs = data.neuralTimesMs{1};
    D.DS = cat(3,data.DS{1},data.DS{2});
    D.trialStimSample = data.trialStimSample{1};
    D.stimSamplesOrig = horzcat(data.stimSamplesOrig{1},data.stimSamplesOrig{2});
    D.analog = cat(3,data.analog{1},data.analog{2});
    D.analogTimes = data.analogTimes{1};
    D.idx_redcell = data.idx_redcell;
    D.idx_notredcell = data.idx_notredcell;
    D.suite2pDir = data.suite2pDir;
else
    D = data;
end