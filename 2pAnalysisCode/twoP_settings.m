function S = twoP_settings

S = struct;

if convertCharsToStrings(getenv('COMPUTERNAME')) == convertCharsToStrings('MANHASSET')
    S.dir.imagingRootDir = 'G:\2PData';
elseif convertCharsToStrings(getenv('COMPUTERNAME')) == convertCharsToStrings('SUNHP')
    %     imagingRootDir = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\suite2p';
    S.dir.imagingRootDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy';
    S.dir.codeRootDir = 'C:\Users\Xiaonan Richard Sun\Dropbox\Users\Richard\matlab';
    S.dir.bhvRootDir = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon';
    S.isUnix = false;
elseif convertCharsToStrings(computer) == convertCharsToStrings('GLNXA64') || isunix == 1
    S.dir.imagingRootDir = '/grid/churchland/data/data/richard_s2p_npy';
    S.dir.codeRootDir = '/grid/churchland/home/xisun/matlab';
    S.dir.bhvRootDir = '/grid/churchland/data/data/Behavior_Simon';
    S.isUnix = true;
end

S.dir.imagingSubDir = fullfile('suite2p', 'plane0');
S.dir.bhvSubDir = fullfile('SpatialDisc', 'Session Data');

S.segIdx = [1 0.75 1.25 0.5 1];
S.cellTypes = {'CSP';'Plex';'Fez'};
S.expertise = {'Naive','Trained','Expert'};
S.nShuffle = 20;
S.msPerFrame = 32.3638;
S.trialStimFrame = 93;

S.opts.preStim = S.trialStimFrame*S.msPerFrame/1000; % Duration of the data (in seconds) before the stimulus occurs
S.opts.frameRate = 1000/S.msPerFrame; % Frame rate of imaging
S.sRate = S.opts.frameRate;

S.segFrames = cumsum(floor(S.segIdx * S.sRate)); %max nr of frames per segment

S.epoches(1).idx = [S.segFrames(1) S.segFrames(2)];
S.epoches(1).name = 'Init';
S.epoches(2).idx = [S.segFrames(2)+1 S.segFrames(3)];
S.epoches(2).name = 'Stim';
S.epoches(3).idx = [S.segFrames(2)+1 floor(mean(S.epoches(2).idx))];
S.epoches(3).name = 'StimEarly';
S.epoches(4).idx = [ceil(mean(S.epoches(2).idx)) S.segFrames(3)];
S.epoches(4).name = 'StimLate';
S.epoches(5).idx = [S.segFrames(3)+1 S.segFrames(4)];
S.epoches(5).name = 'Delay';
S.epoches(6).idx = [S.segFrames(4)+1 S.segFrames(5)];
S.epoches(6).name = 'Response';
S.epoches(7).idx = [S.segFrames(4)+1 floor(mean(S.epoches(6).idx))];
S.epoches(7).name = 'ResponseEarly';
S.epoches(8).idx = [ceil(mean(S.epoches(6).idx)) S.segFrames(5)];
S.epoches(8).name = 'ResponseLate';
S.allEpoches.idx = [S.epoches(1).idx; S.epoches(3).idx; S.epoches(4).idx;...
    S.epoches(5).idx; S.epoches(7).idx; S.epoches(8).idx];
S.allEpoches.name = {S.epoches(1).name; S.epoches(3).name; S.epoches(4).name;...
    S.epoches(5).name; S.epoches(7).name; S.epoches(8).name};