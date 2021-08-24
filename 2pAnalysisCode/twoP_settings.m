function S = twoP_settings

S = struct;


if convertCharsToStrings(getenv('COMPUTERNAME')) == convertCharsToStrings('MANHASSET')
    S.dir.imagingRootDir = 'G:\2PData';
elseif convertCharsToStrings(getenv('COMPUTERNAME')) == convertCharsToStrings('SUNHP')
%     imagingRootDir = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\suite2p';
    S.dir.imagingRootDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy';
    S.dir.codeRootDir = 'C:\Users\Xiaonan Richard Sun\Dropbox\Users\Richard\matlab';
    S.dir.bhvRootDir = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon';
elseif convertCharsToStrings(computer) == convertCharsToStrings('GLNXA64') || isunix == 1
    S.dir.imagingRootDir = '/grid/churchland/data/data/richard_s2p_npy';
    S.dir.codeRootDir = '/grid/churchland/home/xisun/matlab';
    S.dir.bhvRootDir = '/grid/churchland/data/data/Behavior_Simon';
end

S.dir.imagingSubDir = fullfile('suite2p', 'plane0');
S.dir.bhvSubDir = fullfile('SpatialDisc', 'Session Data');

S.segIdx = [1 0.75 1.25 0.5 1]; 