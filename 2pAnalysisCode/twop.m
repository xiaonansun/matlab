classdef twop
    properties
        dir_path_imaging;
        dir;
        animal;
        session;
    end
    methods
        function obj = s2p_output_dir(obj)
%             class constructor
            %%
            obj.dir.imagingSubDir1 = 'imaging';
            obj.dir.imagingSubDir2 = fullfile('suite2p', 'plane0');
            obj.dir.bhvSubDir = fullfile('SpatialDisc', 'Session Data');
            if convertCharsToStrings(getenv('COMPUTERNAME')) == convertCharsToStrings('MANHASSET')
                obj.dir.imagingRootDir = 'G:\2PData';
            elseif convertCharsToStrings(getenv('COMPUTERNAME')) == convertCharsToStrings('SUNHP')
                %     imagingRootDir = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\suite2p';
                obj.dir.imagingRootDir = '\\grid-hs\churchland_nlsas_data\data\richard_s2p_npy';
                obj.dir.codeRootDir = 'C:\Users\Xiaonan Richard Sun\Dropbox\Users\Richard\matlab';
                obj.dir.bhvRootDir = '\\grid-hs\churchland_nlsas_data\data\Behavior_Simon';
                obj.dir.bhvVidRootDir = '\\grid-hs\churchland_nlsas_data\BehaviorVideo';
                obj.dir.bhvVidSubDir = 'SpatialDisc\Session Data';
                obj.dir.isUnix = false;
            elseif convertCharsToStrings(computer) == convertCharsToStrings('GLNXA64') || isunix == 1
                obj.dir.imagingRootDir = '/grid/churchland/data/data/richard_s2p_npy';
                obj.dir.codeRootDir = '/grid/churchland/home/xisun/matlab';
                obj.dir.bhvRootDir = '/grid/churchland/data/data/Behavior_Simon';
                obj.dir.bhvVidRootDir = '/grid/churchland/data/BehaviorVideo';
                obj.dir.bhvVidSubDir = 'SpatialDisc/Session Data';
                obj.dir.isUnix = true;
            end
            obj.dir_path_imaging = fullfile(obj.dir.imagingRootDir,obj.animal,obj.dir.imagingSubDir1,obj.session,obj.dir.imagingSubDir2);
%             obj.dir_analysis = fullfile(obj.dir.imagingRootDir,'analysis');
            obj.dir = obj.dir;
        end
    end
end