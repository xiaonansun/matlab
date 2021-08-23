function rateDisc_moveToGrid(animal)
% code to move data from nlsas back to grid. This is to perform blockwise
% SVD over all recordings or run modeling / cross-validation code on HPC.

%% some variables
cPath = ['\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\' animal filesep]; %path to nlsas on windows
tPath = ['\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\' animal filesep]; %path to grid on windows
trainingRange = 'allAudio'; %use this to run analysis only in a certain range of training

%% get list of recordings and reference image

recs = rateDisc_getRecsForAnimal(animal, trainingRange);
for iRecs = 1 : length(recs)
    try
        % get data and cut to size
        fPath = [cPath 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
        ftPath = [tPath 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
        if ~exist(ftPath, 'dir')
            mkdir(ftPath);
        end
        if ~exist([ftPath 'BehaviorVideo'], 'dir')
            mkdir([ftPath 'BehaviorVideo']);
        end        
%         copyfile([fPath 'Vc.mat'], [ftPath 'Vc.mat']);
%         copyfile([fPath 'hemoV.mat'], [ftPath 'hemoV.mat']);

%         copyfile([fPath 'Vc.mat'], [ftPath 'Vc.mat']);
%         copyfile([fPath 'opts2.mat'], [ftPath 'opts2.mat']);
%         copyfile([fPath 'blueAvg.mat'], [ftPath 'blueAvg.mat']);
% %         copyfile([fPath 'rsVc.mat'], [ftPath 'rsVc.mat']);
%         copyfile([fPath 'mask.mat'], [ftPath 'mask.mat']);
%         
%         % move behavior data for model
        analogFiles = dir([fPath 'Analog_*']);
        for x = 1 : length(analogFiles)
            copyfile([fPath analogFiles(x).name], [ftPath analogFiles(x).name]);
        end
        bhvFile = dir([fPath filesep animal '_SpatialDisc*.mat']);
        copyfile([fPath filesep bhvFile(1).name], [ftPath filesep bhvFile(1).name]);
%         copyfile([fPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'], [ftPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat']);
%         copyfile([fPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat'], [ftPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat']);
        copyfile([fPath 'BehaviorVideo' filesep 'bhvOpts.mat'], [ftPath 'BehaviorVideo' filesep 'bhvOpts.mat']);
        copyfile([fPath 'BehaviorVideo' filesep 'FilteredPupil.mat'], [ftPath 'BehaviorVideo' filesep 'FilteredPupil.mat']);
        copyfile([fPath 'BehaviorVideo' filesep 'segInd1.mat'], [ftPath 'BehaviorVideo' filesep 'segInd1.mat']);
        copyfile([fPath 'BehaviorVideo' filesep 'segInd2.mat'], [ftPath 'BehaviorVideo' filesep 'segInd2.mat']);
        copyfile([fPath 'BehaviorVideo' filesep 'SVD_Cam1-Seg1.mat'], [ftPath 'BehaviorVideo' filesep 'SVD_Cam1-Seg1.mat']);
        copyfile([fPath 'BehaviorVideo' filesep 'SVD_Cam1-Seg2.mat'], [ftPath 'BehaviorVideo' filesep 'SVD_Cam1-Seg2.mat']);
                
        %         try
        %             copyfile([fPath 'interpVc.mat'], [ftPath 'interpVc.mat']);
        %             copyfile([fPath 'dimBeta.mat'], [ftPath 'dimBeta.mat']);
        %             copyfile([fPath 'regData.mat'], [ftPath 'regData.mat']);
        %             copyfile([fPath 'orgdimBeta.mat'], [ftPath 'orgdimBeta.mat']);
        %             copyfile([fPath 'orgregData.mat'], [ftPath 'orgregData.mat']);
        %             copyfile([fPath 'opts.mat'], [ftPath 'opts.mat']);
        %             copyfile([fPath 'opts2.mat'], [ftPath 'opt2s.mat']);
        %         catch
        %             disp('Failed to copy model results to grid.')
        %  be       end
        disp([animal '-' num2str(iRecs) '-' recs(iRecs).name]);
    catch
        disp(['Error in rec: ' animal '-' num2str(iRecs) '-' recs(iRecs).name]);
    end
end

