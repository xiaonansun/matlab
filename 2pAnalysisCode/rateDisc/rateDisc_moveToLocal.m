function rateDisc_moveToLocal(animal)
% code to move data from nlsas back to grid. This is to perform blockwise
% SVD over all recordings or run modeling / cross-validation code on HPC.

%% some variables
% cPath = ['\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\' animal filesep]; %path to nlsas on windows
% tPath = ['\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\' animal filesep]; %path to grid on windows
% cPath = ['\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\' animal filesep]; %path to grid on windows
cPath = ['\\churchlandNAS\homes\DOMAIN=CSHL\smusall\BpodImager\Animals\' animal filesep ]; %path to churchlandNAS
tPath = ['Q:\BpodImager\Animals\' animal filesep ]; %local path
trainingRange = 'allAudio'; %use this to run analysis only in a certain range of training

%% get list of recordings and reference image

recs = rateDisc_getRecsForAnimal(animal, trainingRange);
% for iRecs = 1 : length(recs)
for iRecs = 1 : length(recs)
    try
        % get data and cut to size
        fPath = [cPath 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
        ftPath = [tPath 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
        if ~exist(ftPath, 'dir')
            mkdir(ftPath);
        end
%         if ~exist([ftPath 'BehaviorVideo'], 'dir')
%             mkdir([ftPath 'BehaviorVideo']);
%         end        

        delete([ftPath 'AC.mat']);
        delete([ftPath 'largeAC.mat']);
        delete([ftPath 'largeAC_20_50.mat']);
        copyfile([fPath 'newAC_20_50.mat'], [ftPath 'newAC_20_50.mat']);
        
%         copyfile([fPath 'largeAC_20_50.mat'], [ftPath 'largeAC_20_50.mat']);
%         copyfile([fPath 'largeAC.mat'], [ftPath 'largeAC.mat']);
        copyfile([fPath 'alignU.mat'], [ftPath 'alignU.mat']);
        copyfile([fPath 'QR.mat'], [ftPath 'QR.mat']);
%         copyfile([fPath 'Vc.mat'], [ftPath 'Vc.mat']);
        copyfile([fPath 'mask.mat'], [ftPath 'mask.mat']);
%         copyfile([fPath 'opts2.mat'], [ftPath 'opts2.mat']);
%         copyfile([fPath 'orgdimBeta.mat'], [ftPath 'orgdimBeta.mat']);
%         copyfile([fPath 'vidRegData.mat'], [ftPath 'vidRegData.mat']);
%         copyfile([fPath 'orgVidBeta.mat'], [ftPath 'orgVidBeta.mat']);
%         copyfile([fPath 'interpVc.mat'], [ftPath 'interpVc.mat']);
%         copyfile([fPath 'motorBeta.mat'], [ftPath 'motorBeta.mat']);
%         copyfile([fPath 'spontMotorBeta.mat'], [ftPath 'spontMotorBeta.mat']);
%         copyfile([fPath 'motorregData.mat'], [ftPath 'motorregData.mat']);
%         copyfile([fPath 'orgregData.mat'], [ftPath 'orgregData.mat']);
%         copyfile([fPath 'spontMotorregData.mat'], [ftPath 'spontMotorregData.mat']);


        %         % move behavior data for model
%         analogFiles = dir([fPath 'Analog_*']);
%         for x = 1 : length(analogFiles)
%             copyfile([fPath analogFiles(x).name], [ftPath analogFiles(x).name]);
%         end
%         bhvFile = dir([fPath filesep animal '_SpatialDisc*.mat']);
%         copyfile([fPath filesep bhvFile(1).name], [ftPath filesep bhvFile(1).name]);
%         copyfile([fPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'], [ftPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat']);
%         copyfile([fPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat'], [ftPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat']);
%         copyfile([fPath 'BehaviorVideo' filesep 'bhvOpts.mat'], [ftPath 'BehaviorVideo' filesep 'bhvOpts.mat']);
%         copyfile([fPath 'BehaviorVideo' filesep 'FilteredPupil.mat'], [ftPath 'BehaviorVideo' filesep 'FilteredPupil.mat']);
%         copyfile([fPath 'BehaviorVideo' filesep 'segInd1.mat'], [ftPath 'BehaviorVideo' filesep 'segInd1.mat']);
%         copyfile([fPath 'BehaviorVideo' filesep 'segInd2.mat'], [ftPath 'BehaviorVideo' filesep 'segInd2.mat']);
%         copyfile([fPath 'BehaviorVideo' filesep 'SVD_Cam1-Seg1.mat'], [ftPath 'BehaviorVideo' filesep 'SVD_Cam1-Seg1.mat']);
%         copyfile([fPath 'BehaviorVideo' filesep 'SVD_Cam1-Seg2.mat'], [ftPath 'BehaviorVideo' filesep 'SVD_Cam1-Seg2.mat']);
                
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
    catch WE
        disp(WE.message);
        disp(['Error in rec: ' animal '-' num2str(iRecs) '-' recs(iRecs).name]);
    end
end

