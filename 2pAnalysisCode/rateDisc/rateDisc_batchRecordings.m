function rateDisc_batchRecordings
% code to run over all recordings and do something.


% cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %path to grid on windows
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server


% get animal info
[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings;
trainingRange = 'audioDisc'; %use this to run analysis only in a certain range of training
animals = dataOverview(:,1);
animals = animals(1:4);

disp(cPath);
for iAnimals = 1 : length(animals)
    
    cAnimal = animals{iAnimals};

    % get recs
    recs = dir([cPath cAnimal filesep 'SpatialDisc' filesep]);
    recs = rateDisc_filterRecs(trainDates{ismember(dataOverview(:,1), cAnimal)}, recs, trainingRange); %this sorts recordings by date

    for iRecs = 1 : length(recs)
        disp([cAnimal ' - ' recs(iRecs).name]);
        
        try
            
            rateDisc_choiceModel(cPath, cAnimal, recs(iRecs).name, 'Widefield');
            
%         fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs(iRecs).name filesep]; %session data path
%         
%         if ~exist([fPath 'opts.mat'], 'file')
%             copyfile([fPath 'opts2.mat'], [fPath 'opts.mat']);
%         end
%         
%         load([fPath 'opts.mat'], 'opts');
%         
%         if opts.frameRate == 30 && ~exist([fPath 'rsVc.mat'],'file')
%             
%             load([fPath 'hemoV.mat'], 'hemoV');
%             load([fPath 'blueV.mat'], 'blueV','U','blueFrametimes','Sv','totalVar','trials','bTrials');
%             
%             U = U(:, :, 1:200);
%             blueV = blueV(1:200, :, :);
%             hemoV = hemoV(1:200, :, :);
%             
%             % downsample to 15Hz
%             blueV = rateDisc_downsampV(blueV,2);
%             hemoV = rateDisc_downsampV(hemoV,2);
%             opts.frameRate = 15;
%             
%             [Vc, regC, T, hemoVar] = Widefield_SvdHemoCorrect(U, blueV, hemoV, opts.frameRate);
%             
%             save([fPath 'rsVc.mat'],'Vc','U','blueFrametimes','Sv','totalVar','trials','bTrials','-v7.3');
%             save([fPath 'rsHemoCorrection.mat'],'regC','T', 'hemoVar')
%             
%         end

        catch
            disp('Resampling failed');
        end
    end
end