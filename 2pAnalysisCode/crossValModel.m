    function [Vm, cBeta, cR, subIdx, cRidge, cLabels, cMap, cMovie] =  crossValModel(cLabels)
        
        cIdx = ismember(regIdx(~rejIdx), find(ismember(regLabels,cLabels))); %get index for task regressors
        cLabels = regLabels(sort(find(ismember(regLabels,cLabels)))); %make sure motorLabels is in the right order
        
        %create new regressor index that matches motor labels
        subIdx = regIdx(~rejIdx);
        subIdx = subIdx(cIdx);
        temp = unique(subIdx);
        for x = 1 : length(temp)
            subIdx(subIdx == temp(x)) = x;
        end
        cR = fullR(:,cIdx);
        
        Vm = zeros(size(Vc),'single'); %pre-allocate motor-reconstructed V
        randIdx = randperm(size(Vc,2)); %generate randum number index
        foldCnt = floor(size(Vc,2) / ridgeFolds);
        cBeta = cell(1,ridgeFolds);
        
        for iFolds = 1:ridgeFolds
            dataIdx = true(1,size(Vc,2));
            
            if ridgeFolds > 1
                dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
                if iFolds == 1
                    [cRidge, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', cR(dataIdx,:), true); %get beta weights and ridge penalty for task only model
                else
                    [~, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', cR(dataIdx,:), true, cRidge); %get beta weights for task only model. ridge value should be the same as in the first run.
                end
                Vm(:,~dataIdx) = (cR(~dataIdx,:) * cBeta{iFolds})'; %predict remaining data: 
                % cR(~dataIdx,:) is the held-out design matrix 
                %  cBeta{iFolds} are the betas of the training data
                
                if rem(iFolds,ridgeFolds/5) == 0
                    fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
                end
            else
                [cRidge, cBeta{iFolds}] = ridgeMML(Vc', cR, true); %get beta weights for task-only model.
                Vm = (cR * cBeta{iFolds})'; %predict remaining data
                disp('Ridgefold is <= 1, fit to complete dataset instead');
            end
        end
        
        % computed all predicted variance
%         Vc = reshape(Vc,size(Vc,1),[]);
%         Vm = reshape(Vm,size(Vm,1),[]);
%         if strcmpi(dType, 'Widefield')
%             if length(size(U)) == 3
%                 U = arrayShrink(U, squeeze(isnan(U(:,:,1))));
%             end
%             covVc = cov(Vc');  % S x S
%             covVm = cov(Vm');  % S x S
%             cCovV = bsxfun(@minus, Vm, mean(Vm,2)) * Vc' / (size(Vc, 2) - 1);  % S x S
%             covP = sum((U * cCovV) .* U, 2)';  % 1 x P
%             varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
%             varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
%             stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
%             cMap = gather((covP ./ stdPxPy)');
%         elseif strcmpi(dType, 'twoP')
%             for x = 1 : size(Vc,1)
%                 cMap(x) = corr2(Vc(x,:),Vm(x,:)); % get correlation between data and model
%             end
%         end
%             
%         % movie for predicted variance
%         cMovie = NaN;
%         if strcmpi(dType, 'Widefield')
%             cMovie = zeros(size(U,1),frames, 'single');
%             for iFrames = 1:frames
%                 
%                 frameIdx = iFrames:frames:size(Vc,2); %index for the same frame in each trial
%                 cData = bsxfun(@minus, Vc(:,frameIdx), mean(Vc(:,frameIdx),2));
%                 cModel = bsxfun(@minus, Vm(:,frameIdx), mean(Vm(:,frameIdx),2));
%                 covVc = cov(cData');  % S x S
%                 covVm = cov(cModel');  % S x S
%                 cCovV = cModel * cData' / (length(frameIdx) - 1);  % S x S
%                 covP = sum((U * cCovV) .* U, 2)';  % 1 x P
%                 varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
%                 varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
%                 stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
%                 cMovie(:,iFrames) = gather(covP ./ stdPxPy)';
%                 clear cData cModel
%                 
%             end
%         end
%         fprintf('Run finished. RMSE: %f\n', median(cMovie(:).^2));
        
        
    end
