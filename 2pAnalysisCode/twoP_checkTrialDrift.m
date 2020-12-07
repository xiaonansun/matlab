function [rejIdx, rampVal] = twoP_checkTrialDrift(data, lowLim, highLim, trialIdx)

if ~exist('trialIdx','var') || isinf(sum(trialIdx)) || isempty(trialIdx)
    trialIdx = 1:size(data.neural,3);
end
    
% Firing rates over trials
traces = squeeze(mean(data.neural(:,:,trialIdx), 2))';
temp1 = mean(traces(1:round(size(traces,1)/2), :)); %mean amp in the first half of the session
temp2 = mean(traces(round(size(traces,1)/2):end,:)); %mean amp in the second half of the session
rampVal = temp1 ./ (temp1 + temp2);
lowIdx = rampVal < lowLim;
highIdx = rampVal > highLim;
rejIdx = highIdx | lowIdx;

figure; 
subplot(2,3,1); plot(traces ./ max(traces)); title('Norm. FR per trial');
xlim([0 size(traces,1)]); axis square

subplot(2,3,2); hold on;
plot(smooth(mean(traces, 2)), 'linewidth', 2); 
plot(smooth(mean(traces(:, ~highIdx & ~lowIdx), 2)), 'linewidth', 2); 
title('Mean FR over neurons');
xlim([0 size(traces,1)]); axis square

subplot(2,3,3); histogram(temp1 ./ (temp1 + temp2)); title('Activity ratio of 1st v. 2nd half of session')
xlim([0 1]); vline(lowLim); vline(highLim); axis square

subplot(2,3,4); plot((traces(:, (lowIdx)))); title('Mean FR below lower limit');
xlim([0 size(traces,1)]); axis square

subplot(2,3,5); plot((traces(:, (highIdx)))); title('Mean FR above higher limit');
xlim([0 size(traces,1)]); axis square;

subplot(2,3,6); 
plot(smooth(mean(traces(:, (lowIdx)), 2)), 'linewidth', 2); hold on;
plot(smooth(mean(traces(:, (highIdx)), 2)), 'linewidth', 2); 
plot(smooth(mean(traces(:, ~highIdx & ~lowIdx), 2)), 'linewidth', 2); 
title('Mean FR over trials'); ylim([0 1]);
xlim([0 size(traces,1)]); axis square
legend({'Late high' 'Early high' 'Remaining cells'})
fprintf('Rejected %d/%d neurons\n', sum(rejIdx), size(traces,2))

figure
plot((traces(:, (~rejIdx)))); title('Mean FR , non-rejected neurons');
xlim([0 size(traces,1)]); axis square