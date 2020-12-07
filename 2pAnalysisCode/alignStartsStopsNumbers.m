function [starts, stops, numbers] = alignStartsStopsNumbers(stimOnSamples, codeStarts, trialNumbers, trialStopSamples)

starts = NaN(1, length(stimOnSamples));
numbers = NaN(1, length(stimOnSamples));
stops = NaN(1, length(stimOnSamples));

startI = 1;
stopI = 1;

for tr = 1:length(stimOnSamples)
  if startI == length(codeStarts)
    error('Ran off end of codeStarts too early');
  end
  
  while startI < length(codeStarts) && codeStarts(startI + 1) < stimOnSamples(tr)
    startI = startI + 1;
  end
  
  while stopI < length(trialStopSamples) && trialStopSamples(stopI) < stimOnSamples(tr)
    stopI = stopI + 1;
  end
  
  if trialStopSamples(stopI) > stimOnSamples(tr)
      starts(tr) = codeStarts(startI);
      numbers(tr) = trialNumbers(startI);
      stops(tr) = trialStopSamples(stopI);
  end
end

if any(numbers < 0)
  error('%d bad trial codes! Likely means the corresponding trial start times are wrong', sum(numbers < 0));
end

if any(diff(numbers) < 0)
    warning('Non-monotonic trial codes! All trials before last non-monotic rise will be rejected!');
    idx = find(diff(numbers) < 0);
    starts = starts(idx(end)+1 : end);
    stops = stops(idx(end)+1 : end);
    numbers = numbers(idx(end)+1 : end);
end