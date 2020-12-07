function meanVoltByFrame = binMOMAnalog(volt)
% meanVoltByFrame = binMOMAnalog(volt)
% 
% Input is the voltages from the Sutter MOM analog channels (returned by
% readMOMAnalog() ). Parses the slow galvo trace (second analog input,
% recorded at 1 kHz) to figure out how many analog samples there are per
% frame, then averages each voltage trace over the frame. If the input is
% nChannels x nSamples, the result is nChannels x nFrames.

% Find out how many analog samples there were per frame by parsing the slow
% galvo trace
slowGalvo = volt(2, :);
[samplesPerFrame, offsetSamples] = samplesPerFrameFromSlowGalvo(slowGalvo);

% Find the start and end of each frame, in samples
sampEnds = cumsum(samplesPerFrame);
sampStarts = [1 sampEnds(1:end-1)+1];
sampEnds = sampEnds + offsetSamples;
sampStarts = sampStarts + offsetSamples;

% The last frame is occasionally a sample too long, because of how the
% system stops. If so, fix it.
if sampEnds(end) > size(volt, 2)
  % Make sure the error isn't too bad. Value "5" chosen arbitrarily.
  if sampEnds(end) - size(volt, 2) > 5
    warning('binMOMAnalog: some issue overrunning the end of the trace with the last frame, %d samples', sampEnds(end) - size(volt, 2));
  end
  sampEnds(end) = size(volt, 2);
end

nFrames = length(sampStarts);

% Bin by taking the mean
meanVoltByFrame = NaN(size(volt, 1), nFrames);
for fr = 1:nFrames
  meanVoltByFrame(:, fr) = mean(volt(:, sampStarts(fr):sampEnds(fr)), 2);
end
