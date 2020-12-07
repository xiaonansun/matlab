function hf = laserTrigPETHM(M)
% hf = laserTrigPETHM(M)
%
% Make a PETH with 1-ms bins, with the locking event being the laser
% pulses.

%% Error checking

% Make a dummy figure if there's no PETH
if isempty(M.laserPETH)
  hf = figure;
  textTitle('No PETH');
  return;
end

%% Make the histogram

hf = figure;
% histc is annoying; the last value contains matches to the right edge of
% the highest bin. We'll exclude those spikes.

hollowHist(M.laserPETH, M.times.laserPETHEdges, 'k', 1,1)
%bar(M.times.laserPETHEdges, M.laserPETH, 'histc');
%plot(M.times.laserPETHEdges, M.laserPETH, 'c-');
set(gca, 'Box', 'off');
set(gca, 'TickDir', 'out');
xlim(M.times.laserPETHEdges([1 end]));
xlabel('Time from laser pulse (ms)');

titl = sprintf('Tetrode %d, cluster %d, rating %0.1f, SNR %1.1f', M.tetrode, M.cluster, M.rating, M.SNR);
textTitle(titl);
