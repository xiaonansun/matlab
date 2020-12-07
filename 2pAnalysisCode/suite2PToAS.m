function [A, S, neuropilCoeffs, kernels, recon] = suite2PToAS(dat, showPlots)
% [A, S, neuropilCoeffs, kernels, reconstruction] = suite2PToAS(dat [, showPlots])
% 
% Convert Suite2P's dat data structure to something easier to use. 
% A is the spatial footprints of the neurons, pixels x neurons. 
% S is the inferred "spikes", neurons x time. 
% neuropilCoeffs is the scaling coefficient for the neuropil.  (i.e. 40x1
% double)
% kernels is the calcium convolution kernel. 
% recon is the reconstructed traces for the neurons, obtained by re-convolving the
% inferred spiking with the kernel and regressing to get the baseline and
% scaling. I don't understand the absolute scaling of these values, so it
% should be converted to dF/F later.
% 
% Optional showPlots shows some sanity-check plots. Default false.


%% Optional arguments

if ~exist('showPlots', 'var')
  showPlots = 0;
end


%% Retrieve neuropil coefficients

roiIs = find(dat.cl.iscell(dat.cl.isroi));

neuropilCoeffs = cellfun(@(d) d.B(3), dat.cl.dcell);
neuropilCoeffs = neuropilCoeffs(roiIs);

if showPlots
  figure;
  histogram(neuropilCoeffs);
  title('Neuropil coefficients -- should be <1!');
end

cellIs = find(dat.cl.iscell & dat.cl.isroi);


%% Display background image

CaIm = dat.mimg(:, :, dat.maxmap);

if showPlots
  figure; imagesc(CaIm);
  axis off
  axis equal
  colormap('gray');
end


%% Produce A matrix

A = zeros(numel(CaIm), length(cellIs));
for i = 1:length(cellIs)
  these = (dat.res.iclust(:) == cellIs(i));
  A(these, i) = dat.res.lambda(these);
end
A = sparse(A);

if showPlots
  COMs = fastCOMsA(A, size(CaIm));
  plotCaMDS([], CaIm, COMs, A, randperm(size(A, 2)));
end


%% Produce S matrix

S = zeros(length(cellIs), size(dat.Fcell{1}, 2));
for i = 1:length(roiIs)
  r = roiIs(i);
  S(i, dat.cl.dcell{r}.st) = dat.cl.dcell{r}.c;
end


%% Pull out kernels

if nargout > 3
  kernelLen = length(dat.cl.dcell{roiIs(1)}.kernel);
  kernels = NaN(kernelLen, length(roiIs));
  for r = 1:length(roiIs)
    kernels(:, r) = dat.cl.dcell{r}.kernel;
  end
end


%% Reconstruct fluorescence traces

if nargout > 4
  T = size(dat.Fcell{1}, 2);

  % Convolve spikes with kernels
  recon = NaN(length(cellIs), T);
  for c = 1:length(cellIs)
    thisConv = conv(S(c, :)', kernels(:, c));
    recon(c, :) = thisConv(1:T);
  end
  
  % Find scalings
  baselines = NaN(length(roiIs), 1);
  scalings = NaN(length(roiIs), 1);
  for c = 1:length(cellIs)
    b = regress(dat.Fcell{1}(cellIs(c), :)', [ones(T, 1), recon(c, :)']);
    baselines(c) = b(1);
    scalings(c) = b(2);
  end
  
  recon = baselines + recon .* scalings;
end


