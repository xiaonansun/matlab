function [normingImg, softConst] = brightenFilter2DGauss(inputImg,params)
% Takes inputImg and computes a low-pass filtered version of it plus a
% constant value (softConst). The output (normingImg) can be used to
% normalize inputImg (or a movie) to make its pixel intensities more
% unifrom. 
% Run the following to assess pixel intensities after
% normalization: figure; imagesc(inputImg./(normingImg))
% 
% Example inputs:
% inputImg = medImage{2};
% params.brightConstant = 'Quantile'; % 3000; % whether to use a constant value for softConst or to use inputImg's quantile to set softConst.
% params.brightConstantQuantile = .5; % if using quantile, what quantile?
% params.brightGaussSDs = 30 * [512 512] ./ [585 585]; % 30 microns; scan voltage: 4 (FN); % Guassian width in pixels. 50 * [512 402] ./ [710 690]; % 50 microns; scan voltage: 4.9 (MTK)


%%
gaussDev = params.brightGaussSDs;

% Filter image using gauss2D
normingImg = imfilter(inputImg, gauss2D(gaussDev .* 4, gaussDev), 'symmetric');

if ischar(params.brightConstant) && strcmp(params.brightConstant, 'Quantile')
    softConst = quantile(inputImg(:), params.brightConstantQuantile);
else
    softConst = params.brightConstant;
end
normingImg = normingImg + softConst;


