function rateDisc_qrGlobalV(animal,ndims)
%code to compute QR decomposition on whole-frame V

%% some variables
if ~exist('ndims', 'var') || isempty(ndims)
    ndims = 500;
end
if ispc
    cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\';
%     cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\';
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
end
bPath = [cPath animal filesep 'blockData' filesep]; % path for global dimensions

cFile = matfile([bPath 'wV.mat']); %check wV file to load subset of components
wV = cFile.wV(1:ndims,:,:);

wV = reshape(wV, ndims, []); %make sure wV is in 2D
nanIdx = isnan(wV(1,:)); %keep index for NaN frames
wV(:, nanIdx) = []; %remove NaNs
[q, r] = qr(wV', 0); %get qr results

save([bPath 'wQR.mat'], 'q', 'r', 'nanIdx', '-v7.3');
