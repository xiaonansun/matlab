function [corrMat, varP1, varP2] = modelCorr(Vc,Vm,U)
% 2022-12-19: copy-pasted from Simon's repository.
% - Obviously for 2P data U does not exist
%short code to compute the correlation between lowD data Vc and modeled
%lowD data Vm. Vc and Vm are temporal components, U is the spatial
%components. corrMat is a the correlation between Vc and Vm in each pixel.

%%
Vc = reshape(Vc,size(Vc,1),[]);
Vm = reshape(Vm,size(Vm,1),[]);

if ~exist('U','var') || isempty(U)
    U = [];
else
    if length(size(U)) == 3
        U = arrayShrink(U, squeeze(isnan(U(:,:,1))));
    end
end

if ~exist('U','var') || isempty(U)
    covVc = cov(Vc');  % S x S
    covVm = cov(Vm');  % S x S
    cCovV = bsxfun(@minus, Vm, mean(Vm,2)) * Vc' / (size(Vc, 2) - 1);  % S x S
    covP = diag(cCovV)';  % 1 x P
    varP1 = diag(covVc)';  % 1 x P
    varP2 = diag(covVm)';  % 1 x P
    stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
    corrMat = gather((covP ./ stdPxPy)');
else
    covP = sum((U * cCovV) .* U, 2)';  % 1 x P
    varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
    varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
    stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
    corrMat = gather((covP ./ stdPxPy)');
end