for iRecs = 1 : length(recs)
fPath = [cPath cAnimal filesep 'SpatialDisc' filesep recs(iRecs).name filesep]; %Widefield data path
load([fPath 'blueAvg.mat'],'blueAvg')
load([fPath 'Vc.mat'],'Vc','U')

subplot(1,2,1); imagesc(blueAvg); colormap(cMap); axis image;
title(fPath);

Vc = reshape(Vc,200,[]);
Vc = Vc(:,~isnan(Vc(1,:)));

[~, varP1] = rateDisc_modelCorr(Vc, Vc, U);
subplot(1,2,2); imageScale(arrayShrink(varP1, squeeze(isnan(U(:,:,1))), 'split'));
axis image;colormap(gca,colormap_blueblackred(256)); title(iRecs);

pause
end
