function cOpts = rateDisc_alignSingleSession(blueAvg, alignPic, modelPath, cOpts, cropRange)
%% align image using LEAP estimates
% reduce and reshape image
blueAvg = arrayResize(blueAvg, 2);
blueAvg(272, :, :) = 0;
blueAvg = reshape(blueAvg, size(blueAvg,1), size(blueAvg,2), 1, []);

%make prediction
preds = predict_box(uint8(blueAvg ./ 65535 .* 255), modelPath);
lamb2breg = squeeze(preds.positions_pred(4, :, :) - preds.positions_pred(5, :, :));
lamb2front = squeeze(preds.positions_pred(2, :, :) - preds.positions_pred(5, :, :));
bregfront = squeeze(preds.positions_pred(2, :, :) - preds.positions_pred(4, :, :));

cRotate(1) = atan2(lamb2breg(2), lamb2breg(1)) .* 180 ./ pi .* -1;
cRotate(2) = atan2(lamb2front(2), lamb2front(1)) .* 180 ./ pi .* -1;
cRotate(3) = atan2(bregfront(2), bregfront(1)) .* 180 ./ pi .* -1;
cRotate = 90 - mean(cRotate); %mean angle of center line

bregmaPos = preds.positions_pred(4, :); %estimated position of bregma
lambdaPos = preds.positions_pred(5, :); %estimated position of lambda
frontPos = preds.positions_pred(2, :); %estimated position of front center

%% use predictions to identify anatomical markers
temp = blueAvg; %keep this to add index before rotation
dSize = size(temp);
[xx,yy] = meshgrid(1:dSize(2),1:dSize(1)); %isolate index for selected area
mask = false(dSize(1:2));
mask = mask | hypot(xx - bregmaPos(1), yy - bregmaPos(2)) <= 2;
temp(mask) = Inf;

mask = false(dSize(1:2));
mask = mask | hypot(xx - lambdaPos(1), yy - lambdaPos(2)) <= 2;
mask = mask | hypot(xx - frontPos(1), yy - frontPos(2)) <= 2;
temp(mask) = NaN;

temp = imrotate(temp, cRotate, 'crop');
xCenter = round(mean(ceil(find(isnan(temp) | isinf(temp))/size(temp,1) ))); %find center line
yCenter = rem(find(isinf(temp),1), size(temp,1)); %find bregma
cTranslate = [round(size(temp,2)/2) - xCenter, round(size(temp,1)/2) - yCenter];

%% apply to image
blueAvg = imrotate(blueAvg, cRotate, 'crop');
blueAvg = imtranslate(blueAvg, cTranslate);
blueAvg = imresize(blueAvg, 2);
blueAvg = blueAvg(1:size(alignPic,1), 1:size(alignPic,2));

%% align to alignPic to reduce variability
options.TolX = 0.0001;
nRotate = fminsearch(@(u) Widefield_alignImage(alignPic, blueAvg, u, cropRange), 1, options);
[~, nTranslate] = Widefield_alignImage(alignPic, blueAvg, nRotate, 0);

% add alignment values to opts file for output
cOpts.transParams.angleD = cRotate + nRotate;
cOpts.transParams.tC = round(cTranslate + flipud(nTranslate));
