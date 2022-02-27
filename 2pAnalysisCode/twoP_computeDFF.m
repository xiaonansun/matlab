function dF_traces = twoP_computeDFF(animal,session,interval)

%%
% animal = ; session = ;

S = twoP_settings;
npy = twoP_importSuite2pNPY(animal,session);
fs = 1000/S.msPerFrame;
% npy.ops.sig_baseline;
%%
idx = [40000 80000];
idxCell = logical(npy.iscell(:,1));
Fneu = npy.ops.neucoeff*smoothdata(npy.Fneu(idxCell,:),2,'gaussian',npy.ops.win_baseline*fs);

dF = bsxfun(@minus,npy.F(idxCell,:),Fneu);
    dF_traces = bsxfun(@rdivide,dF,Fneu);

if ~exist('interval','var') || isempty(interval)
    save_path = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,'dF_traces.mat');
else
    dF_traces = dF_traces(:,interval(1):interval(2));
    save_path = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir,...
        ['dF_traces_' num2str(interval(1)) '_' num2str(interval(2)) '.mat']);
end

save(save_path,'dF_traces');

% figure
% plot(npy.F(:,1:idxEnd)');
% figure
% plot(dFF(:,1:idxEnd)');

imagesc(dF_traces);