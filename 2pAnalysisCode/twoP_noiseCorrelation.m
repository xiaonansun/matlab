%% Noise correlation
% n1_idx=70; n2_idx=70; % example neuron indices
% rho=[0.25 0.75];
noise_int = [find(lickWinIdx==-5) find(lickWinIdx==-1)]; % time interval (in indices relative to lick, or zero) over which a mean of inferred spiking activity will be computed
noise = permute(mean(dataLick(:,noise_int(1):noise_int(end),:),2),[1 3 2])'; 
noise_red=[noise(:,data.idx_redcell) noise(:,data.idx_notredcell)];
C = cov(noise); C_red=cov(noise_red);
CU=tril(C,-1); CU=CU(:); CU(CU==0)=[];
C_redU=tril(C_red,-1); C_redU=C_redU(:); C_redU(C_redU==0)=[];

hF1= figure(1);
image(C);

%Sorted by tdTomato+ (upper left) versus tdTomato- neurons
image(C_red); hold on;
hLineHor=line([length(data.idx_redcell) length(data.idx_redcell)], [0 length(noise)]); % Draw a line between red and non-red cells 
hold on;
hLineVert=line([0 length(noise)], [length(data.idx_redcell) length(data.idx_redcell)]); % Draw a line between red and non-red cells
set(hLineHor,'LineWidth',1.5,'Color',[1 1 1],'LineStyle',':');
set(hLineVert,'LineWidth',1.5,'Color',[1 1 1],'LineStyle',':');

coord=zeros(length(npy.stat),2);
for i = 1:length(npy.stat)
    coord(i,:)=npy.stat{i}.med;
end
coord(npy.iscell(:,1)==0,:)=[];

dist=zeros(length(coord));
for i = 1:length(coord)
    for j = 1:length(coord)
    dist(i,j) = norm(coord(i,:)-coord(j,:));
    end
end
distU=tril(dist,-1); distU=distU(:); distU(distU==0)=[];
hF2=figure(2);
image(dist);

hF3=figure;
plot(distU,CU,'.k');
% plot(distU,C_redU,'.k');

% R=corrcov(C);
% image(R);
% imagesc(C,min(C(:))+rho*(max(C(:))-min(C(:))));

% noise2 = mean(dataLick(n2_idx,noise_int(1):noise_int(end),:),2);
% plot(noise1,noise2,'.k');

%% Cross-correlation
allCells=npy.spks(npy.iscell(:,1)==1,:);
size(allCells)

r = xcorr(allCells');
plot(r,'.-k');