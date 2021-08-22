function behavior_PlotPerformanceMatrix(animal)
% find(~cellfun(@isempty,strfind({bhv.SessionName},'Jan'))) % To find all sessions for a given month

bhv = behavior_LoadAllSessions(animal);

%%

trialCounts = zeros(length(bhv),5);
tic
for i = 1:length(bhv)
    trialCounts(i,1)=bhv(i).SessionData.cTrial;
    trialCounts(i,2)=length(bhv(i).SessionData.Rewarded);
    trialCounts(i,3)=length(bhv(i).SessionData.AutoReward);
    trialCounts(i,4)=length(bhv(i).SessionData.SingleSpout);
    trialCounts(i,5)=length(bhv(i).SessionData.Assisted);
end
toc

%%

eventRange = [10 30];
vEvents = eventRange(1):eventRange(end);
matEvents = zeros(length(vEvents),length(vEvents));
events = zeros(sum(trialCounts(:,2)),2 );

tic
j=0;
for i = 1:length(bhv)
    nTrials = length(bhv(i).SessionData.Rewarded);
    E = bhv(i).SessionData.stimEvents;
    E = cellfun(@(x) x{:},[cellfun(@(x) x(1),E,'UniformOutput',false)' cellfun(@(x) x(2),E,'UniformOutput',false)'],'UniformOutput',false);
    events(j+1:sum(trialCounts(1:i,2)),1:2)= cellfun(@numel,E);
    events(j+1:sum(trialCounts(1:i,2)),3)= bhv(i).SessionData.ResponseSide';
    events(j+1:sum(trialCounts(1:i,2)),4)= bhv(i).SessionData.Rewarded';
    events(j+1:sum(trialCounts(1:i,2)),5)= bhv(i).SessionData.AutoReward';
    events(j+1:sum(trialCounts(1:i,2)),6)= bhv(i).SessionData.SingleSpout';
    events(j+1:sum(trialCounts(1:i,2)),7)= round(cellfun(@max,{bhv(i).SessionData.TrialSettings.DistFractions}))';
%     events(j+1:sum(trialCounts(1:i,2)),6)= bhv(i).SessionData.Assisted';
    j=j+nTrials;
end
toc

mEvents = events(~events(:,5),:);
mEvents = mEvents(~mEvents(:,6),:);
mEvents = mEvents(mEvents(:,7)==1,:);
% mEvents = mEvents(~mEvents(:,6),:);

for i = 1:length(vEvents)
    iEvents = mEvents(mEvents(:,1)==vEvents(i),:);
    for j = 1:length(vEvents)
        jEvents = iEvents(iEvents(:,2)==vEvents(j),:);
        matEvents(i,j) = sum(jEvents(:,4))./numel(jEvents(:,4));
        clear jEvents;
    end
    clear iEvents;
end

figure('Name',[animal ': all sessions']);
imagesc(matEvents); colormap jet;
h=gca; h.XAxis.TickLength = [0 0]; h.YAxis.TickLength = [0 0];
xticks([1:length(vEvents)]);
xticklabels(cellstr(num2str(vEvents'))');
yticks([1:length(vEvents)]);
yticklabels(cellstr(num2str(vEvents'))');
xlabel('Number of right-sided events');
ylabel('Number of left-sided events');

