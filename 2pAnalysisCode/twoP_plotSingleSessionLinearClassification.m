function twoP_plotSingleSessionLinearClassification(meta,lr)

%%
S = twoP_settings;

time_vec = 1:length(lr.cvAcc);

hFig = figure('Position',[500 500 900 500],...
    'visible','off');
hold on;

zero_vec = zeros(1,size(lr.cvAcc,2));
zero_vec(isnan(lr.cvAcc))=nan;

% bl = boundedline(time_vec, lr.cvAcc, ones(1,size(lr.cvAcc,2))*0,'nan','remove')

% boundedline(time_vec,mean(lr.cvAccNR),zeros(1,size(lr.cvAccNR,2)), ...
%     'g','nan','gap',...
%     'transparency',0);
line_label = {'tdT- (subsamp)','Mixed','tdT- (2X subsamp)','All','Shuf','tdT+'};
line_color = [0 1 0; 1 0 1; 0 0 1; 0 0 0; 0.5 0.5 0.5; 1 0 0];

hLine(1) = boundedline(time_vec,mean(lr.cvAccNR),std(lr.cvAccNR,0,1,'omitnan')./sqrt(size(lr.cvAccNR,1)),...
    'nan','gap',...
    'transparency',0.1);
hLine(2) = boundedline(time_vec,mean(lr.cvAccMixedUR),std(lr.cvAccMixedUR,0,1,'omitnan')./sqrt(size(lr.cvAccMixedUR,1)),...
    'nan','gap',...
    'transparency',0.1);
hLine(3) = boundedline(time_vec,mean(lr.cvAccMixedUU),std(lr.cvAccMixedUU,0,1,'omitnan')./sqrt(size(lr.cvAccMixedUU,1)),...
    'nan','gap',...
    'transparency',0.1);
hLine(4) = line(time_vec,lr.cvAcc);
hLine(5) = line(time_vec,lr.cvAccShuf,...
    'linewidth',2);
hLine(6) = line(time_vec,lr.cvAccRed);

for i = 1:length(hLine)
set(hLine(i),'color',line_color(i,:));
text(hLine(i).XData(end)+1,hLine(i).YData(end),line_label{i},...
    'color',line_color(i,:),...
    'fontweight','bold');
end

str_title_label = strjoin(string(horzcat(meta(:,1),meta(:,6),meta(:,3),meta(:,5),meta(:,end))));
title(str_title_label,'FontSize',20)
xlabel('Frame'); ylabel('Classifier accuracy (%)');

ax = gca;
yT = 0.3:0.1:1; yTL = mat2cell(yT*100,size(yT,1),size(yT,2));
set(ax,'ytick',0.3:0.1:1,...
    'yticklabels',yTL);
xline_labels = {S.epoches(1).name,S.epoches(2).name,S.epoches(5).name,S.epoches(6).name};
xline(S.segFrames(1:end-1)+1,'-',xline_labels,...
    'labelverticalalignment','bottom',...
    'labelhorizontalalignment','right',...
    'labelorientation','horizontal',...
    'fontsize',8,...
    'fontweight','bold');

offsetAxes(ax);
fig_configAxis(ax);
exportgraphics(hFig,fullfile(S.dir.imagingRootDir,'LogisticRegression',[meta{:,1} '_' meta{:,6} '.pdf']));
disp(['Figure saved as ' fullfile(S.dir.imagingRootDir,'LogisticRegression',[meta{:,1} '_' meta{:,6} '.pdf'])]);