function [hFigure, hTiles] = twoP_histEpoches(input_matrix,input_index,num_bins)
% This function plots the histogram of each column, up to 10 columns
%
% INPUTS
% input_matrix: rows represent values of cells, columnes
% represent epoches, NaNs are allowed
% input_index: optional input variable.

num_max_epoches = 10;
inMat = input_matrix;

if (size(inMat,2) > size(inMat,1)) && (size(inMat,2) <= num_max_epoches || size(inMat,1) <= num_max_epoches)
    inMat = inMat';
elseif (size(inMat,2) >= num_max_epoches) && (size(inMat,1) >= num_max_epoches)
    disp(['Number of plotted epoches must be small than 10. Please check the input size.']);
    return
end

if ~exist('input_index','var') || isempty(input_index)
    inIndex = ones(length(inMat),1,'logical');
else
    inIndex = input_index;
end

if ~exist('num_bins','var') || isempty(num_bins)
    num_bins = 20;
end

limX = [0 1];
hFigure = figure('Position',[500 500 1200 250]);
hTiles = tiledlayout(1,size(inMat,2));
edgecolor = {'k', 'r', 'g'};
histlegends = {'\color{black} All','\color{red} tdT+','\color{green} tdT-'};

for j = 1:size(inMat,2)
        nexttile
        hold on
        for jj = 1:size(inIndex,2)
            hH = histogram(inMat(~isnan(inMat(inIndex(:,jj),j)),j), num_bins,...
                'BinLimits',limX,...
                'DisplayStyle','stairs',...
                'EdgeColor',edgecolor{jj},...
                'FaceColor','none',...
                'LineWidth',1,...
                'Normalization','probability');
        end
        xlabel('AUC');
        if j == 1
           ylabel('Fraction of cells (%)');
        end
        if j == size(inMat,2)
        hL = legend(histlegends,...
            'Location','northeast',...
            'Box','off');
        end
        offsetAxes(gca);
        fig_configAxis(gca);
end


%     title(hTiles,[cell2mat(sAnimal') ' ' sExpertise{:} ' ' sLocation{:} ' ' sDepth{:} ' depth']);
%     exportgraphics(hFigure,fullfile(S.dir.imagingRootDir,'AUC',[cell2mat(sAnimal') '_' sExpertise{:} '_' sLocation{:} '_' sDepth{:} '_AUC_combined.pdf']));