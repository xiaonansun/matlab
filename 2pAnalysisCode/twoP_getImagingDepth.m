function colDepthCategory = twoP_getImagingDepth(colDepthMicrons)
% Input: colDepthMicrons is a column of cells imported from the google
% sheets documents were the recording depth was recorded manually
% Output: colDepthCategory is a column of cells with three categories of
% recording depths: Superficial, Intermediate, and Deep
%%
superficial = 200;
deep = 400;

colDepthCategory = cell(length(colDepthMicrons),1);
depth = cellfun(@str2num,colDepthMicrons,'UniformOutput',false);
for i = 1:length(depth)
    if ~isempty(depth{i})
            if depth{i} < superficial
                colDepthCategory{i} = 'Superficial';
            elseif depth{i} > deep
                colDepthCategory{i} = 'Deep';
            elseif depth{i} > superficial && depth{i} < deep
                colDepthCategory{i} = 'Intermediate';
            end
    else
        colDepthCategory{i} = [];
    end
end

        