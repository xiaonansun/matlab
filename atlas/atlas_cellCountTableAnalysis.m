% https://docs.google.com/spreadsheets/d/16IUkIXvKGEBwxvWFq-6DR1IjNmWZDRVS/edit?usp=sharing&ouid=115298005536272348157&rtpof=true&sd=true

animal_docid = GetGoogleSpreadsheet('1P04X7NLNY71Yg9aITjBGWQpAAIOm0uAYm4yyF9_gfBk');

idx_animal = 1;
disp(['Loading ' animal_docid{idx_animal,1}])
docid = animal_docid{idx_animal,2};
Z = GetGoogleSpreadsheet(docid); Z = Z(2:end,:);

docid_dorsal_areas= '1UNq5yuJfFgn2O8qYvjuxxOYtSKECb6Go9EfDTMWt9bw';
dorsal_areas = GetGoogleSpreadsheet(docid_dorsal_areas);
rowColLabels = 1;
strAllenStrctID = 'Allen_structure_ID';
colAllenStructID = find(contains(dorsal_areas(rowColLabels,:),strAllenStrctID));
AllenStructID = cellfun(@str2double,dorsal_areas(rowColLabels+1:end,colAllenStructID));

%%
strStructIDPath = 'structure_id_path';
strStructID = 'structure ID';
strFullStructName = 'full structure name';
strTotalVoxelCounts = 'total_voxel_counts';
strIpsiCount = 'cell count ipsi';
strContraCount = 'cell count contra';
strBilatCount = 'cell count all';

colStructIDPath = find(contains(Z(rowColLabels,:),strStructIDPath));
colStructID = find(contains(Z(rowColLabels,:),strStructID));
colFullStructName = find(contains(Z(rowColLabels,:),strFullStructName));
colTotalVoxelCounts = find(contains(Z(rowColLabels,:),strTotalVoxelCounts));
colIpsiCount = find(contains(Z(rowColLabels,:),strIpsiCount));
colContraCount = find(contains(Z(rowColLabels,:),strContraCount));
colBilatCount = find(contains(Z(rowColLabels,:),strBilatCount));

full_struct_name = Z(2:end,colFullStructName);
struct_ID_all = Z(2:end,colStructID); struct_ID_all = cellfun(@str2double,struct_ID_all);
ipsi_count = cellfun(@str2double,Z(2:end,colIpsiCount)); % extracts cell count from csv table sent by Katie Matho
contra_count = cellfun(@str2double,Z(2:end,colContraCount));
bilat_count = cellfun(@str2double,Z(2:end,colBilatCount));
voxel_count = cellfun(@str2double,Z(2:end,colTotalVoxelCounts));

vecStructIDPath = cellfun(@(x) regexp(x,'/','split'),Z(2:end,colStructIDPath),'UniformOutput',false); % vectorizes the Allen structure ID path
vecStructIDPath = cellfun(@str2double,vecStructIDPath,'UniformOutput',false);
vecStructIDPath = cellfun(@(x) x(~isnan(x)),vecStructIDPath,'UniformOutput',false);

idx_layer23 = contains(cellfun(@lower,full_struct_name,'UniformOutput',false), 'layer 2/3'); 
idx_layer5 = contains(cellfun(@lower,full_struct_name,'UniformOutput',false), 'layer 5');
idx_layer6a = contains(cellfun(@lower,full_struct_name,'UniformOutput',false), 'layer 6a');
idx_layer6b = contains(cellfun(@lower,full_struct_name,'UniformOutput',false), 'layer 6b');

%% Compute the per-area, per-layer cell count or average fluorescence

area_and_layers = cell(length(AllenStructID),1);
area_struct_ids = cell(length(AllenStructID),1);
dorsal_area_names = dorsal_areas(2:end,3);
count.area_bilat_all_layer = zeros(length(AllenStructID),1);
count.area_bilat_layer23 = zeros(length(AllenStructID),1);
count.area_bilat_layer5 = zeros(length(AllenStructID),1);
count.area_bilat_layer6a6b = zeros(length(AllenStructID),1);

voxel.area_bilat_all_layer = zeros(length(AllenStructID),1);
voxel.area_bilat_layer23 = zeros(length(AllenStructID),1);
voxel.area_bilat_layer5 = zeros(length(AllenStructID),1);
voxel.area_bilat_layer6a6b = zeros(length(AllenStructID),1);

area_bilat_all_layer_F = zeros(length(AllenStructID),2); % column 1 is right, column 2 is left
area_bilat_layer23_F = zeros(length(AllenStructID),2);
area_bilat_layer5_F = zeros(length(AllenStructID),2);
area_bilat_layer6a6b_F = zeros(length(AllenStructID),2);

for i = 1:length(AllenStructID)
    tic
    idx_area_all_layers = cellfun(@(x) find(x == AllenStructID(i)),vecStructIDPath,'UniformOutput',false);
    idx_area_all_layers = logical(~cellfun(@isempty,idx_area_all_layers));
    area_and_layers{i} = full_struct_name(idx_area_all_layers);
    disp(['Working on ' area_and_layers{i}{1} '.'])
    area_struct_ids{i} = struct_ID_all(idx_area_all_layers);
    
    count.area_bilat_all_layer(i) = sum(bilat_count(idx_area_all_layers));
    count.area_bilat_layer23(i) = sum(bilat_count(idx_area_all_layers & idx_layer23));
    count.area_bilat_layer5(i) = sum(bilat_count(idx_area_all_layers & idx_layer5));
    count.area_bilat_layer6a6b(i) = sum(bilat_count(idx_area_all_layers & (idx_layer6a | idx_layer6b)));
    
    voxel.area_bilat_all_layer(i) = sum(voxel_count(idx_area_all_layers));
    voxel.area_bilat_layer23(i) = sum(voxel_count(idx_area_all_layers & idx_layer23));
    voxel.area_bilat_layer5(i) = sum(voxel_count(idx_area_all_layers & idx_layer5));
    voxel.area_bilat_layer6a6b(i) = sum(voxel_count(idx_area_all_layers & (idx_layer6a | idx_layer6b)));

%     area_bilat_all_layer_F(i,1) = mean(sample(ismember(atlasRight,struct_ID_all(idx_area_all_layers))));
%     area_bilat_all_layer_F(i,2) = mean(sample(ismember(atlasLeft,struct_ID_all(idx_area_all_layers))));
%     area_bilat_layer23_F(i,1) = mean(sample(ismember(atlasRight,struct_ID_all(idx_area_all_layers & idx_layer23))));
%     area_bilat_layer23_F(i,2) = mean(sample(ismember(atlasLeft,struct_ID_all(idx_area_all_layers & idx_layer23))));
%     area_bilat_layer5_F(i,1) = mean(sample(ismember(atlasRight,struct_ID_all(idx_area_all_layers & idx_layer5))));
%     area_bilat_layer5_F(i,2) = mean(sample(ismember(atlasLeft,struct_ID_all(idx_area_all_layers & idx_layer5))));
%     area_bilat_layer6a6b_F(i,1) = mean(sample(ismember(atlasRight,struct_ID_all(idx_area_all_layers & (idx_layer6a | idx_layer6b)))));
%     area_bilat_layer6a6b_F(i,2) = mean(sample(ismember(atlasLeft,struct_ID_all(idx_area_all_layers & (idx_layer6a | idx_layer6b)))));

    disp(['Completed ' area_and_layers{i}{1} ' in ' num2str(toc) ' seconds.']);
end

    F.area_bilat_all_layer = mean(area_bilat_all_layer_F,2);
    F.area_bilat_layer23 = mean(area_bilat_layer23_F,2);
    F.area_bilat_layer5 = mean(area_bilat_layer5_F,2);
    F.area_bilat_layer6a6b = mean(area_bilat_layer6a6b_F,2);
%% Save data files

if ispc
    save_dir = 'E:\My Drive\_manuscript_Simon_Anne_cell_type\cell_count_analysis\cell_count_output';
elseif ismac
    save_dir = '/Volumes/GoogleDrive/My Drive/_manuscript_Simon_Anne_cell_type/cell_count_analysis/cell_count_output';
end

save(fullfile(save_dir,[animal_docid{idx_animal,1} '.mat']),'dorsal_area_names', 'area_and_layers',...
    'count',...
    'voxel');

% save(fullfile(save_dir,[dataset_name '.mat']), 'dorsal_area_names', 'area_and_layers','F');

%%
% layer_label = {'layer 1',...
%     'layer 2/3',...
%     'layer 4',...
%     'layer 5',...
%     'layer 6a',...
%     'layer 6b'};

