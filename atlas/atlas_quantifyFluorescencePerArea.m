%% 
docid = '1eiv9ApIrbnldjnpThiz4CENReczkeWomky-xpTgWHV0';
Z = GetGoogleSpreadsheet(docid); C = Z(2:end,:);
table_path = "C:\Users\Xiaonan Richard Sun\.brainglobe\allen_mouse_10um_v1.2\structures.csv";

cell_type = 'Fezf2';
% cell_type = 'PlexinD1';
tiffDir = 'I:\brainreg-register_output\fez\2022-05-04_1';
hemiFileName = 'registered_hemispheres.tiff';
atlasFileName = 'registered_atlas.tiff';
sampleFileName = 'downsampled.tiff';

% T = readtable(table_path);
% C = table2cell(T);
acronym = C(:,1);ec
id = C(:,2);

areas = {'SS','MO','VIS','AUD','RSP','ACA'};
% layers = {'2/3','4','5'};
layers = {'4','5','6'};
exclude = {'1','2','3','4','5','6','7','8','9','0'};

%%
for i = 1:length(areas)
    idxAreas(:,i) = contains(acronym,areas{i});
end


for i = 1:length(layers)
    idxLayers(:,i) = contains(acronym,layers{i});
end


for i = 1:length(exclude)
    idxExclude(:,i) = contains(acronym,exclude{i});
end

%%
idxSelected = sum(idxAreas,2) & sum(idxLayers,2);
idxSelectedName = acronym(find(diff(idxSelected) == 1)-2);
idxSup = find(diff(idxSelected) == 1)+1;
idxDeep = find(diff(idxSelected) == -1);

cnt = 1;
for i = 1:length(idxSelectedName)
    idxAllLayers{i} = (cnt:1:(cnt+idxDeep(i)-idxSup(i)))';
    cnt = cnt+1+idxDeep(i)-idxSup(i);
end
% selectedID = uint32(zeros(size(idxSelected,1),size(idxSelected,2)));
% z = id(idxSelected);
selectedID = cellfun(@str2num,id(idxSelected));

% idxSelectedName = sum(idxAreas,2) & ~sum(idxExclude,2);

% plot(1:1:length(idxSelected),idxSelected,'.k');
% plot(1:1:length(idxSelectedName),idxSelectedName,'.k');

%% 
hemisphere = loadtiff(fullfile(tiffDir,hemiFileName)); % 1 = right, 2 = left
hemisphere = hemisphere(:);
atlas = loadtiff(fullfile(tiffDir,atlasFileName)); atlas = atlas(:);
atlasRight = atlas.*(uint32(hemisphere == 1));
atlasLeft = atlas.*(uint32(hemisphere == 2));

sample = loadtiff(fullfile(tiffDir,sampleFileName)); sample = sample(:);
% idAtlas = ismember(atlas,selectedID);

for i = 1:length(selectedID)
    idxAtlasRight = find(atlasRight == selectedID(i));
    idxAtlasLeft = find(atlasLeft == selectedID(i));
    sampleArea{1,i} = sample(idxAtlasRight);
    sampleArea{2,i} = sample(idxAtlasLeft);
end

%%
meanF = cellfun(@mean,sampleArea)';
for i = 1:length(idxAllLayers)
meanF_allLayers(i,:) = sum(meanF(idxAllLayers{i},:));

% meanF_L23(i,:) = meanF(idxAllLayers{i}(1),:); % use for plexin

% meanF_L5(i,:) = meanF(idxAllLayers{i}(end),:); % use for plexin
meanF_L5(i,:) = meanF(idxAllLayers{i}(end-2),:);
meanF_L6(i,:) = mean(meanF(idxAllLayers{i}(end-1:end),:));
end

save_dir = 'C:\Users\Xiaonan Richard Sun\Dropbox\Users\Richard\whole_brain_expression_data';
save_file_name = [cell_type '_' datestr(datetime,'yyyy-mm-dd_HH-MM-SS')];
% save(fullfile(save_dir,cell_type,save_file_name),'idxSelectedName','meanF_allLayers','meanF_L23','meanF_L5');
save(fullfile(save_dir,save_file_name),'idxSelectedName','meanF_allLayers','meanF_L6','meanF_L5');