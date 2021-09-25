function lightsheet_runZProjectBatch
nImg = 20; zProjectMethod = 'maximum';
baseDir = 'F:\LightsheetBatch2020-10-18';
samples = dir(baseDir);
samples(ismember({samples.name}, {'.', '..'}))=[];
notDir = ~cell2mat({samples.isdir});
samples(notDir)=[];
samples = {samples.name};
for i = 1:length(samples)
    res = lightsheet_averageStack(samples{i}, nImg, zProjectMethod);
end
