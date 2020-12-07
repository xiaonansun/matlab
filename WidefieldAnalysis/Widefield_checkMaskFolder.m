function Widefield_checkMaskFolder(cPath)

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

%% check for subfolders
allRecs = dir(cPath);
allRecs(1:2) = [];

for iRecs = 1:length(allRecs)
    try
        Widefield_checkMask([cPath allRecs(iRecs).name])
        disp(['Completed folder: ' cPath allRecs(iRecs).name]);
    catch
        disp(['Error in folder: ' cPath allRecs(iRecs).name]);
    end
end
    
%%    

