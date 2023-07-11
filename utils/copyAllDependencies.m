function copyAllDependencies(inFunName)
% This function copies all dependencies into the same directory as the
% input .m file.
% To best use this function, it would be ideal for the .m file to reside in
% its own directory (e.g. a subfolder within the github folder)

%%
inFunName = 'C:\Users\xiaonansun\Documents\GitHub\FUS-imaging\fus_importSonicationImages.m';
[outDir,~,~] = fileparts(inFunName);
[fList,pList] = matlab.codetools.requiredFilesAndProducts(inFunName);

for i = 1:length(fList)
    %% i = 1
    if ~strcmp(inFunName,fList{i})
        copyfile(fList{i},outDir);
        disp([fList{i} ' copied to ' outDir]);
    end
end
