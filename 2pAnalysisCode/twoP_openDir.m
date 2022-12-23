function twoP_openDir(animal,session)

S = twoP_settings;
open_dir = fullfile(S.dir.imagingRootDir,animal,'imaging',session,S.dir.imagingSubDir);
winopen(open_dir);