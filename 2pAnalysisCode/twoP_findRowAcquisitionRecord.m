function idx = twoP_findRowAcquisitionRecord(animal,session)

exps = twoP_getAcquisitionRecord;

colAnimal = contains(exps(1,:),'Animal'); colSession = contains(exps(1,:),'Folder');
idxAnimal = strcmp(exps(:,colAnimal),animal);
idxSession = strcmp(exps(:,colSession),session);

idx = find(idxAnimal.*idxSession);