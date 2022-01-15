exps = twoP_getAcquisitionRecord;

colAnimal = 1;
colSession = 6;
colLocation = 3;
errors = cell(size(exps,1),1);

parfor i = 2:size(exps,1)
    animal = exps{i,1};
    session = exps{i,6};
    try
        twoP_plotSingleSessionIS(animal,session)
    catch ME
        errors{i} = ME;
    end
    
    disp([animal ' ' session ' plotted!']);
end