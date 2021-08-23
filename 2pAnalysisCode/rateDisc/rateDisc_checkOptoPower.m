function cBhv = rateDisc_checkOptoPower(cBhv, cAnimal)
% check for missing power value of optogenetic stimulation. in some of the
% older recordings, those were not captured properly and will be added
% here.

if length(unique(cBhv.AnimalID)) > 1
    disp('Multiple animals detected. CheckOptoPower is only meant for single animals. No changes made');
else
    oInd = cBhv.optoDur > 0 & isnan(cBhv.optoPower); %optogenetic trials without power value

    if any(ismember({'mSM80' 'mSM81' 'mSM82'}, cAnimal))
        cBhv.optoPower(oInd & cBhv.date < datenum('6/14/19')) = 10; %high power dates
        cBhv.optoPower(oInd & cBhv.date >= datenum('6/14/19') & cBhv.date <= datenum('7/1/19')) = 1; %low power dates
        
    elseif any(ismember({'mSM83' 'mSM84'}, cAnimal))
        cBhv.optoPower(oInd) = 10; %high power dates
        
    elseif ismember({'Fez7'}, cAnimal)
        cBhv.optoPower(oInd & cBhv.date <= datenum('5/6/19')) = 0; %these are control trials before induction
        cBhv.optoPower(oInd & cBhv.date > datenum('5/6/19')) = 10; %rest is high power
        
    elseif any(ismember({'Plex05' 'Plex06'}, cAnimal))
        cBhv.optoPower(oInd) = 10; %these are high power trials
    end
end
        
        
        