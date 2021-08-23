function [cTrainDates, cTrainLabels] = rateDisc_convertTrainDates(trainDates)

if ischar(trainDates{1})
    animalCnt = 1;
    trainDates{1} = trainDates;
else
    animalCnt = length(trainDates);
end

for iAnimals = 1 : animalCnt
    
    Cnt = 1;
    a = textscan(trainDates{iAnimals}{2}, '%s%s%s','Delimiter', ' '); %onset of detection learning
    cTrainDates{iAnimals}{Cnt} = a{3}{1};
    cTrainLabels{Cnt} = 'AudioDetect-Learn-5thPrc';
    
    Cnt = Cnt + 1;
    a = textscan(trainDates{iAnimals}{3}, '%s%s%s','Delimiter', ' '); %animal in the middle of detection performance
    cTrainDates{iAnimals}{Cnt} = a{3}{1};
    cTrainLabels{Cnt} = 'AudioDetect-Learn-50thPrc';
    
    Cnt = Cnt + 1;
    a = textscan(trainDates{iAnimals}{4}, '%s%s%s','Delimiter', ' '); %animal close to top detection performance
    cTrainDates{iAnimals}{Cnt} = a{3}{1};
    cTrainLabels{Cnt} = 'AudioDetect-Learn-95thPrc';

    Cnt = Cnt + 1;
    a = textscan(trainDates{iAnimals}{5}, '%s%s%s','Delimiter', ' '); %last detection session
    cTrainDates{iAnimals}{Cnt} = a{3}{1};
    cTrainLabels{Cnt} = 'AudioDetect-Last';

    Cnt = Cnt + 1;
    a = textscan(trainDates{iAnimals}{7}, '%s%s%s%s','Delimiter', ' '); %first discrimination session
    cTrainDates{iAnimals}{Cnt} = a{4}{1};    
    cTrainLabels{Cnt} = 'AudioDisc-first';

    Cnt = Cnt + 1;
    a = textscan(trainDates{iAnimals}{8}, '%s%s%s%s','Delimiter', ' '); %last discrimination session
    cTrainDates{iAnimals}{Cnt} = a{4}{1};
    cTrainLabels{Cnt} = 'AudioDisc-last';

    if ~isempty(trainDates{iAnimals}{9})
    Cnt = Cnt + 1;
    a = textscan(trainDates{iAnimals}{9}, '%s%s%s%s%s','Delimiter', ' '); %last discrimination session
    cTrainDates{iAnimals}{Cnt} = a{5}{1};
    cTrainLabels{Cnt} = 'TacDisc-Novice-first';
    
    Cnt = Cnt + 1;
    a = textscan(trainDates{iAnimals}{10}, '%s%s%s%s%s','Delimiter', ' '); %last discrimination session
    cTrainDates{iAnimals}{Cnt} = a{5}{1};
    cTrainLabels{Cnt} = 'TacDisc-Novice-last';

    Cnt = Cnt + 1;
    a = textscan(trainDates{iAnimals}{12}, '%s%s%s','Delimiter', ' '); %%onset of tactile detection learning
    cTrainDates{iAnimals}{Cnt} = a{3}{1};
    cTrainLabels{Cnt} = 'TacDetect-Learn-5thPrc';
    
    Cnt = Cnt + 1;
    a = textscan(trainDates{iAnimals}{13}, '%s%s%s','Delimiter', ' '); %%onset of tactile detection learning
    cTrainDates{iAnimals}{Cnt} = a{3}{1};
    cTrainLabels{Cnt} = 'TacDetect-Learn-50thPrc';

    Cnt = Cnt + 1;
    a = textscan(trainDates{iAnimals}{14}, '%s%s%s','Delimiter', ' '); %animal close to top tactile detection performance
    cTrainDates{iAnimals}{Cnt} = a{3}{1};
    cTrainLabels{Cnt} = 'TacDetect-Learn-95thPrc';

    Cnt = Cnt + 1;
    a = textscan(trainDates{iAnimals}{16}, '%s%s%s%s','Delimiter', ' '); %first tactile discrimination session
    cTrainDates{iAnimals}{Cnt} = a{4}{1};
    cTrainLabels{Cnt} = 'TacDisc-first';

    Cnt = Cnt + 1;
    a = textscan(trainDates{iAnimals}{17}, '%s%s%s%s','Delimiter', ' '); %last tactile discrimination session
    cTrainDates{iAnimals}{Cnt} = a{4}{1};
    cTrainLabels{Cnt} = 'TacDisc-last';
    end
end

if animalCnt == 1
    cTrainDates = cTrainDates{1};
end
    