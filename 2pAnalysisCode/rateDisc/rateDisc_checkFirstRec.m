function [firstDate, firstFile] = rateDisc_checkFirstRec(bhvPath, animal)
% function to check date of first recording for a given animal. This is to
% know when an animal was in the setup for the first time.

paradigm = 'SpatialDisc';
minTrials = 50; %minimal amount of trials to count as a true session
for iChecks = 1:10 %check for files repeatedly. Sometimes the server needs a moment to be indexed correctly
    Files = ls([bhvPath '\' animal '\' paradigm '\Session Data\*' animal '*']); %behavioral files in correct cPath
    if ~isempty(Files)
        break;
    end
    pause(0.1);
end

cPath = [bhvPath animal '\' paradigm '\Session Data\']; %folder with behavioral data
cDate = datenum(Files(:,length([animal '_' paradigm '_'])+1:length([animal '_' paradigm '_'])+10)); %isolate dates from Filenames
cDate = cDate + (str2num(Files(:,length([animal '_' paradigm '_'])+19:end-4))*0.01); %add session nr to timestamp

[cDate,ind] = sort(cDate,'ascend'); %sort counts to get the order of files to days correct. Oldest file should be first in the list.
cDate = floor(cDate);
Files = Files(ind,:); %adjust order of filenames to get it to be chronological

%check trialcount to find the first 'real' session
for iFiles = 1:size(Files,1)

    load([cPath Files(iFiles,:)], 'SessionData'); %load current bhv file
    if length(SessionData.Rewarded) > minTrials
        break
    end
end

% return name of first behavioral file and its date number
firstFile = Files(iFiles,:);
firstDate = datenum(Files(iFiles,length([animal '_' paradigm '_'])+1:length([animal '_' paradigm '_'])+10)); %isolate dates from Filenames

